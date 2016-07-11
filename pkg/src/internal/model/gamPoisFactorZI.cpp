// Copyright 2016-05 Ghislain Durif
//
// This file is part of the `countMatrixFactor' library for R and related languages.
// It is made available under the terms of the GNU General Public
// License, version 2, or at your option, any later version,
// incorporated herein by reference.
//
// This program is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public
// License along with this program; if not, write to the Free
// Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
// MA 02111-1307, USA

/*!
* \file gamPoisFactorZI.cpp
* \brief class definition for penalized Gamma Poisson Factor Model
* \author Ghislain Durif
* \version 0.1
* \date 09/05/2016
*/

#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <cstdio>
#include <vector>
#include "gamPoisFactorZI.h"
#include "gamDistrib.h"
#include "intermediate.h"

#define mdirac() unaryExpr(std::ptr_fun<double,double>(intermediate::dirac))
#define mexp() unaryExpr(std::ptr_fun<double,double>(std::exp))
#define mlgamma() unaryExpr(std::ptr_fun<double,double>(lgamma))
#define mlog() unaryExpr(std::ptr_fun<double,double>(std::log))
#define mlogit() unaryExpr(std::ptr_fun<double,double>(intermediate::logit))
#define mlogitinv() unaryExpr(std::ptr_fun<double,double>(intermediate::logitinv))

// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision

using std::vector;

// SRC
namespace countMatrixFactor {

    // CONSTRUCTOR
    gamPoisFactorZI::gamPoisFactorZI(int n, int p, int K, const MatrixXi &X,
                                       const MatrixXd &phi1, const MatrixXd &phi2,
                                       const MatrixXd &theta1, const MatrixXd &theta2,
                                       const MatrixXd &alpha1, const MatrixXd &alpha2,
                                       const MatrixXd &beta1, const MatrixXd &beta2)
    : gamPoisFactor::gamPoisFactor(n, p, K, X, phi1, phi2, theta1, theta2,
      alpha1, alpha2, beta1, beta2)
    {
        m_ZI = true;

        // ZI probabilities and frequencies
        m_prob0 = VectorXd::Zero(p);
        m_prob = VectorXd::Zero(p);
        m_freq = VectorXd::Zero(p);

        // sufficient stats
        m_EZ_logU_k = MatrixXd::Zero(n,p);
        m_EZ_logV_k = MatrixXd::Zero(n,p);
        m_EU_EV_k = MatrixXd::Zero(n,p);
        m_ElogU_ElogV_k = MatrixXd::Zero(n,p);
        m_ElgamZ_k = MatrixXd::Zero(n,p);

        // Rcpp::Rcout << "dim of m_EZ_logU_k = " << m_EZ_logU_k.rows() << "x" << m_EZ_logU_k.cols() << std::endl;
        // Rcpp::Rcout << "dim of m_EZ_logV_k = " << m_EZ_logV_k.rows() << "x" << m_EZ_logV_k.cols() << std::endl;
        // Rcpp::Rcout << "dim of m_ElgamZ_k = " << m_ElgamZ_k.rows() << "x" << m_ElgamZ_k.cols() << std::endl;
        //
        // Rcpp::Rcout << "m_EZ_logU_k = " << m_EZ_logU_k << std::endl << std::endl;
        // Rcpp::Rcout << "m_EZ_logV_k = " << m_EZ_logV_k << std::endl << std::endl;
        // Rcpp::Rcout << "m_ElgamZ_k = " << m_ElgamZ_k << std::endl << std::endl;
    }

    // DESTRUCTOR
    gamPoisFactorZI::~gamPoisFactorZI() {}

    // member functions: documented in src

    /*!
    * \brief Initialization of sufficient statistics
    */
    void gamPoisFactorZI::Init() {

        // Gamma variational parameter
        Rcpp::Rcout << "Init: Gamma variational parameter" << std::endl;
        for(int k=0; k<m_K; k++) {
            // local parameters
            for(int i=0; i<m_N; i++) {
                double param1 = 0;
                double param2 = 0;
                estimParam(1000, m_alpha1cur(i,k), m_alpha2cur(i,k), param1, param2);
                m_phi1cur(i,k) = param1;
                m_phi1old(i,k) = param1;
                m_phi2cur(i,k) = param2;
                m_phi2old(i,k) = param2;
            }
            // global parameters
            for(int j=0; j<m_P; j++) {
                double param1 = 0;
                double param2 = 0;
                estimParam(1000, m_beta1cur(j,k), m_beta2cur(j,k), param1, param2);
                m_theta1cur(j,k) = param1;
                m_theta1old(j,k) = param1;
                m_theta2cur(j,k) = param2;
                m_theta2old(j,k) = param2;
            }
        }

        // Bernoulli proba
        for(int j=0; j<m_P; j++) {
            for(int i=0; i<m_N; i++) {
                if(m_X(i,j) != 0) {
                    m_freq(j) = m_freq(j)+1;
                }
            }
        }

        m_freq = m_freq.array() / m_N;
        m_prob0 = m_freq;
        m_prob = m_freq;

        // sufficient statistics
        Rcpp::Rcout << "Init: sufficient statistics" << std::endl;
        Egam(m_phi1cur, m_phi2cur, m_EU);
        Elgam(m_phi1cur, m_phi2cur, m_ElogU);
        Egam(m_theta1cur, m_theta2cur, m_EV);
        Elgam(m_theta1cur, m_theta2cur, m_ElogV);

        // Multinomial Z parameters
        this->multinomParam();

    }

    //-------------------//
    // parameter updates //
    //-------------------//

    /*!
     * \brief update rule for multinomial parameters in variational inference
     *
     * m_EZ_i_{jk} = sum_i pi_j * E[Z_{ijk}]
     * m_EZ_j_{ik} = sum_j pi_j * E[Z_{ijk}]
     */
    void gamPoisFactorZI::multinomParam() {

        intermediate::checkExp(m_ElogU);
        intermediate::checkExp(m_ElogV);

        // sum_k exp(E[log(U_{ik})]) * exp(E[log(V_{jk})])
        m_ElogU_ElogV_k = m_ElogU.mexp() * m_ElogV.mexp().transpose();

        m_EZ_j = m_ElogU.mexp().array() * ( ((m_X.cast<double>().array().rowwise() * m_prob.transpose().array()) / m_ElogU_ElogV_k.array() ).matrix() * m_ElogV.mexp() ).array();
        m_EZ_i = m_ElogV.mexp().array() * ( ((m_X.cast<double>().array().rowwise() * m_prob.transpose().array()) / m_ElogU_ElogV_k.array() ).matrix().transpose() * m_ElogU.mexp() ).array();

        // test
        // for(int i=0; i<m_N; i++) {
        //     for(int k = 0; k<m_K; k++) {
        //         double test = 0;
        //         for(int j=0; j<m_P; j++) {
        //             test += m_prob(j) * m_X(i,j) * std::exp(m_ElogV(j,k)) / m_ElogU_ElogV_k(i,j);
        //         }
        //         test *= std::exp(m_ElogU(i,k));
        //
        //         if(test != m_EZ_j(i,k)) {
        //             Rcpp::Rcout << "test = " << test << " m_EZ_j = " <<  m_EZ_j(i,k) << std::endl;
        //         }
        //     }
        // }
        //
        // for(int j=0; j<m_P; j++) {
        //     for(int k = 0; k<m_K; k++) {
        //         double test = 0;
        //         for(int i=0; i<m_N; i++) {
        //             test += m_prob(j) * m_X(i,j) * std::exp(m_ElogU(i,k)) / m_ElogU_ElogV_k(i,j);
        //         }
        //         test *= std::exp(m_ElogV(j,k));
        //         if(test != m_EZ_i(j,k)) {
        //             Rcpp::Rcout << "test = " << test << " m_EZ_i = " <<  m_EZ_i(j,k) << std::endl;
        //         }
        //     }
        // }
    }

    /*!
     * \brief update rule for local parameters phi (factor U) in variational inference
     */
    void gamPoisFactorZI::localParam() {
        m_phi1cur = m_alpha1cur.array() + m_EZ_j.array();
        m_phi2cur = m_alpha2cur.array().rowwise() + (m_EV.array().colwise() * m_prob.array()).colwise().sum();

        // test
        // for(int i=0; i<m_N; i++) {
        //    for(int k = 0; k<m_K; k++) {
        //        double test = m_alpha2cur(i,k) + m_prob.dot(m_EV.col(k));
        //        if(test != m_phi2cur(i,k)) {
        //            Rcpp::Rcout << "test = " << test << " m_phi2cur = " <<  m_phi2cur(i,k) << std::endl;
        //        }
        //    }
        // }

        // expectation and log-expectation
        Egam(m_phi1cur, m_phi2cur, m_EU);
        Elgam(m_phi1cur, m_phi2cur, m_ElogU);
    }

    /*!
     * \brief update rule for global parameters theta (factor V) in variational inference
     */
    void gamPoisFactorZI::globalParam() {
        m_theta1cur = m_beta1cur.array() + m_EZ_i.array();
        m_theta2cur = ((MatrixXd::Zero(m_N, m_K).rowwise() + m_EU.colwise().sum()).array().colwise() * m_prob.array()) + m_beta2cur.array();

        // test
        // for(int j=0; j<m_P; j++) {
        //     for(int k = 0; k<m_K; k++) {
        //         double test = m_beta2cur(j,k) + m_prob(j) * m_EU.col(k).sum();
        //         if(test != m_theta2cur(j,k)) {
        //             Rcpp::Rcout << "test = " << test << " m_theta2cur = " <<  m_theta2cur(j,k) << std::endl;
        //         }
        //     }
        // }

        // expectation and log-expectation
        Egam(m_theta1cur, m_theta2cur, m_EV);
        Elgam(m_theta1cur, m_theta2cur, m_ElogV);
    }

    /*!
     * \brief update rule for Bernoulli parameter (of ZI indicator) in variational inference
     */
    void gamPoisFactorZI::ZIproba() {

        // sum_k E[U_{ik}] * E[V_{jk}]
        // Rcpp::Rcout << "sum_k E[U_{ik}] * E[V_{jk}]" << std::endl;
        m_EU_EV_k = m_EU * m_EV.transpose();

        // multinomial probabilities
        VectorXd omega = VectorXd::Zero(m_K);
        VectorXd logGamZ  = VectorXd::Zero(m_K);

        // Rcpp::Rcout << "m_ElogU = " <<  m_ElogU << std::endl;
        // Rcpp::Rcout << "m_ElogV = " <<  m_ElogV << std::endl;
        // Rcpp::Rcout << "m_X = " <<  m_X << std::endl;
        // Rcpp::Rcout << "m_N = " <<  m_N << std::endl;
        // Rcpp::Rcout << "m_P = " <<  m_P << std::endl;
        // Rcpp::Rcout << "m_K = " <<  m_K << std::endl;

        // for(int i = 0; i < m_N; i++) {
        //     for(int j=0; j < m_P; j++) {
        //         // Rcpp::Rcout << "i = " <<  i << std::endl;
        //         // Rcpp::Rcout << "j = " <<  j << std::endl;
        //         // Rcpp::Rcout << "m_X(i,j) = " << m_X(i,j) << std::endl;
        //         // omega(i,j,k) = exp(E[log U_{ik}] + E[log V_{jk}]) / sum_l exp(E[log U_{il}] + E[log V_{jl}])
        //         // Rcpp::Rcout << "omega(i,j,k)" << std::endl;
        //         for(int k=0; k < m_K; k++) {
        //             omega(k) = std::exp(m_ElogU(i,k) + m_ElogV(j,k)) / m_ElogU_ElogV_k(i,j);
        //             logGamZ(k) = intermediate::lgamBinom(m_X(i,j), omega(k));
        //         }
        //         // Rcpp::Rcout << "omega(i,j,.) = " << omega.transpose() << std::endl;
        //         // Rcpp::Rcout << "logGamZ(i,j,.) = " << logGamZ.transpose() << std::endl;
        //         // sum_k E[Z_{ijk}] * E[log U_{ik}]
        //         // Rcpp::Rcout << "sum_k E[Z_{ijk}] * E[log U_{ik}]" << std::endl;
        //         m_EZ_logU_k(i,j) = (double) m_X(i,j) * (m_ElogU.row(i).array() * omega.array()).sum();
        //         // sum_k E[Z_{ijk}] * E[log V_{jk}]
        //         // Rcpp::Rcout << "sum_k E[Z_{ijk}] * E[log V_{jk}]" << std::endl;
        //         m_EZ_logV_k(i,j) = (double) m_X(i,j) * (m_ElogV.row(j).array() * omega.array()).sum();
        //         // sum_k E[log(Z_{ijk}!)]
        //         // Rcpp::Rcout << "sum_k E[log(Z_{ijk}!)]" << std::endl;
        //         m_ElgamZ_k(i,j) = logGamZ.sum();
        //     }
        // }

        Rcpp::Rcout << "m_prob" << std::endl;
        for(int j=0; j<m_P; j++) {
            if((m_prob0(j) < 0.9) && (m_prob0(j) > 0)) {
                // Rcpp::Rcout << "sum_K E[U_{ik}] E[V_{jk}] = " << m_EU_EV_k.col(j).transpose() << std::endl;
                // double Epois = m_EZ_logU_k.col(j).sum() + m_EZ_logV_k.col(j).sum() - m_EU_EV_k.col(j).sum() - m_ElgamZ_k.col(j).sum();
                // double tmp = std::log(m_prob0(j)) - ((double) 1 / (double) m_N) * Epois - m_freq(j) * std::log(1-m_prob0(j));
                // // Rcpp::Rcout << "Epois = " << ((double) 1 / (double) m_N) * Epois << std::endl;
                // // Rcpp::Rcout << "log(m_prob0(j)) = " << std::log(m_prob0(j)) << std::endl;
                // // Rcpp::Rcout << "m_freq(j) = " << m_freq(j) << std::endl;
                // // Rcpp::Rcout << "log(1-m_prob0(j)) = " << std::log(1-m_prob0(j)) << std::endl;
                // // Rcpp::Rcout << "tmp = " << tmp << std::endl;
                // m_prob(j) = intermediate::logitinv(tmp);

                // m_prob(j) = m_freq(j) / (1 - std::exp((-1) * (m_EU_EV_k.col(j).sum())));

                double tmp = intermediate::logit(m_prob0(j)) - intermediate::logit(std::exp((-1) * (m_EU_EV_k.col(j).sum())));
                m_prob(j) = intermediate::logitinv(tmp);
            }
        }
        Rcpp::Rcout << " m_prob = " <<  m_prob.transpose() << std::endl;
    }

    /*!
     * \brief update rule for Bernoulli parameter (of ZI indicator) in prior
     */
    void gamPoisFactorZI::priorZIproba() {
        m_prob0 = m_prob;
    }

    /*!
     * \brief parameter update in standard variational
     */
    void gamPoisFactorZI::updateVarational() {

        // Multinomial parameters
        Rcpp::Rcout << "algorithm: Multinomial parameters" << std::endl;
        this->multinomParam();

        // ZI proba
        Rcpp::Rcout << "algorithm: ZI proba" << std::endl;
        this->ZIproba();

        // local parameters
        // U : param phi
        Rcpp::Rcout << "algorithm: local parameters" << std::endl;
        this->localParam();

        // global parameters
        // V : param theta
        Rcpp::Rcout << "algorithm: global parameters" << std::endl;
        this->globalParam();

        // Poisson rate
        Rcpp::Rcout << "algorithm: Poisson rate" << std::endl;
        this->poissonRate();
    }

    /*!
     * \brief parameter update in variational EM (E-step)
     */
    void gamPoisFactorZI::updateEstep() {
        // Multinomial parameters
        // Rcpp::Rcout << "algorithm: Multinomial parameters" << std::endl;
        this->multinomParam();

        // ZI proba
        //Rcpp::Rcout << "algorithm: ZI proba" << std::endl;
        this->ZIproba();

        // local parameters
        // U : param phi
        // Rcpp::Rcout << "algorithm: local parameters" << std::endl;
        this->localParam();

        // global parameters
        // V : param theta
        // Rcpp::Rcout << "algorithm: global parameters" << std::endl;
        this->globalParam();

        // Poisson rate
        // Rcpp::Rcout << "algorithm: Poisson rate" << std::endl;
        this->poissonRate();
    }

    /*!
     * \brief parameter update in variational EM (M-step)
     */
    void gamPoisFactorZI::updateMstep() {
        // ZI proba
        //Rcpp::Rcout << "algorithm: ZI proba prior" << std::endl;
        this->priorZIproba();

        // local parameters
        // U : param phi
        // Rcpp::Rcout << "algorithm: local parameters" << std::endl;
        this->localPriorParam();

        // global parameters
        // V : param theta
        // Rcpp::Rcout << "algorithm: global parameters" << std::endl;
        this->globalPriorParam();
    }

    //-------------------//
    //       return      //
    //-------------------//

    /*!
    * \brief create list with results to be return
    *
    * @param[out] list containing output
    */
    void gamPoisFactorZI::returnObject(Rcpp::List &results) {

        Rcpp::List returnObj1;
        gamPoisFactor::returnObject(returnObj1);

        Rcpp::List ZI = Rcpp::List::create(Rcpp::Named("prob") = m_prob,
                                           Rcpp::Named("freq") = m_freq,
                                           Rcpp::Named("prob0") = m_prob0);

        Rcpp::List returnObj2 = Rcpp::List::create(Rcpp::Named("ZIparams") = ZI);

        SEXP tmp1 = Rcpp::Language("c", returnObj1, returnObj2).eval();

        SEXP tmp2 = Rcpp::Language("c", results, tmp1).eval();

        results = tmp2;

    }
}