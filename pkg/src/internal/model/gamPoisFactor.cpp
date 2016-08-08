// Copyright 2016-04 Ghislain Durif
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
 * \file gamPoisFactor.cpp
 * \brief class definition for standard Gamma Poisson Factor Model
 * \author Ghislain Durif
 * \version 0.2
 * \date 04/05/2016
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <cstdio>
#include <vector>
#include "explainedVariance.h"
#include "gamPoisFactor.h"
#include "gamDistrib.h"
#include "intermediate.h"
#include "loglikelihood.h"

#define mexp() unaryExpr(std::ptr_fun<double,double>(std::exp))
#define mlgamma() unaryExpr(std::ptr_fun<double,double>(lgamma))
#define mlog() unaryExpr(std::ptr_fun<double,double>(std::log))
#define mpsiInv() unaryExpr(std::bind2nd(std::pointer_to_binary_function<double,int,double>(intermediate::psiInv),6))
#define msquare() unaryExpr(std::bind2nd(std::pointer_to_binary_function<double,int,double>(std::pow),2))
#define msqrt() unaryExpr(std::ptr_fun<double,double>(std::sqrt))
#define mthreshold() unaryExpr(std::bind2nd(std::pointer_to_binary_function<double,double,double>(intermediate::threshold),1E-12))

// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision

using std::vector;

// SRC
namespace countMatrixFactor {

    // CONSTRUCTOR
    gamPoisFactor::gamPoisFactor(int n, int p, int K, const MatrixXi &X,
                                 const MatrixXd &phi1, const MatrixXd &phi2,
                                 const MatrixXd &theta1, const MatrixXd &theta2,
                                 const MatrixXd &alpha1, const MatrixXd &alpha2,
                                 const MatrixXd &beta1, const MatrixXd &beta2) {

        // dimensions
        m_N = n;
        m_P = p;
        m_K = K;
        m_ZI = false;

        // data
        m_X = MatrixXi(X);
        m_lambda0 = MatrixXd(X.cast<double>());
        intermediate::eraseZero(m_lambda0);

        m_lambda = MatrixXd::Zero(n, p);

        // variational parameters
        m_phi1cur = MatrixXd::Zero(n,K); //MatrixXd(phi1);
        m_phi2cur = MatrixXd::Zero(n,K); //MatrixXd(phi2);

        m_phi1old = MatrixXd::Zero(n,K); //MatrixXd(phi1);
        m_phi2old = MatrixXd::Zero(n,K); //MatrixXd(phi2);

        m_theta1cur = MatrixXd::Zero(p,K); //MatrixXd(theta1);
        m_theta2cur = MatrixXd::Zero(p,K); //MatrixXd(theta2);

        m_theta1old = MatrixXd::Zero(p,K); //MatrixXd(theta1);
        m_theta2old = MatrixXd::Zero(p,K); //MatrixXd(theta2);

        // sufficient statistics
        m_EU = MatrixXd::Zero(n,K);
        m_ElogU = MatrixXd::Zero(n,K);

        m_EV = MatrixXd::Zero(p,K);
        m_ElogV = MatrixXd::Zero(p,K);

        m_EZ_i = MatrixXd::Zero(p,K);
        m_EZ_j = MatrixXd::Zero(n,K);

        m_exp_ElogU_ElogV_k = MatrixXd::Zero(n,p);

        // prior parameter
        m_alpha1cur = MatrixXd(alpha1);
        m_alpha2cur = MatrixXd(alpha1);
        m_beta1cur = MatrixXd(beta1);
        m_beta2cur = MatrixXd(beta2);

        m_alpha1old = MatrixXd(alpha1);
        m_alpha2old = MatrixXd(alpha1);
        m_beta1old = MatrixXd(beta1);
        m_beta2old = MatrixXd(beta2);

        // order of factors
        m_orderDeviance = VectorXi::Zero(K);
        m_orderExpVar0 = VectorXi::Zero(K);
        m_orderExpVarU = VectorXi::Zero(K);
        m_orderExpVarV = VectorXi::Zero(K);

        // criterion following order of factors
        m_kDeviance = VectorXd::Zero(K);
        m_kExpVar0 = VectorXd::Zero(K);
        m_kExpVarU = VectorXd::Zero(K);
        m_kExpVarV = VectorXd::Zero(K);

    }

    // DESTRUCTOR
    gamPoisFactor::~gamPoisFactor() {}

    // GETTER
    // getter for ZI boolean
    bool gamPoisFactor::getZI() {
        return m_ZI;
    }

    // member functions: documented in src

    /*!
     * \brief Initialization of sufficient statistics
     */
    void gamPoisFactor::Init() {

        // Gamma prior parameter (to avoid scaling issue)
        m_alpha1cur = m_alpha1cur.array() / m_K;
        m_beta1cur = m_beta1cur.array() / m_K;

        // Rcpp::Rcout << "X = " << m_X << std::endl << std::endl;

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

        // sufficient statistics
        Rcpp::Rcout << "Init: sufficient statistics" << std::endl;
        Egam(m_phi1cur, m_phi2cur, m_EU);
        Elgam(m_phi1cur, m_phi2cur, m_ElogU);
        Egam(m_theta1cur, m_theta2cur, m_EV);
        Elgam(m_theta1cur, m_theta2cur, m_ElogV);

        // Multinomial Z parameters
        this->multinomParam();

        // Poisson rate
        this->poissonRate();
    }

    //-------------------//
    //      criteria     //
    //-------------------//

    /*!
     * \brief compute conditonal log-likelihood
     *
     * @return the current value (double)
     */
    double gamPoisFactor::computeCondLogLike() {
        double res = poisLogLike(m_X, m_lambda);
        return res;
    }

    /*!
    * \brief compute Prior log-likelihood
    *
    * @return the current value (double)
    */
    double gamPoisFactor::computePriorLogLike() {
        double res = gammaLogLike(m_EU, m_alpha1cur, m_alpha2cur) + gammaLogLike(m_EV, m_beta1cur, m_beta2cur);
        return res;
    }

    /*!
    * \brief compute Posterior log-likelihood
    *
    * @return the current value (double)
    */
    double gamPoisFactor::computePostLogLike() {
        double res = gammaLogLike(m_EU, m_phi1cur, m_phi2cur) + gammaLogLike(m_EV, m_theta1cur, m_theta2cur);
        return res;
    }

    /*!
    * \brief compute complete log-likelihood
    *
    * @return the current value (double)
    */
    double gamPoisFactor::computeCompLogLike() {
        double res = poisLogLike(m_X, m_lambda) + gammaLogLike(m_EU, m_phi1cur, m_phi2cur) + gammaLogLike(m_EV, m_theta1cur, m_theta2cur);
        return res;
    }

    /*!
    * \brief compute complete log-likelihood
    *
    * @return the current value (double)
    */
    double gamPoisFactor::computeMargLogLike() {
        // return m_condLogLike(iter) + m_priorLogLike(iter) - m_postLogLike(iter);
        double res = poisLogLike(m_X, m_lambda) + gammaLogLike(m_EU, m_alpha1cur, m_alpha2cur) + gammaLogLike(m_EV, m_beta1cur, m_beta2cur)
                    - gammaLogLike(m_EU, m_phi1cur, m_phi2cur) - gammaLogLike(m_EV, m_theta1cur, m_theta2cur);
        return res;
    }

    /*!
    * \brief compute evidence lower bound
    *
    * @return the current value (double)
    */
    double gamPoisFactor::computeELBO() {
        intermediate::checkExp(m_ElogU);
        intermediate::checkExp(m_ElogV);

        double res1 = (-1) * ( ( (m_X.cast<double>().array() + 1).mlgamma() ).sum() + ( m_lambda ).sum() );
        //Rcpp::Rcout << "ELBO: res1 = " << res1 << std::endl;
        double res2 = ( m_X.cast<double>().array() * (m_exp_ElogU_ElogV_k).mlog().array()).sum();
        //Rcpp::Rcout << "ELBO: res2 = " << res2 << std::endl;

        double res3 = ( (m_alpha1cur.array() - 1) * m_ElogU.array() + m_alpha1cur.array() * m_alpha2cur.mlog().array()
                            - m_alpha2cur.array() * m_EU.array() - m_alpha1cur.mlgamma().array() ).sum();
        //Rcpp::Rcout << "ELBO: res3 = " << res3 << std::endl;
        double res4 = (-1) * ( (m_phi1cur.array() - 1) * m_ElogU.array() + m_phi1cur.array() * m_phi2cur.mlog().array()
                                   - m_phi2cur.array() * m_EU.array() - m_phi1cur.mlgamma().array() ).sum();
        //Rcpp::Rcout << "ELBO: res4 = " << res4 << std::endl;

        double res5 = ( (m_beta1cur.array() - 1) * m_ElogV.array() + m_beta1cur.array() * m_beta2cur.mlog().array()
                            - m_beta2cur.array() * m_EV.array() - m_beta1cur.mlgamma().array() ).sum();
        //Rcpp::Rcout << "ELBO: res5 = " << res5 << std::endl;
        double res6 = (-1) * ( (m_theta1cur.array() - 1) * m_ElogV.array() + m_theta1cur.array() * m_theta2cur.mlog().array()
                                   - m_theta2cur.array() * m_EV.array() - m_theta1cur.mlgamma().array() ).sum();
        //Rcpp::Rcout << "ELBO: res6 = " << res6 << std::endl;

        double res = res1 + res2 + res3 + res4 + res5 + res6;

        return res;
    }

    /*!
    * \brief compute deviance between estimated and saturated model
    *
    * @return the current value (double)
    */
    double gamPoisFactor::computeDeviance() {
        double res = poisDeviance(m_X, m_lambda, m_lambda0);
        return res;
    }

    /*!
     * \brief compute explained variance regarding residuals sum of squares
     *
     * @return the current value (double)
     */
    double gamPoisFactor::computeExpVar0() {
        double res = expVar0(m_X, m_EU, m_EV);
        return res;
    }

    /*!
     * \brief compute explained variance regarding U
     *
     * @return the current value (double)
     */
    double gamPoisFactor::computeExpVarU() {
        double res = expVarU(m_X, m_EU);
        return res;
    }

    /*!
     * \brief compute explained variance regarding V
     *
     * @return the current value (double)
     */
    double gamPoisFactor::computeExpVarV() {
        double res = expVarV(m_X, m_EV);
        return res;
    }

    //--------------------------------------------//
    // parameter updates for standard variational //
    //--------------------------------------------//

    /*!
     * \brief update rule for poisson rates in variational inference
     */
    void gamPoisFactor::poissonRate() {
        m_lambda = m_EU * m_EV.transpose();
    }

    /*!
     * \brief update rule for multinomial parameters in variational inference
     *
     * m_EZ_i_{jk} = sum_i E[Z_{ijk}]
     * m_EZ_j_{ik} = sum_j E[Z_{ijk}]
     */
    void gamPoisFactor::multinomParam() {

        intermediate::checkExp(m_ElogU);
        intermediate::checkExp(m_ElogV);

        // sum_k exp(E[log(U_{ik})]) * exp(E[log(V_{jk})])
        m_exp_ElogU_ElogV_k = m_ElogU.mexp() * m_ElogV.mexp().transpose();

        // Rcpp::Rcout << "dim = " << sum_k.rows() << " x " << sum_k.cols()  << std::endl;

        m_EZ_j = m_ElogU.mexp().array() * ( (m_X.cast<double>().array() / m_exp_ElogU_ElogV_k.array() ).matrix() * m_ElogV.mexp() ).array();
        m_EZ_i = m_ElogV.mexp().array() * ( (m_X.cast<double>().array() / m_exp_ElogU_ElogV_k.array() ).matrix().transpose() * m_ElogU.mexp() ).array();

        // test
        for(int i=0; i<m_N; i++) {
            for(int k = 0; k<m_K; k++) {
                double test = 0;
                for(int j=0; j<m_P; j++) {
                    test += m_X(i,j) * std::exp(m_ElogU(i,k)) * std::exp(m_ElogV(j,k)) / m_exp_ElogU_ElogV_k (i,j);
                    // Rcpp::Rcout << "X(i,j) = " << m_X(i,j) << " / Z_{ijk} = " << m_X(i,j) * std::exp(m_ElogU(i,k)) * std::exp(m_ElogV(j,k)) / m_exp_ElogU_ElogV_k (i,j) << std::endl;
                }

                if(test != m_EZ_j(i,k)) {
                    // Rcpp::Rcout << "test = " << test << " m_EZ_j = " <<  m_EZ_j(i,k) << std::endl;
                }
            }
        }

        for(int j=0; j<m_P; j++) {
            for(int k = 0; k<m_K; k++) {
                double test = 0;
                for(int i=0; i<m_N; i++) {
                    test += m_X(i,j) * std::exp(m_ElogU(i,k)) * std::exp(m_ElogV(j,k)) / m_exp_ElogU_ElogV_k (i,j);
                }

                if(test != m_EZ_i(j,k)) {
                    // Rcpp::Rcout << "test = " << test << " m_EZ_i = " <<  m_EZ_i(j,k) << std::endl;
                }
            }
        }
    }

    /*!
     * \brief update rule for local parameters phi (factor U) in variational inference
     */
    void gamPoisFactor::localParam() {
        m_phi1cur = m_alpha1cur.array() + m_EZ_j.array();
        m_phi2cur = m_alpha2cur.rowwise() + m_EV.colwise().sum();

        // expectation and log-expectation
        Egam(m_phi1cur, m_phi2cur, m_EU);
        Elgam(m_phi1cur, m_phi2cur, m_ElogU);
    }

    /*!
     * \brief update rule for global parameters theta (factor V) in variational inference
     */
    void gamPoisFactor::globalParam() {
        m_theta1cur = m_beta1cur.array() + m_EZ_i.array();
        m_theta2cur = m_beta2cur.rowwise() + m_EU.colwise().sum();

        // expectation and log-expectation
        Egam(m_theta1cur, m_theta2cur, m_EV);
        Elgam(m_theta1cur, m_theta2cur, m_ElogV);
    }

    /*!
     * \brief parameter update in standard variational
     */
    void gamPoisFactor::updateVarational() {

        // Multinomial parameters
        // Rcpp::Rcout << "variational: Multinomial parameters" << std::endl;
        this->multinomParam();

        // Rcpp::Rcout << "sum_i E[Z_{ijk}] = " << m_EZ_i << std::endl;
        // Rcpp::Rcout << "sum_j E[Z_{ijk}] = " << m_EZ_j << std::endl << std::endl;

        // Rcpp::Rcout << "sum_k exp(E[log U_{ik}]) * exp(E[log V_{jk}]) = " << m_exp_ElogU_ElogV_k << std::endl << std::endl;

        // local parameters
        // U : param phi
        // Rcpp::Rcout << "variational: local parameters" << std::endl;
        this->localParam();

        // Rcpp::Rcout << "phi1 = " << m_phi1cur << std::endl;
        // Rcpp::Rcout << "phi2 = " << m_phi2cur << std::endl;

        // Rcpp::Rcout << "E[U_{ik}] = " << m_EU << std::endl;
        // Rcpp::Rcout << "E[log U_{ik}] = " << m_ElogU << std::endl << std::endl;

        // global parameters
        // V : param theta
        // Rcpp::Rcout << "variational: global parameters" << std::endl;
        this->globalParam();

        // Rcpp::Rcout << "theta1 = " << m_theta1cur << std::endl;
        // Rcpp::Rcout << "theta2 = " << m_theta2cur << std::endl;

        // Rcpp::Rcout << "E[V_{jk}] = " << m_EV << std::endl;
        // Rcpp::Rcout << "E[log V_{jk}] = " << m_ElogV << std::endl << std::endl;

        // Poisson rate
        // Rcpp::Rcout << "variationaml: Poisson rate" << std::endl;
        this->poissonRate();

        // Rcpp::Rcout << "sum_k E[U_{ik}] E[V_{jk}] = " << m_lambda << std::endl << std::endl;
    }

    /*!
     * \brief update parameters between iterations
     */
    void gamPoisFactor::nextIterate() {
        m_phi1old = m_phi1cur;
        m_phi2old = m_phi2cur;

        m_theta1old = m_theta1cur;
        m_theta2old = m_theta2cur;
    }

    //--------------------------------------//
    // parameter updates for variational EM //
    //--------------------------------------//

    /*!
     * \brief local parameter update: alpha (factor U)
     */
    void gamPoisFactor::localPriorParam() {
        m_alpha1cur = (m_alpha2cur.mlog().rowwise() + m_ElogU.colwise().mean()).mpsiInv();
        // m_alpha1cur = ((m_alpha2cur.mlog().rowwise() + m_ElogU.colwise().mean()).array() - 0.5).mthreshold().mpsiInv();

        m_alpha2cur = m_alpha1cur.array().rowwise() / m_EU.colwise().mean().array();
        // m_alpha2cur = (m_N * m_alpha1cur.array()).mthreshold().array().rowwise() / (m_EU.colwise().sum().array() + 10).array();

        // VectorXd ElogU(m_ElogU.colwise().mean());
        // VectorXd EU(m_EU.colwise().mean());
        // // m_beta1cur = ((m_beta2cur.mlog().array() - ((double) 0.5 / (double) m_P) * m_beta2cur.cwiseInverse().array()).rowwise() + m_ElogU.colwise().mean()).mpsiInv();
        // // m_beta2cur = (m_P * m_beta1cur.array() + ((std::pow(m_P,2) * m_beta1cur.msquare().array()) + 8 * 1 * (m_beta1cur.array().rowwise() * EU).array() ).msqrt().array()).array().rowwise() / EU.array();
        //
        // double lambda = 0;
        // int i = 0;
        // for(int k = 0; k<m_K; k++) {
        //     m_beta1cur(i,k) = intermediate::psiInv(std::log(m_beta2cur(i,k)) - lambda / ((double) m_P * m_beta2cur(i,k)) + ElogU(k), 6);
        //     m_beta2cur(i,k) = (m_beta1cur(i,k) + std::sqrt( std::pow(m_beta1cur(i,k), 2) + (double) 8.0 * lambda * EU(k) * m_beta1cur(i,k) / (double) m_P)) / (2.0 * EU(k));
        // }
        //
        // m_beta1cur = m_beta1cur.row(1).replicate(m_P,1);
        // m_beta2cur = m_beta2cur.row(1).replicate(m_P,1);
    }

    /*!
     * \brief global parameter update: beta (factor V)
     */
    void gamPoisFactor::globalPriorParam() {
        m_beta1cur = (m_beta2cur.mlog().rowwise() + m_ElogV.colwise().mean()).mpsiInv();
        // m_beta1cur = ((m_beta2cur.mlog().rowwise() + m_ElogV.colwise().mean()).array() - 0.5).mthreshold().mpsiInv();

        m_beta2cur = m_beta1cur.array().rowwise() / m_EV.colwise().mean().array();
        // m_beta2cur = (m_P * m_beta1cur.array()).mthreshold().array().rowwise() / (m_EV.colwise().sum().array() + 10).array();

        // VectorXd ElogV(m_ElogV.colwise().mean());
        // VectorXd EV(m_EV.colwise().mean());
        // m_beta1cur = ((m_beta2cur.mlog().array() - ((double) 0.5 / (double) m_P) * m_beta2cur.cwiseInverse().array()).rowwise() + m_ElogV.colwise().mean()).mpsiInv();
        // m_beta2cur = (m_P * m_beta1cur.array() + ((std::pow(m_P,2) * m_beta1cur.msquare().array()) + 8 * 0.5 * (m_beta1cur.array().rowwise() * EV).array() ).msqrt().array()).array().rowwise() / EV.array();

        // double lambda = 0;
        // int i = 0;
        // for(int k = 0; k<m_K; k++) {
        //     m_alpha1cur(i,k) = intermediate::psiInv(std::log(m_alpha2cur(i,k)) - lambda / ((double) m_N * m_alpha2cur(i,k)) + ElogV(k), 6);
        //     m_alpha2cur(i,k) = (m_alpha1cur(i,k) + std::sqrt( std::pow(m_alpha1cur(i,k), 2) + (double) 8.0 * lambda * EV(k) * m_alpha1cur(i,k) / (double) m_N)) / (2.0 * EV(k));
        // }
        //
        // m_alpha1cur = m_alpha1cur.row(1).replicate(m_N,1);
        // m_alpha2cur = m_alpha2cur.row(1).replicate(m_N,1);
    }

    /*!
     * \brief parameter update in variational EM (E-step)
     */
    void gamPoisFactor::updateEstep() {
        this->updateVarational();
    }

    /*!
     * \brief parameter update in variational EM (M-step)
     */
    void gamPoisFactor::updateMstep() {
        // local parameters
        // U : param phi
        // Rcpp::Rcout << "algorithm: local parameters" << std::endl;
        this->localPriorParam();

        // Rcpp::Rcout << "alpha1 = " << m_alpha1cur << std::endl;
        // Rcpp::Rcout << "alpha2 = " << m_alpha2cur << std::endl << std::endl;

        // global parameters
        // V : param theta
        // Rcpp::Rcout << "algorithm: global parameters" << std::endl;
        this->globalPriorParam();

        // Rcpp::Rcout << "beta1 = " << m_beta1cur << std::endl;
        // Rcpp::Rcout << "beta2 = " << m_beta2cur << std::endl << std::endl;


    }

    // update parameters between iterations in Estep
    void gamPoisFactor::nextIterateEstep() {
        m_phi1old = m_phi1cur;
        m_phi2old = m_phi2cur;

        m_theta1old = m_theta1cur;
        m_theta2old = m_theta2cur;
    }

    // update parameters between iterations in Mstep
    void gamPoisFactor::nextIterateMstep() {
        m_alpha1old = m_alpha1cur;
        m_alpha2old = m_alpha2cur;

        m_beta1old = m_beta1cur;
        m_beta2old = m_beta2cur;
    }


    //-------------------//
    //     algorithm     //
    //-------------------//

    /*!
     * \brief compute normalized gap between two iterates
     */
    double gamPoisFactor::normGap() {
        double paramNorm = sqrt(intermediate::parameterNorm2(m_phi1old, m_phi2old)
                                    + intermediate::parameterNorm2(m_theta1old, m_theta2old));
        double diffNorm = sqrt(intermediate::differenceNorm2(m_phi1old, m_phi2old, m_phi1cur, m_phi2cur)
                                   + intermediate::differenceNorm2(m_theta1old, m_theta2old, m_theta1cur, m_theta2cur));

        double res = diffNorm / paramNorm;
        return res;
    }

    /*!
     * \brief compute normalized gap between two iterates in E-step of EM algo
     */
    double gamPoisFactor::normGapEstep() {
        double paramNorm = sqrt(intermediate::parameterNorm2(m_phi1old, m_phi2old)
                                    + intermediate::parameterNorm2(m_theta1old, m_theta2old));
        double diffNorm = sqrt(intermediate::differenceNorm2(m_phi1old, m_phi2old, m_phi1cur, m_phi2cur)
                                   + intermediate::differenceNorm2(m_theta1old, m_theta2old, m_theta1cur, m_theta2cur));

        double res = diffNorm / paramNorm;
        return res;
    }

    /*!
     * \brief compute normalized gap between two iterates in M-step of EM algo
     */
    double gamPoisFactor::normGapMstep() {
        double paramNorm = sqrt(intermediate::parameterNorm2(m_alpha1old, m_alpha2old)
                                   + intermediate::parameterNorm2(m_beta1old, m_beta2old));
        double diffNorm = sqrt(intermediate::differenceNorm2(m_alpha1old, m_alpha2old, m_alpha1cur, m_alpha2cur)
                                   + intermediate::differenceNorm2(m_beta1old, m_beta2old, m_beta1cur, m_beta2cur));

        double res = diffNorm / paramNorm;
        return res;
    }

    /*!
     * \brief compute normalized gap between two iterates in EM algo
     */
    double gamPoisFactor::normGapEM() {
        double paramNorm = sqrt(intermediate::parameterNorm2(m_phi1old, m_phi2old)
                                    + intermediate::parameterNorm2(m_theta1old, m_theta2old)
                                    + intermediate::parameterNorm2(m_alpha1old, m_alpha2old)
                                    + intermediate::parameterNorm2(m_beta1old, m_beta2old));
        double diffNorm = sqrt(intermediate::differenceNorm2(m_phi1old, m_phi2old, m_phi1cur, m_phi2cur)
                                   + intermediate::differenceNorm2(m_theta1old, m_theta2old, m_theta1cur, m_theta2cur)
                                   + intermediate::differenceNorm2(m_alpha1old, m_alpha2old, m_alpha1cur, m_alpha2cur)
                                   + intermediate::differenceNorm2(m_beta1old, m_beta2old, m_beta1cur, m_beta2cur));

        double res = diffNorm / paramNorm;
        return res;
    }

    //-------------------//
    //   order factors   //
    //-------------------//

    /*!
    * \brief order factors according to expVar0
    *
    * @param[out] vector of factor order
    */
    void gamPoisFactor::orderExpVar0(VectorXi &order) {
        vector<int> leftIndex;
        vector<int> chosenIndex;

        MatrixXd Utmp = MatrixXd::Zero(m_EU.rows(), m_EU.cols());
        MatrixXd Vtmp = MatrixXd::Zero(m_EV.rows(), m_EV.cols());

        // init
        for(int k=0; k<m_K; k++) {
            leftIndex.push_back(k);
        }

        // ordering
        for(int k=0; k<m_K; k++) {
            double val_max = -1E6;
            int ind_max = 0;
            int ind_left = 0;
            for(int ind=0; ind<leftIndex.size(); ind++) {
                MatrixXd tmpU(m_EU.rows(), k+1);
                MatrixXd tmpV(m_EV.rows(), k+1);

                tmpU.col(k) = m_EU.col(leftIndex[ind]);
                tmpV.col(k) = m_EV.col(leftIndex[ind]);

                if(k!=0) {
                    tmpU.leftCols(k) = Utmp.leftCols(k);
                    tmpV.leftCols(k) = Vtmp.leftCols(k);
                }

                double res = expVar0(m_X, tmpU, tmpV);
                if(res > val_max) {
                    val_max = res;
                    ind_max = leftIndex[ind];
                    ind_left = ind;
                }
            }
            // max found
            chosenIndex.push_back(ind_max);
            leftIndex.erase(leftIndex.begin() + ind_left);

            Utmp.col(k) = m_EU.col(ind_max);
            Vtmp.col(k) = m_EV.col(ind_max);

            // store it
            m_kExpVar0(k) = val_max;
        }

        // return
        for(int k=0; k<m_K; k++) {
            order(k) = chosenIndex[k] + 1;
        }
    }

    /*!
     * \brief order factors according to expVarU
     *
     * @param[out] vector of factor order
     */
    void gamPoisFactor::orderExpVarU(VectorXi &order) {
        vector<int> leftIndex;
        vector<int> chosenIndex;

        MatrixXd Utmp = MatrixXd::Zero(m_EU.rows(), m_EU.cols());

        // init
        for(int k=0; k<m_K; k++) {
            leftIndex.push_back(k);
        }

        // ordering
        for(int k=0; k<m_K; k++) {
            double val_max = -1E6;
            int ind_max = 0;
            int ind_left = 0;
            for(int ind=0; ind<leftIndex.size(); ind++) {
                MatrixXd tmpU(m_EU.rows(), k+1);
                tmpU.col(k) = m_EU.col(leftIndex[ind]);

                if(k!=0) {
                    tmpU.leftCols(k) = Utmp.leftCols(k);
                }

                double res = expVarU(m_X, tmpU);
                if(res > val_max) {
                    val_max = res;
                    ind_max = leftIndex[ind];
                    ind_left = ind;
                }
            }
            // max found
            chosenIndex.push_back(ind_max);
            leftIndex.erase(leftIndex.begin() + ind_left);

            Utmp.col(k) = m_EU.col(ind_max);

            // store it
            m_kExpVarU(k) = val_max;
        }

        // return
        for(int k=0; k<m_K; k++) {
            order(k) = chosenIndex[k] + 1;
        }
    }

    /*!
     * \brief order factors according to expVarV
     *
     * @param[out] vector of factor order
     */
    void gamPoisFactor::orderExpVarV(VectorXi &order) {
        vector<int> leftIndex;
        vector<int> chosenIndex;

        MatrixXd Vtmp = MatrixXd::Zero(m_EV.rows(), m_EV.cols());

        // init
        for(int k=0; k<m_K; k++) {
            leftIndex.push_back(k);
        }

        // ordering
        for(int k=0; k<m_K; k++) {
            double val_max = -1E6;
            int ind_max = 0;
            int ind_left = 0;
            for(int ind=0; ind<leftIndex.size(); ind++) {
                MatrixXd tmpV(m_EV.rows(), k+1);
                tmpV.col(k) = m_EV.col(leftIndex[ind]);

                if(k!=0) {
                    tmpV.leftCols(k) = Vtmp.leftCols(k);
                }

                double res = expVarV(m_X, tmpV);
                if(res > val_max) {
                    val_max = res;
                    ind_max = leftIndex[ind];
                    ind_left = ind;
                }
            }
            // max found
            chosenIndex.push_back(ind_max);
            leftIndex.erase(leftIndex.begin() + ind_left);

            Vtmp.col(k) = m_EV.col(ind_max);

            // store it
            m_kExpVarV(k) = val_max;
        }

        // return
        for(int k=0; k<m_K; k++) {
            order(k) = chosenIndex[k] + 1;
        }
    }

    /*!
     * \brief order factors according to deviance
     *
     * @param[out] vector of factor order
     */
    void gamPoisFactor::orderDeviance(VectorXi &order) {
        vector<int> leftIndex;
        vector<int> chosenIndex;

        MatrixXd Utmp = MatrixXd::Zero(m_EU.rows(), m_EU.cols());
        MatrixXd Vtmp = MatrixXd::Zero(m_EV.rows(), m_EV.cols());

        // init
        for(int k=0; k<m_K; k++) {
            leftIndex.push_back(k);
        }

        // ordering
        for(int k=0; k<m_K; k++) {
            double val_min = 0;
            int ind_min = 0;
            int ind_left = 0;
            for(int ind=0; ind<leftIndex.size(); ind++) {

                MatrixXd tmpU(m_EU.rows(), k+1);
                MatrixXd tmpV(m_EV.rows(), k+1);

                tmpU.col(k) = m_EU.col(leftIndex[ind]);
                tmpV.col(k) = m_EV.col(leftIndex[ind]);

                if(k!=0) {
                    tmpU.leftCols(k) = Utmp.leftCols(k);
                    tmpV.leftCols(k) = Vtmp.leftCols(k);
                }

                MatrixXd lambda = tmpU * tmpV.transpose();

                //Rcpp::Rcout << "min lambda = " << lambda.minCoeff() << std::endl;
                //Rcpp::Rcout << "max lambda = " << lambda.maxCoeff() << std::endl;

                //Rcpp::Rcout << "min lambda0 = " << m_lambda0.minCoeff() << std::endl;
                //Rcpp::Rcout << "max lambda0 = " << m_lambda0.maxCoeff() << std::endl;

                double res = poisDeviance(m_X, lambda, m_lambda0);
                if(ind==0) {
                    val_min = res;
                    ind_min = leftIndex[ind];
                    ind_left = ind;
                }
                if(res < val_min) {
                    val_min = res;
                    ind_min = leftIndex[ind];
                    ind_left = ind;
                }
            }
            // max found
            chosenIndex.push_back(ind_min);
            leftIndex.erase(leftIndex.begin() + ind_left);

            Utmp.col(k) = m_EU.col(ind_min);
            Vtmp.col(k) = m_EV.col(ind_min);

            // store it
            m_kDeviance(k) = val_min;

        }

        // return
        for(int k=0; k<m_K; k++) {
            order(k) = chosenIndex[k] + 1;
        }
    }

    /*!
     * \brief compute factor order
     */
    void gamPoisFactor::computeOrder() {
        this->orderExpVar0(m_orderExpVar0);
        this->orderExpVarU(m_orderExpVarU);
        this->orderExpVarV(m_orderExpVarV);
        this->orderDeviance(m_orderDeviance);
    }


    //-------------------//
    //       return      //
    //-------------------//

    /*!
     * \brief create list with results to be return
     *
     * @param[out] list containing output
     */
    void gamPoisFactor::returnObject(Rcpp::List &results) {

        Rcpp::List params = Rcpp::List::create(Rcpp::Named("phi1") = m_phi1cur,
                                               Rcpp::Named("phi2") = m_phi2cur,
                                               Rcpp::Named("theta1") = m_theta1cur,
                                               Rcpp::Named("theta2") = m_theta2cur,
                                               Rcpp::Named("alpha1") = m_alpha1cur,
                                               Rcpp::Named("alpha2") = m_alpha2cur,
                                               Rcpp::Named("beta1") = m_beta1cur,
                                               Rcpp::Named("beta2") = m_beta2cur);

        Rcpp::List stats = Rcpp::List::create(Rcpp::Named("EU") = m_EU,
                                              Rcpp::Named("EV") = m_EV,
                                              Rcpp::Named("ElogU") = m_ElogU,
                                              Rcpp::Named("ElogV") = m_ElogV);

        Rcpp::List order = Rcpp::List::create(Rcpp::Named("orderDeviance") = m_orderDeviance,
                                              Rcpp::Named("orderExpVar0") = m_orderExpVar0,
                                              Rcpp::Named("orderExpVarU") = m_orderExpVarU,
                                              Rcpp::Named("orderExpVarV") = m_orderExpVarV);

        Rcpp::List criteria_k = Rcpp::List::create(Rcpp::Named("kDeviance") = m_kDeviance,
                                                  Rcpp::Named("kExpVar0") = m_kExpVar0,
                                                  Rcpp::Named("kExpVarU") = m_kExpVarU,
                                                  Rcpp::Named("kExpVarV") = m_kExpVarV);

        Rcpp::List returnObj = Rcpp::List::create(Rcpp::Named("U") = m_EU,
                                                  Rcpp::Named("V") = m_EV,
                                                  Rcpp::Named("lambda") = m_lambda,
                                                  Rcpp::Named("params") = params,
                                                  Rcpp::Named("stats") = stats,
                                                  Rcpp::Named("order") = order,
                                                  Rcpp::Named("criteria_k") = criteria_k,
                                                  Rcpp::Named("ZI") = m_ZI);

        SEXP tmp = Rcpp::Language("c", results, returnObj).eval();

        results = tmp;
    }
}