// Copyright 2016-08 Ghislain Durif
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
* \file gamPoisFactorSparse.cpp
* \brief class implementation for sparse Gamma Poisson Factor Model
* \author Ghislain Durif
* \version 0.1
* \date 02/08/2016
*/

#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <cstdio>
#include <vector>
#include "gamPoisFactorSparse.h"
#include "gamDistrib.h"
#include "intermediate.h"

#define mdirac() unaryExpr(std::ptr_fun<double,double>(intermediate::dirac))
#define mexp() unaryExpr(std::ptr_fun<double,double>(intermediate::safeExp))
#define mlgamma() unaryExpr(std::ptr_fun<double,double>(lgamma))
#define mlog() unaryExpr(std::ptr_fun<double,double>(std::log))
#define mlogit() unaryExpr(std::ptr_fun<double,double>(intermediate::logit))
#define mlogitinv() unaryExpr(std::ptr_fun<double,double>(intermediate::logitinv))
#define mpsiInv() unaryExpr(std::bind2nd(std::pointer_to_binary_function<double,int,double>(intermediate::psiInv),6))

// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision

using std::vector;

// SRC
namespace countMatrixFactor {

    // CONSTRUCTOR
    gamPoisFactorSparse::gamPoisFactorSparse(int n, int p, int K, const MatrixXi &X,
                                             const MatrixXd &phi1, const MatrixXd &phi2,
                                             const MatrixXd &theta1, const MatrixXd &theta2,
                                             const MatrixXd &alpha1, const MatrixXd &alpha2,
                                             const MatrixXd &beta1, const MatrixXd &beta2)
        : gamPoisFactor::gamPoisFactor(n, p, K, X, phi1, phi2, theta1, theta2,
          alpha1, alpha2, beta1, beta2)
        {

            // sparse probabilities and frequencies
            m_probSparsePrior = VectorXd::Zero(p);
            m_probSparse = MatrixXd::Zero(p,K);
            m_S = MatrixXd::Zero(p,K);

            // sufficient stats (sum on i)
            m_EZ_logU_i = MatrixXd::Zero(p,K);
            m_EZ_logV_i = MatrixXd::Zero(p,K);
            m_EU_EV_i = MatrixXd::Zero(p,K);
            m_ElgamZ_i = MatrixXd::Zero(p,K);
        }

    // DESTRUCTOR
    gamPoisFactorSparse::~gamPoisFactorSparse() {}

    // member functions: documented in src

    /*!
    * \brief Initialization of sufficient statistics
    */
    void gamPoisFactorSparse::Init() {

        // // Gamma prior parameter (to avoid scaling issue)
        // m_alpha2cur = m_alpha2cur.array() * std::sqrt(m_K);
        // m_beta2cur = m_beta2cur.array() * std::sqrt(m_K);

        // Gamma prior parameter (to avoid scaling issue)
        m_alpha1cur = m_alpha1cur.array() / std::sqrt(m_K);
        m_beta1cur = m_beta1cur.array() / std::sqrt(m_K);

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
            // Rcpp::NumericVector sample = Rcpp::runif(m_P);
            // VectorXd sample2 = Rcpp::as<Map<VectorXd> >(sample);
            // m_probSparse.col(k) = sample2;
            for(int k=0; k<m_K; k++) {
                m_probSparse(j,k) = 0.9;
                m_S(j,k) = 1;
            }
            m_probSparsePrior(j) = 0.9;
        }

        // Rcpp::Rcout << "phi1 = " <<  m_phi1cur << std::endl << std::endl;
        // Rcpp::Rcout << "phi2 = " <<  m_phi2cur << std::endl << std::endl;

        // Rcpp::Rcout << "theta1 = " <<  m_theta1cur << std::endl << std::endl;
        // Rcpp::Rcout << "theta2 = " <<  m_theta2cur << std::endl << std::endl;

        // sufficient statistics
        Rcpp::Rcout << "Init: sufficient statistics" << std::endl;
        Egam(m_phi1cur, m_phi2cur, m_EU);
        Elgam(m_phi1cur, m_phi2cur, m_ElogU);
        Egam(m_theta1cur, m_theta2cur, m_EV);
        Elgam(m_theta1cur, m_theta2cur, m_ElogV);

        // Rcpp::Rcout << "EU = " <<  m_EU << std::endl << std::endl;
        // Rcpp::Rcout << "ElogU = " <<  m_ElogU << std::endl << std::endl;
        // Rcpp::Rcout << "EV = " <<  m_EV << std::endl << std::endl;
        // Rcpp::Rcout << "ElogV = " <<  m_ElogV << std::endl << std::endl;

        // Multinomial Z parameters
        this->multinomParam();

    }

    //-------------------//
    //      criteria     //
    //-------------------//

    /*!
     * \brief compute evidence lower bound
     *
     * @return the current value (double)
     */
    double gamPoisFactorSparse::computeELBO() {
        // intermediate::checkExp(m_ElogU);
        // intermediate::checkExp(m_ElogV);

        VectorXd omega = VectorXd::Zero(m_K);

        double resFinal = 0;

        // sum_{ijk} E[log p(Z_{ijk} | U_{ik}, V_{jk}, S_{jk})] - E[q(Z_{ijk})]
        double res1 = 0;
        double res2 = 0;
        for(int i=0; i<m_N; i++) {
            for(int j=0; j<m_P; j++) {
                // for(int k=0; k<m_K; k++) {
                //     if(m_exp_ElogU_ElogV_k(i,j) > 0) {
                //         omega(k) = m_S(j,k) * std::exp(m_ElogU(i,k) + m_ElogV(j,k)) / m_exp_ElogU_ElogV_k(i,j);
                //     } else {
                //         omega(k) = 0;
                //     }
                // }
                // for(int k=0; k<m_K; k++) {
                //     res1 += m_probSparse(j,k) * (m_X(i,j) * omega(k) * (m_ElogU(i,k) + m_ElogV(j,k))
                //                                      - m_EU(i,k) * m_EV(j,k));
                //     res2 += m_X(i,j) * omega(k) * std::log(omega(k)>0 ? omega(k) : 1);
                // }
                for(int k=0; k<m_K; k++) {
                    res1 += m_probSparse(j,k) * ( m_X(i,j) * std::log(m_exp_ElogU_ElogV_k(i,j)>0 ? m_exp_ElogU_ElogV_k(i,j) : 1) - m_EU(i,k) * m_EV(j,k));
                }
            }
        }
        resFinal += res1 - res2;

        // regarding U
        double res3 = ( (m_alpha1cur.array() - 1) * m_ElogU.array() + m_alpha1cur.array() * m_alpha2cur.mlog().array()
                            - m_alpha2cur.array() * m_EU.array() - m_alpha1cur.mlgamma().array() ).sum();
        double res4 = ( (m_phi1cur.array() - 1) * m_ElogU.array() + m_phi1cur.array() * m_phi2cur.mlog().array()
                                   - m_phi2cur.array() * m_EU.array() - m_phi1cur.mlgamma().array() ).sum();
        resFinal += res3 - res4;

        // regarding V
        double res5 = 0;
        for(int j=0; j<m_P; j++) {
            for(int k=0; k<m_K; k++) {
                // res5 += m_probSparse(j,k) * ( (m_beta1cur(j,k) - 1) * m_ElogV(j,k) + m_beta1cur(j,k) * std::log(m_beta2cur(j,k))
                //                                   - m_beta2cur(j,k) * m_EV(j,k) - lgamma(m_beta1cur(j,k)));
                res5 += ( (m_beta1cur(j,k) - 1) * m_ElogV(j,k) + m_beta1cur(j,k) * std::log(m_beta2cur(j,k))
                                                  - m_beta2cur(j,k) * m_EV(j,k) - lgamma(m_beta1cur(j,k)));
            }
        }
        double res6 = ( (m_theta1cur.array() - 1) * m_ElogV.array() + m_theta1cur.array() * m_theta2cur.mlog().array()
                                   - m_theta2cur.array() * m_EV.array() - m_theta1cur.mlgamma().array() ).sum();
        resFinal += res5 - res6;

        // regarding S
        double res7 = 0;
        double res8 = 0;
        for(int j=0; j<m_P; j++) {
            for(int k=0; k<m_K; k++) {
                if((m_probSparsePrior(j)>0) && (m_probSparsePrior(j)<1)) {
                    res7 += m_probSparse(j,k) * std::log(m_probSparsePrior(j)) + (1-m_probSparse(j,k)) * std::log(1-m_probSparsePrior(j));
                }
                if((m_probSparse(j,k)>0) && (m_probSparse(j,k)<1)) {
                    res8 += m_probSparse(j,k) * std::log(m_probSparse(j,k)) + (1-m_probSparse(j,k)) * std::log(1-m_probSparse(j,k));
                }
            }
        }
        resFinal += res7 - res8;

        return resFinal;
    }

    //-------------------//
    // parameter updates //
    //-------------------//

    /*!
     * \brief update rule for poisson rates in variational inference
     */
    void gamPoisFactorSparse::poissonRate() {
        m_lambda = m_EU * (m_S.array() * m_EV.array()).matrix().transpose();
    }

    /*!
    * \brief update rule for multinomial parameters in variational inference
    *
    * m_EZ_i_{jk} = sum_i pi_j * E[Z_{ijk}]
    * m_EZ_j_{ik} = sum_j pi_j * E[Z_{ijk}]
    */
    void gamPoisFactorSparse::multinomParam() {

        // intermediate::checkExp(m_ElogU);
        // intermediate::checkExp(m_ElogV);

        // MatrixXd test = MatrixXd::Zero(m_N,m_P);

        // sum_k exp(E[log(U_{ik})]) * exp(E[log(V_{jk})])
        for(int i=0; i<m_N; i++) {
            for(int j=0; j<m_P; j++) {
                // double rhoMin = ( m_ElogU.row(i) + m_ElogV.row(j) ).minCoeff();
                // double rhoMax = ( m_ElogU.row(i) + m_ElogV.row(j) ).maxCoeff();
                double res = 0;
                for(int k=0; k<m_K; k++) {
                    if(m_S(j,k) > 0) {
                        res += std::exp(m_ElogU(i,k) + m_ElogV(j,k));
                    }
                    // res += m_probSparse(j,k) * std::exp(m_ElogU(i,k) + m_ElogV(j,k));
                    // test(i,j) += std::exp(m_ElogU(i,k) + m_ElogV(j,k));
                }
                m_exp_ElogU_ElogV_k(i,j) = res;
            }
        }


        // Rcpp::Rcout << "sum exp_ElogU_ElogV_k = " <<  m_exp_ElogU_ElogV_k << std::endl << std::endl;

        // sum_j E[Z_{ijk}]
        // sum_i E[Z_{ijk}]
        // m_EZ_j = m_ElogU.mexp().array() * ( (m_X.cast<double>().array() / m_exp_ElogU_ElogV_k.array() ).matrix() * m_ElogV.mexp() ).array();
        // m_EZ_i = m_ElogV.mexp().array() * ( (m_X.cast<double>().array() / m_exp_ElogU_ElogV_k.array() ).matrix().transpose() * m_ElogU.mexp() ).array();

        // sum_j E[Z_{ijk}]
        for(int i=0; i<m_N; i++) {
            for(int k = 0; k<m_K; k++) {
                double res = 0;
                double test0 = 0;
                for(int j=0; j<m_P; j++) {
                    // if(m_S(j,k) > 0) {
                    //     if(m_exp_ElogU_ElogV_k(i,j) > 0) {
                    //         res += m_X(i,j) * std::exp(m_ElogU(i,k) + m_ElogV(j,k)) / m_exp_ElogU_ElogV_k(i,j);
                    //     }
                    // }
                    // test0 += m_X(i,j) * std::exp(m_ElogU(i,k) + m_ElogV(j,k)) / test(i,j);
                    if(m_exp_ElogU_ElogV_k(i,j) > 0) {
                        // res += m_probSparse(j,k) * m_X(i,j) * std::exp(m_ElogU(i,k) + m_ElogV(j,k)) / m_exp_ElogU_ElogV_k(i,j);
                        res += m_probSparse(j,k) * m_X(i,j) * m_S(j,k) * std::exp(m_ElogU(i,k) + m_ElogV(j,k)) / m_exp_ElogU_ElogV_k(i,j);
                    }
                }
                // Rcpp::Rcout << "value computed = " << res << std::endl;
                // Rcpp::Rcout << "test = " << test0 << std::endl << std::endl;
                m_EZ_j(i,k) = res;
            }
        }

        // sum_i E[Z_{ijk}]
        for(int j=0; j<m_P; j++) {
            for(int k = 0; k<m_K; k++) {
                double res = 0;
                for(int i=0; i<m_N; i++) {
                    // if(m_S(j,k) > 0) {
                    //     if(m_exp_ElogU_ElogV_k(i,j) > 0) {
                    //         res += m_X(i,j) * std::exp(m_ElogU(i,k) + m_ElogV(j,k)) / m_exp_ElogU_ElogV_k(i,j);
                    //     }
                    // }
                    if(m_exp_ElogU_ElogV_k(i,j) > 0) {
                        res += m_probSparse(j,k) * m_X(i,j) * m_S(j,k) * std::exp(m_ElogU(i,k) + m_ElogV(j,k)) / m_exp_ElogU_ElogV_k(i,j);
                    }
                }
                m_EZ_i(j,k) = res;
            }
        }

        // Rcpp::Rcout << "sum EZ_i = " <<  m_EZ_i << std::endl << std::endl;
        // Rcpp::Rcout << "sum EZ_j = " <<  m_EZ_j << std::endl << std::endl;
    }

    /*!
    * \brief update rule for local parameters phi (factor U) in variational inference
    */
    void gamPoisFactorSparse::localParam() {
        m_phi1cur = m_alpha1cur.array() + m_EZ_j.array();
        // m_phi2cur = m_alpha2cur.array().rowwise() + (m_EV.array() * m_probSparse.array()).colwise().sum();

        // test
        for(int i=0; i<m_N; i++) {
            for(int k = 0; k<m_K; k++) {
                m_phi2cur(i,k) = m_alpha2cur(i,k) + m_probSparse.col(k).dot(m_EV.col(k));
            }
        }

        // Rcpp::Rcout << "phi1 = " <<  m_phi1cur << std::endl << std::endl;
        // Rcpp::Rcout << "phi2 = " <<  m_phi2cur << std::endl << std::endl;

        // expectation and log-expectation
        Egam(m_phi1cur, m_phi2cur, m_EU);
        Elgam(m_phi1cur, m_phi2cur, m_ElogU);

        // Rcpp::Rcout << "EU = " <<  m_EU << std::endl << std::endl;
        // Rcpp::Rcout << "ElogU = " <<  m_ElogU << std::endl << std::endl;
    }

    /*!
    * \brief update rule for global parameters theta (factor V) in variational inference
    */
    void gamPoisFactorSparse::globalParam() {
        m_theta1cur = m_beta1cur.array() + m_EZ_i.array();
        // m_theta2cur = m_beta2cur.array() + (m_probSparse.array().rowwise() * m_EU.colwise().sum().array()).array();

        // test
        for(int j=0; j<m_P; j++) {
            for(int k = 0; k<m_K; k++) {
                m_theta2cur(j,k) = m_beta2cur(j,k) + m_probSparse(j,k) * m_EU.col(k).sum();
            }
        }

        // Rcpp::Rcout << "theta1 = " <<  m_theta1cur << std::endl << std::endl;
        // Rcpp::Rcout << "theta2 = " <<  m_theta2cur << std::endl << std::endl;


        // expectation and log-expectation
        Egam(m_theta1cur, m_theta2cur, m_EV);
        Elgam(m_theta1cur, m_theta2cur, m_ElogV);

        // Rcpp::Rcout << "EV = " <<  m_EV << std::endl << std::endl;
        // Rcpp::Rcout << "ElogV = " <<  m_ElogV << std::endl << std::endl;
    }

    /*!
    * \brief update rule for Bernoulli parameter (of ZI indicator) in variational inference
    */
    void gamPoisFactorSparse::Sproba() {

        // m_EZ_logU_i = MatrixXd::Zero(m_P,m_K);
        // m_EZ_logV_i = MatrixXd::Zero(m_P,m_K);
        // m_EU_EV_i = MatrixXd::Zero(m_P,m_K);
        // m_ElgamZ_i = MatrixXd::Zero(m_P,m_K);

        // multinomial probabilities
        VectorXd omega = VectorXd::Zero(m_N);

        // Rcpp::Rcout << "m_probSparsePrior = " <<  m_probSparsePrior << std::endl;

        for(int j=0; j<m_P; j++) {
            for(int k=0; k<m_K; k++) {

                // m_EZ_logU_i;       /*!< p x K, \sum_i E[Z_{ijk}] * E[log U_{ik}] */
                // m_EZ_logV_i;       /*!< p x K, \sum_i E[Z_{ijk}] * E[log V_{jk}] */
                // m_EU_EV_i;         /*!< p x K, \sum_i E[U_{ik}] * E[V_{jk}] */
                // m_ElgamZ_i;        /*!< p x K, \sum_i E[log(Z_{ijk}!)] */
                m_EZ_logU_i(j,k) = 0;
                m_EZ_logV_i(j,k) = 0;
                m_EU_EV_i(j,k) = 0;
                // m_ElgamZ_i(j,k) = 0;

                for(int i=0; i<m_N; i++) {
                    if(m_exp_ElogU_ElogV_k(i,j)>0) {
                        omega(i) = m_S(j,k) * std::exp(m_ElogU(i,k) + m_ElogV(j,k)) / m_exp_ElogU_ElogV_k(i,j);
                    } else {
                        omega(i) = 0;
                    }
                    m_EZ_logU_i(j,k) += m_X(i,j) * omega(i) * m_ElogU(i,k);
                    m_EZ_logV_i(j,k) += m_X(i,j) * omega(i) * m_ElogV(j,k);
                    m_EU_EV_i(j,k) += m_EU(i,k) * m_EV(j,k);
                    // m_ElgamZ_i(j,k) += intermediate::lgamBinom(m_X(i,j), omega(k));
                    // Rcpp::Rcout << "omega(i) = " << omega(i) << std::endl;
                }
                //
                // Rcpp::Rcout << "m_EU_EV_i(j,k) = " << m_EU_EV_i(j,k) << std::endl;
                // Rcpp::Rcout << "m_EZ_logU_i(j,k) = " << m_EZ_logU_i(j,k) << std::endl;
                // Rcpp::Rcout << "m_EZ_logV_i(j,k) = " << m_EZ_logV_i(j,k) << std::endl;
                // Rcpp::Rcout << "m_ElgamZ_i(j,k) = " << m_ElgamZ_i(j,k) << std::endl  << std::endl;

                // Rcpp::Rcout << "####### (j,k) = " << j << " , " << k << std::endl << std::endl;

                if(m_probSparsePrior(j) == 1) {
                    m_probSparse(j,k) = 1;
                } else if(m_probSparsePrior(j) == 0) {
                    m_probSparse(j,k) = 0;
                } else {
                    double res1 = - m_EU_EV_i(j,k) + m_EZ_logU_i(j,k) + m_EZ_logV_i(j,k); // - m_ElgamZ_i(j,k);
                    // Rcpp::Rcout << "(m_beta1cur(j,k) - 1) * m_ElogV(j,k) = " << (m_beta1cur(j,k) - 1) * m_ElogV(j,k) << std::endl;
                    // Rcpp::Rcout << "m_beta1cur(j,k) * std::log(m_beta2cur(j,k)) = " << m_beta1cur(j,k) * std::log(m_beta2cur(j,k)) << std::endl;
                    // Rcpp::Rcout << "- m_beta2cur(j,k) * m_EV(j,k) = " << - m_beta2cur(j,k) * m_EV(j,k) << std::endl;
                    // Rcpp::Rcout << "- lgamma(m_beta1cur(j,k)) = " << - lgamma(m_beta1cur(j,k)) << std::endl  << std::endl;
                    double res2 = (m_beta1cur(j,k) - 1) * m_ElogV(j,k)
                                    + m_beta1cur(j,k) * std::log(m_beta2cur(j,k))
                                    - m_beta2cur(j,k) * m_EV(j,k)
                                    - lgamma(m_beta1cur(j,k));
                    double res = res1; // + res2;
                    // Rcpp::Rcout << "from the Poisson = " <<  res1 << std::endl;
                    // Rcpp::Rcout << "from the Gamma = " <<  res2 << std::endl;
                    // Rcpp::Rcout << "term to correct the expit = " <<  res << std::endl;
                    double tmp = intermediate::expit( intermediate::logit(m_probSparsePrior(j)) + res);
                    // if(tmp==1) {
                    //     m_probSparse(j,k) = 1 - 1E-12;
                    // } else {
                    //     m_probSparse(j,k) = tmp;
                    // }
                    m_probSparse(j,k) = tmp;
                    // Rcpp::Rcout << "proba? = " <<  intermediate::expit( intermediate::logit(m_probSparsePrior(j)) + res) << std::endl;
                    // Rcpp::Rcout << "computed value = " <<  m_probSparse(j,k) << std::endl << std::endl;
                    m_S(j,k) = m_probSparse(j,k) > 0.5 ? 1 : 0;
                }
            }
        }
        // for(int j=0; j<m_P; j++) {
        //     for(int k=0; k<m_K; k++) {
        //         if(m_probSparsePrior(j) == 1) {
        //             m_probSparse(j,k) = 1;
        //         } else if(m_probSparsePrior(j) == 0) {
        //             m_probSparse(j,k) = 0;
        //         } else {
        //             double res = (m_beta1cur(j,k) - 1) * m_ElogV(j,k)
        //                             + m_beta1cur(j,k) * std::log(m_beta2cur(j,k))
        //                             - m_beta2cur(j,k) * m_EV(j,k)
        //                             - lgamma(m_beta1cur(j,k));
        //             // Rcpp::Rcout << "term to correct the expit = " <<  res << std::endl;
        //             // Rcpp::Rcout << "logit prior = " <<  intermediate::logit(m_probSparsePrior(k)) << std::endl;
        //             // Rcpp::Rcout << "proba = " <<  intermediate::expit( intermediate::logit(m_probSparsePrior(j)) - res) << std::endl;
        //             m_probSparse(j,k) = intermediate::threshold(intermediate::expit( intermediate::logit(m_probSparsePrior(j)) - res),1E-12);
        //             m_S(j,k) = m_probSparse(j,k) > 0.5 ? 1 : 0;
        //             // Rcpp::Rcout << "S(j,k) = " <<  m_S(j,k) << std::endl;
        //         }
        //     }
        // }
        // Rcpp::Rcout << "m_probSparse = " <<  m_probSparse << std::endl << std::endl;
    }

    /*!
    * \brief parameter update in standard variational
    */
    void gamPoisFactorSparse::updateVarational() {

        // Multinomial parameters
        Rcpp::Rcout << "algorithm: Multinomial parameters" << std::endl;
        this->multinomParam();

        // local parameters
        // U : param phi
        Rcpp::Rcout << "algorithm: local parameters" << std::endl;
        this->localParam();

        // global parameters
        // V : param theta
        Rcpp::Rcout << "algorithm: global parameters" << std::endl;
        this->globalParam();

        // sparse proba
        Rcpp::Rcout << "algorithm: sparse proba" << std::endl;
        this->Sproba();

        // Poisson rate
        Rcpp::Rcout << "algorithm: Poisson rate" << std::endl;
        this->poissonRate();
    }

    //--------------------------------------//
    // parameter updates for variational EM //
    //--------------------------------------//

    /*!
     * \brief update rule for Bernoulli parameter (of ZI indicator) in prior
     */
    void gamPoisFactorSparse::priorSproba() {
        m_probSparsePrior = m_probSparse.rowwise().mean();
        // Rcpp::Rcout << "m_probSparsePrior = " <<  m_probSparsePrior << std::endl << std::endl;
    }


    /*!
     * \brief local parameter update: alpha (factor U)
     */
    void gamPoisFactorSparse::localPriorParam() {
        m_alpha1cur = (m_alpha2cur.mlog().rowwise() + m_ElogU.colwise().mean()).mpsiInv();
        m_alpha2cur = m_alpha1cur.array().rowwise() / m_EU.colwise().mean().array();
    }

    /*!
     * \brief global parameter update: beta (factor V)
     */
    void gamPoisFactorSparse::globalPriorParam() {
        m_beta1cur = (m_beta2cur.mlog().rowwise() + m_ElogV.colwise().mean()).mpsiInv();
        m_beta2cur = m_beta1cur.array().rowwise() / m_EV.colwise().mean().array();

        // VectorXd probSum(m_probSparse.colwise().sum());
        //
        // int j0=0;
        // for(int k=0; k<m_K; k++) {
        //     double res1 = 0;
        //     double res2 = 0;
        //     for(int j=0; j<m_P; j++) {
        //         res1 += m_probSparse(j,k) * m_ElogV(j,k);
        //         res2 += m_probSparse(j,k) * m_EV(j,k);
        //     }
        //     m_beta1cur(j0,k) = intermediate::psiInv(std::log(m_beta2cur(j0,k)) + res1/probSum(k), 6);
        //     m_beta2cur(j0,k) = m_beta1cur(j0,k) * probSum(k) / res2;
        // }
        // m_beta1cur = m_beta1cur.row(j0).replicate(m_P,1);
        // m_beta2cur = m_beta2cur.row(j0).replicate(m_P,1);
    }

    /*!
    * \brief parameter update in variational EM (E-step)
    */
    void gamPoisFactorSparse::updateEstep() {

        // Multinomial parameters
        // Rcpp::Rcout << "algorithm: Multinomial parameters" << std::endl;
        this->multinomParam();

        // local parameters
        // U : param phi
        // Rcpp::Rcout << "algorithm: local parameters" << std::endl;
        this->localParam();

        // global parameters
        // V : param theta
        // Rcpp::Rcout << "algorithm: global parameters" << std::endl;
        this->globalParam();

        // sparse proba
        //Rcpp::Rcout << "algorithm: sparse proba" << std::endl;
        this->Sproba();

        // Poisson rate
        // Rcpp::Rcout << "algorithm: Poisson rate" << std::endl;
        this->poissonRate();
    }

    /*!
    * \brief parameter update in variational EM (M-step)
    */
    void gamPoisFactorSparse::updateMstep(int iter) {
        m_curIter = iter;

        // sparse proba
        //Rcpp::Rcout << "algorithm: sparse proba prior" << std::endl;
        this->priorSproba();

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
    void gamPoisFactorSparse::returnObject(Rcpp::List &results) {

        Rcpp::List returnObj1;
        gamPoisFactor::returnObject(returnObj1);

        Rcpp::List sparse = Rcpp::List::create(Rcpp::Named("probSparse") = m_probSparse,
                                               Rcpp::Named("probSparsePrior") = m_probSparsePrior,
                                               Rcpp::Named("sparseIndic") = m_S);

        Rcpp::List returnObj2 = Rcpp::List::create(Rcpp::Named("sparseParams") = sparse);

        SEXP tmp1 = Rcpp::Language("c", returnObj1, returnObj2).eval();

        SEXP tmp2 = Rcpp::Language("c", results, tmp1).eval();

        results = tmp2;

    }
}