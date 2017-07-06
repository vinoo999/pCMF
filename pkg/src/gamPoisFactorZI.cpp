// Copyright 2016-05 Ghislain Durif
//
// This file is part of the `pCMF' library for R and related languages.
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
* \brief class implementation for zero-inflated Gamma Poisson Factor Model
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
#include "random.h"

#include <omp.h>
// [[Rcpp::plugins(openmp)]]

#define mdirac() unaryExpr(std::ptr_fun<double,double>(intermediate::dirac))
#define mexp() unaryExpr(std::ptr_fun<double,double>(std::exp))
#define mlgamma() unaryExpr(std::ptr_fun<double,double>(lgamma))
#define mlog() unaryExpr(std::ptr_fun<double,double>(std::log))
#define msquare() unaryExpr(std::bind2nd(std::pointer_to_binary_function<double,double,double>(std::pow),2))

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
        m_probZIprior = VectorXd::Zero(p);
        m_probZI = MatrixXd::Zero(n,p);
        m_freqZI = VectorXd::Zero(p);

    }

    // DESTRUCTOR
    gamPoisFactorZI::~gamPoisFactorZI() {}

    // member functions: documented in src

    /*!
    * \brief Initialization of sufficient statistics
    */
    void gamPoisFactorZI::Init(myRandom::RNGType rng) {

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
                estimParam(1000, m_alpha1cur(i,k), m_alpha2cur(i,k),
                           param1, param2, rng);
                m_phi1cur(i,k) = param1;
                m_phi1old(i,k) = param1;
                m_phi2cur(i,k) = param2;
                m_phi2old(i,k) = param2;
            }
            // global parameters
            for(int j=0; j<m_P; j++) {
                double param1 = 0;
                double param2 = 0;
                estimParam(1000, m_beta1cur(j,k), m_beta2cur(j,k),
                           param1, param2, rng);
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
                    m_freqZI(j) = m_freqZI(j)+1;
                    m_probZI(i,j) = 1;
                } else {
                    m_probZI(i,j) = 0;
                }
            }
        }

        m_freqZI = m_freqZI.array() / m_N;
        m_probZIprior = m_freqZI;

        // Rcpp::Rcout << "m_freqZI = " << m_freqZI << std::endl;
        // Rcpp::Rcout << "m_freqZI.transpose() = " << m_freqZI.head(10).transpose().replicate(m_N,1) << std::endl;

        // m_probZI = m_freqZI.transpose().replicate(m_N,1);

        // Rcpp::Rcout << "m_probZI = " << m_probZI.leftCols(15) << std::endl;

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
     * \brief compute evidence lower bound
     *
     * @return the current value (double)
     */
    double gamPoisFactorZI::computeELBO() {
        // intermediate::checkExp(m_ElogU);
        // intermediate::checkExp(m_ElogV);

        // sum_k exp(E[log(U_{ik})]) * exp(E[log(V_{jk})])
        // m_exp_ElogU_ElogV_k = m_ElogU.mexp() * m_ElogV.mexp().transpose();
        // for(int i=0; i<m_N; i++) {
        //     for(int j=0; j<m_P; j++) {
        //         double res = 0;
        //         for(int k=0; k<m_K; k++) {
        //             res += (m_ElogU(i,k) + m_ElogV(j,k) >= -300 ? std::exp(m_ElogU(i,k) + m_ElogV(j,k)) : 0);
        //         }
        //         m_exp_ElogU_ElogV_k(i,j) = res;
        //     }
        // }
        int i, j, k;
        MatrixXd *ElogU = &m_ElogU;
        MatrixXd *ElogV = &m_ElogV;
        MatrixXd *exp_ElogU_ElogV_k = &m_exp_ElogU_ElogV_k;
        #if defined(_OPENMP)
        #pragma omp parallel for private(j,k)
        #endif
        for(i=0; i<m_N; i++) {
            for(j=0; j<m_P; j++) {
                double res = 0;
                for(k=0; k<m_K; k++) {
                    res += (ElogU->coeffRef(i,k) + ElogV->coeffRef(j,k) >= -300 ? std::exp(ElogU->coeffRef(i,k) + ElogV->coeffRef(j,k)) : 0);
                }
                exp_ElogU_ElogV_k->coeffRef(i,j) = res;
            }
        }

        double resFinal = 0;

        VectorXd omega = VectorXd::Zero(m_K);

        // concerning Z
        double res1 = 0;
        double res2 = 0;

        for(i=0; i<m_N; i++) {
            for(j=0; j<m_P; j++) {
                for(k=0; k<m_K; k++) {
                    if(m_exp_ElogU_ElogV_k(i,j) > 0) {
                        omega(k) = (m_ElogU(i,k) + m_ElogV(j,k) >= -300 ? std::exp(m_ElogU(i,k) + m_ElogV(j,k)) : 0) / m_exp_ElogU_ElogV_k(i,j);
                    } else {
                        omega(k) = 0;
                    }
                }
                for(k=0; k<m_K; k++) {
                    res1 += m_probZI(i,j) * (m_X(i,j) * omega(k) * (m_ElogU(i,k) + m_ElogV(j,k))
                                                     - m_EU(i,k) * m_EV(j,k));
                    res2 += m_X(i,j) * omega(k) * (m_exp_ElogU_ElogV_k(i,j) > 0 ? m_ElogU(i,k) + m_ElogV(j,k) - std::log(m_exp_ElogU_ElogV_k(i,j)) : 0);
                }
            }
        }
        resFinal += res1 - res2;

        // regarding U
        double res3 = ( (m_alpha1cur.array() - 1) * m_ElogU.array() + m_alpha1cur.array() * m_alpha2cur.mlog().array()
                            - m_alpha2cur.array() * m_EU.array() - m_alpha1cur.mlgamma().array() ).sum();
        // Rcpp::Rcout << "ELBO: res3 = " << res3 << std::endl;
        double res4 = ( (m_phi1cur.array() - 1) * m_ElogU.array() + m_phi1cur.array() * m_phi2cur.mlog().array()
                                   - m_phi2cur.array() * m_EU.array() - m_phi1cur.mlgamma().array() ).sum();
        // Rcpp::Rcout << "ELBO: res4 = " << res4 << std::endl;
        resFinal += res3 - res4;

        // regarding V
        double res5 = ( (m_beta1cur.array() - 1) * m_ElogV.array() + m_beta1cur.array() * m_beta2cur.mlog().array()
                            - m_beta2cur.array() * m_EV.array() - m_beta1cur.mlgamma().array() ).sum();
        // Rcpp::Rcout << "ELBO: res5 = " << res5 << std::endl;
        double res6 = (-1) * ( (m_theta1cur.array() - 1) * m_ElogV.array() + m_theta1cur.array() * m_theta2cur.mlog().array()
                                   - m_theta2cur.array() * m_EV.array() - m_theta1cur.mlgamma().array() ).sum();
        // Rcpp::Rcout << "ELBO: res6 = " << res6 << std::endl;
        resFinal += res5 - res6;

        // regarding D
        double res7 = 0;
        double res8 = 0;
        for(j=0; j<m_P; j++) {
            for(k=0; k<m_K; k++) {
                if((m_probZIprior(j)>0) && (m_probZIprior(j)<1)) {
                    res7 += m_probZI(j,k) * std::log(m_probZIprior(j)) + (1-m_probZI(j,k)) * std::log(1-m_probZIprior(j));
                }
                if((m_probZI(j,k)>0) && (m_probZI(j,k)<1)) {
                    res8 += m_probZI(j,k) * std::log(m_probZI(j,k)) + (1-m_probZI(j,k)) * std::log(1-m_probZI(j,k));
                }
            }
        }
        resFinal += res7 - res8;

        return(resFinal);


        // // double res1 = (-1) * ( ( (m_X.cast<double>().array() + 1).mlgamma() ).sum() + ( m_lambda ).sum() );
        // double res1 = (-1) * ( (m_probZI.array() * m_lambda.array()).sum() );
        // // Rcpp::Rcout << "ELBO: res1 = " << res1 << std::endl;
        // double res2 = ( (m_probZI.array() * m_X.cast<double>().array()) * (m_exp_ElogU_ElogV_k).mlog().array()).sum();
        // // Rcpp::Rcout << "ELBO: res2 = " << res2 << std::endl;
        //
        // double res2a =  (-1) * (m_probZI.array() * (m_X.cast<double>().array() + 1).mlgamma().array()).sum();
        //
        // double res3 = ( (m_alpha1cur.array() - 1) * m_ElogU.array() + m_alpha1cur.array() * m_alpha2cur.mlog().array()
        //                     - m_alpha2cur.array() * m_EU.array() - m_alpha1cur.mlgamma().array() ).sum();
        // // Rcpp::Rcout << "ELBO: res3 = " << res3 << std::endl;
        // double res4 = (-1) * ( (m_phi1cur.array() - 1) * m_ElogU.array() + m_phi1cur.array() * m_phi2cur.mlog().array()
        //                            - m_phi2cur.array() * m_EU.array() - m_phi1cur.mlgamma().array() ).sum();
        // // Rcpp::Rcout << "ELBO: res4 = " << res4 << std::endl;
        //
        // double res5 = ( (m_beta1cur.array() - 1) * m_ElogV.array() + m_beta1cur.array() * m_beta2cur.mlog().array()
        //                     - m_beta2cur.array() * m_EV.array() - m_beta1cur.mlgamma().array() ).sum();
        // // Rcpp::Rcout << "ELBO: res5 = " << res5 << std::endl;
        // double res6 = (-1) * ( (m_theta1cur.array() - 1) * m_ElogV.array() + m_theta1cur.array() * m_theta2cur.mlog().array()
        //                            - m_theta2cur.array() * m_EV.array() - m_theta1cur.mlgamma().array() ).sum();
        // // Rcpp::Rcout << "ELBO: res6 = " << res6 << std::endl;
        //
        // // regarding D
        // double res7 = 0;
        // double res8 = 0;
        // for(int j=0; j<m_P; j++) {
        //     for(int k=0; k<m_K; k++) {
        //         if((m_probZIprior(j)>0) && (m_probZIprior(j)<1)) {
        //             res7 += m_probZI(j,k) * std::log(m_probZIprior(j)) + (1-m_probZI(j,k)) * std::log(1-m_probZIprior(j));
        //         }
        //         if((m_probZI(j,k)>0) && (m_probZI(j,k)<1)) {
        //             res8 += m_probZI(j,k) * std::log(m_probZI(j,k)) + (1-m_probZI(j,k)) * std::log(1-m_probZI(j,k));
        //         }
        //     }
        // }
        //
        // double res = res1 + res2 + res2a + res3 + res4 + res5 + res6 + res7 - res8;
        //
        // return res;
    }

    //--------------------------------------------//
    // parameter updates for standard variational //
    //--------------------------------------------//

    /*!
     * \brief update rule for poisson rates in variational inference
     */
    void gamPoisFactorZI::poissonRate() {
        m_lambda = m_probZI.array() * (m_EU * m_EV.transpose()).array();
        // m_lambda = m_EU * m_EV.transpose();
    }


    /*!
     * \brief update rule for multinomial parameters in variational inference
     *
     * m_EZ_i_{jk} = sum_i pi_j * E[Z_{ijk}]
     * m_EZ_j_{ik} = sum_j pi_j * E[Z_{ijk}]
     */
    void gamPoisFactorZI::multinomParam() {

        // intermediate::checkExp(m_ElogU);
        // intermediate::checkExp(m_ElogV);

        // sum_k exp(E[log(U_{ik})]) * exp(E[log(V_{jk})])
        // m_exp_ElogU_ElogV_k = m_ElogU.mexp() * m_ElogV.mexp().transpose();
        // for(int i=0; i<m_N; i++) {
        //     for(int j=0; j<m_P; j++) {
        //         double res = 0;
        //         for(int k=0; k<m_K; k++) {
        //             res += (m_ElogU(i,k) + m_ElogV(j,k) >= -300 ? std::exp(m_ElogU(i,k) + m_ElogV(j,k)) : 0);
        //         }
        //         m_exp_ElogU_ElogV_k(i,j) = res;
        //     }
        // }
        int i, j, k;
        MatrixXi *X = &m_X;
        MatrixXd *ElogU = &m_ElogU;
        MatrixXd *ElogV = &m_ElogV;
        MatrixXd *exp_ElogU_ElogV_k = &m_exp_ElogU_ElogV_k;
        MatrixXd *EZ_i = &m_EZ_i;
        MatrixXd *EZ_j = &m_EZ_j;
        MatrixXd *probZI = &m_probZI;
        #if defined(_OPENMP)
        #pragma omp parallel for private(j,k)
        #endif
        for(i=0; i<m_N; i++) {
            for(j=0; j<m_P; j++) {
                double res = 0;
                for(k=0; k<m_K; k++) {
                    res += (ElogU->coeffRef(i,k) + ElogV->coeffRef(j,k) >= -300 ? std::exp(ElogU->coeffRef(i,k) + ElogV->coeffRef(j,k)) : 0);
                }
                exp_ElogU_ElogV_k->coeffRef(i,j) = res;
            }
        }

        // sum_j E[Z_{ijk}] * p_{i,j}
        // sum_i E[Z_{ijk}] * p_{i,j}
        // m_EZ_j = m_ElogU.mexp().array() * ( ((m_X.cast<double>().array() * m_probZI.msquare().array()) / m_exp_ElogU_ElogV_k.array() ).matrix() * m_ElogV.mexp() ).array();
        // m_EZ_i = m_ElogV.mexp().array() * ( ((m_X.cast<double>().array() * m_probZI.msquare().array()) / m_exp_ElogU_ElogV_k.array() ).matrix().transpose() * m_ElogU.mexp() ).array();
        #if defined(_OPENMP)
        #pragma omp parallel for private(j,k)
        #endif
        for(i=0; i<m_N; i++) {
            for(k = 0; k<m_K; k++) {
                double res = 0;
                for(j=0; j<m_P; j++) {
                    if(exp_ElogU_ElogV_k->coeffRef(i,j)>0) {
                        res += probZI->coeffRef(i,j) * X->coeffRef(i,j) * (ElogU->coeffRef(i,k) + ElogV->coeffRef(j,k) >= -300 ? std::exp(ElogU->coeffRef(i,k) + ElogV->coeffRef(j,k)) : 0) / exp_ElogU_ElogV_k->coeffRef(i,j);
                    }
                }
                EZ_j->coeffRef(i,k) = res;
            }
        }
        #if defined(_OPENMP)
        #pragma omp parallel for private(i,k)
        #endif
        for(j=0; j<m_P; j++) {
            for(k = 0; k<m_K; k++) {
                double res = 0;
                for(i=0; i<m_N; i++) {
                    if(exp_ElogU_ElogV_k->coeffRef(i,j)>0) {
                        res += probZI->coeffRef(i,j) * X->coeffRef(i,j) * (ElogU->coeffRef(i,k) + ElogV->coeffRef(j,k) >= -300 ? std::exp(ElogU->coeffRef(i,k) + ElogV->coeffRef(j,k)) : 0) / exp_ElogU_ElogV_k->coeffRef(i,j);
                    }
                }
                EZ_i->coeffRef(j,k) = res;
            }
        }
    }

    /*!
     * \brief update rule for local parameters phi (factor U) in variational inference
     */
    void gamPoisFactorZI::localParam() {
        m_phi1cur = m_alpha1cur.array() + m_EZ_j.array();
        m_phi2cur = m_alpha2cur.array() + (m_probZI * m_EV).array();

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
        m_theta2cur = m_beta2cur.array() + (m_probZI.transpose() * m_EU).array();

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

        MatrixXd lambda = m_EU * m_EV.transpose();

        // Rcpp::Rcout << "m_probZI" << std::endl;
        int i, j;
        MatrixXi *X = &m_X;
        MatrixXd *ElogU = &m_ElogU;
        MatrixXd *ElogV = &m_ElogV;
        MatrixXd *exp_ElogU_ElogV_k = &m_exp_ElogU_ElogV_k;
        MatrixXd *EZ_i = &m_EZ_i;
        MatrixXd *EZ_j = &m_EZ_j;
        MatrixXd *probZI = &m_probZI;
        VectorXd *probZIprior = &m_probZIprior;
        #if defined(_OPENMP)
        #pragma omp parallel for private(j)
        #endif
        for(i= 0; i<m_N; i++) {
            for(j=0; j<m_P; j++) {
                if(X->coeffRef(i,j) != 0) {
                    probZI->coeffRef(i,j) = 1;
                } else {
                    if(probZIprior->coeffRef(j) == 1) {
                        probZI->coeffRef(i,j) = 1;
                    } else if(probZIprior->coeffRef(j) == 0) {
                        probZI->coeffRef(i,j) = 0;
                    } else {
                        // m_probZI(i,j) = intermediate::threshold(intermediate::expit( intermediate::logit(m_probZIprior(j)) - m_lambda(i,j)),1E-12);
                        probZI->coeffRef(i,j) = intermediate::expit( intermediate::logit(probZIprior->coeffRef(j)) - lambda(i,j));
                    }
                }
            }
        }

        // Rcpp::Rcout << "m_probZI = " << m_probZI.leftCols(15) << std::endl;
    }

    /*!
     * \brief parameter update in standard variational
     */
    void gamPoisFactorZI::updateVarational() {

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

        // Poisson rate
        Rcpp::Rcout << "algorithm: Poisson rate" << std::endl;
        this->poissonRate();

        // ZI proba
        Rcpp::Rcout << "algorithm: ZI proba" << std::endl;
        this->ZIproba();
    }

    //--------------------------------------//
    // parameter updates for variational EM //
    //--------------------------------------//

    /*!
     * \brief update rule for Bernoulli parameter (of ZI indicator) in prior
     */
    void gamPoisFactorZI::priorZIproba() {
        m_probZIprior = m_probZI.colwise().mean();
    }

    /*!
     * \brief parameter update in variational EM (E-step)
     */
    void gamPoisFactorZI::updateEstep() {

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

        // Poisson rate
        // Rcpp::Rcout << "algorithm: Poisson rate" << std::endl;
        this->poissonRate();

        // ZI proba
        //Rcpp::Rcout << "algorithm: ZI proba" << std::endl;
        this->ZIproba();
    }

    /*!
     * \brief parameter update in variational EM (M-step)
     */
    void gamPoisFactorZI::updateMstep(int iter) {
        m_curIter = iter;

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

        Rcpp::List ZI = Rcpp::List::create(Rcpp::Named("prob") = m_probZI,
                                           Rcpp::Named("freq") = m_freqZI,
                                           Rcpp::Named("probPrior") = m_probZIprior);

        Rcpp::List returnObj2 = Rcpp::List::create(Rcpp::Named("ZIparams") = ZI);

        SEXP tmp1 = Rcpp::Language("c", returnObj1, returnObj2).eval();

        SEXP tmp2 = Rcpp::Language("c", results, tmp1).eval();

        results = tmp2;

    }
}
