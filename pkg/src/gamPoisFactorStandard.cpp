// Copyright 2016-04 Ghislain Durif
//
// This file is part of the `hClustering' library for R and related languages.
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
 * \file gamPoisFactorStandard.cpp
 * \brief class definition for standard Gamma Poisson Factor Model
 * \author Ghislain Durif
 * \version 0.1
 * \date 04/05/2016
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include <stdio>
#include <boost/math/special_functions/digamma.hpp>
#include "gamPoisFactorStandard.h"

#define exp() unaryExpr(std::ptr_fun<double,double>(exp))
#define digamma() unaryExpr(std::ptr_fun<double,double>(digamma))
#define lgamma() unaryExpr(std::ptr_fun<double,double>(lgamma))
#define log() unaryExpr(std::ptr_fun<double,double>(log))
#define square() unaryExpr(std::bind2nd(std::pointer_to_binary_function<double,double,double>(pow),2))

// [[Rcpp::depends(BH)]]
using boost::math::digamma;

// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision

// SRC
namespace countMatrixFactor {

    // CONSTRUCTOR
    gamPoisFactorStandard::gamPoisFactorStandard(int n, int p, int K, int iterMax, int order,
                                                 int stabRange, double epsilon, bool verbose,
                                                 const MatrixXi &X,
                                                 const MatrixXd &phi1, const MatrixXd &phi2,
                                                 const MatrixXd &theta1, const MatrixXd &theta2,
                                                 const MatrixXd &alpha1, const MatrixXd &alpha2,
                                                 const MatrixXd &beta1, const MatrixXd &beta2) :
                gamPoisFactor::gamPoisFactor(n, p, K, iterMax, order,
                                              stabRange, epsilon, verbose,
                                              X, phi1, phi2, theta1, theta2,
                                              alpha1, alpha2, beta1, beta2) {}

    // DESTRUCTOR
    gamPoisFactorStandard::~gamPoisFactorStandard() {}

    // member functions: documented in src

    /*!
     * \brief Initialization of sufficient statistics
     */
    void gamPoisFactorStandard::Init() {

        // Gamma variational parameter
        Egam(m_phi1cur, m_phi2cur, m_EU);
        Elgam(m_phi1cur, m_phi2cur, m_ElogU);
        Egam(m_theta1cur, m_theta2cur, m_EV);
        Elgam(m_theta1cur, m_theta2cur, m_ElogV);

        // Multinomial Z parameters
        this->multinomParam();

    }

    //-------------------//
    //      criteria     //
    //-------------------//

    /*!
     * \brief compute log-likelihood
     *
     * @param[in] iter current iteration
     */
    void gamPoisFactorStandard::computeLogLike(int iter) {
        m_condLogLike(iter) = poisLogLike(m_X, m_lambda);
        m_priorLogLike(iter) = gammaLogLike(m_EU, m_alpha1, m_alpha2) + gammaLogLike(m_EV, m_beta1, m_beta2);
        m_postLogLike(iter) = gammaLogLike(m_EU, m_phi1cur, m_phi2cur) + gammaLogLike(m_EV, m_theta1cur, m_theta2cur);
        m_compLogLike(iter) = m_condLogLike(iter) + m_postLogLike(iter);
        m_margLogLike(iter) = m_condLogLike(iter) + m_priorLogLike(iter) - m_postLogLike(iter);
    }

    /*!
     * \brief compute evidence lower bound
     *
     * @param[in] iter current iteration
     */
    void gamPoisFactorStandard::computeELBO(int iter) {
        m_elbo(iter) = this->ELBO();
    }

    /*!
     * \brief evidence lower bound for the specific gamma Poisson factor model
     *
     * @return value of the ELBO for the current values of estimates
     */
    double gamPoisFactorStandard::ELBO() {

        intermediate::checkExp(m_ElogU);
        intermediate::checkExp(m_ElogV);

        double res1 = (-1) * ( ( (m_X.array() + 1).lgamma() ).sum() + ( m_EU * m_EV.transpose() ).sum() );
        double res2 = ( m_X.array() * (m_ElogU.exp() * m_ElogV.exp().transpose()).log().array()).sum();

        double res3 = ( (m_alpha1.array() - 1) * m_ElogU.array() + m_alpha1.array() * m_alpha2.log().array()
                        - m_alpha2.array() * m_EU.array() - m.alpha1.lgamma() ).sum();
        double res4 = (-1) * ( (m_phi1.array() - 1) * m_ElogU.array() + m_phi1.array() * m_phi2.log().array()
                                   - m_phi2.array() * m_EU.array() - m.phi1.lgamma().array() ).sum();

        double res5 = ( (m_beta1.array() - 1) * m_ElogV.array() + m_beta1.array() * m_beta2.log().array()
                            - m_beta2.array() * m_EV.array() - m.beta1.lgamma() ).sum();
        double res6 = (-1) * ( (m_theta1.array() - 1) * m_ElogV.array() + m_theta1.array() * m_theta2.log().array()
                                   - m_theta2.array() * m_EV.array() - m.theta1.lgamma().array() ).sum();

        double res = res1 + res2 + res3 + res4 + res5 + res6;
    }

    /*!
     * \brief compute deviance between estimated and saturated model
     *
     * @param[in] iter current iteration
     */
    void gamPoisFactorStandard::computeDeviance(int iter) {
        m_deviance(iter) = poisDeviance(m_X, m_lambda, m_lambda0);
    }

    //-------------------//
    // parameter updates //
    //-------------------//

    // poisson rate
    void gamPoisFactorStandard::poissonRate() {
        m_lambda = m_EU * m_EV.transpose();
    }

    // Poisson intensity
    void gamPoisFactorStandard::multinomParam() {

        intermediate::checkExp(m_ElogU);
        intermediate::checkExp(m_ElogV);

        m_EZ_j = m_ElogU.exp().array() * ( (m_X.cast<double>().array() / (m_ElogU.exp() * m_ElogV.exp().transpose()).array() ).matrix() * m_ElogV.exp() ).array();
        m_EZ_j = m_ElogU.exp().array() * ( (m_X.cast<double>().array() / (m_ElogU.exp() * m_ElogV.exp().transpose()).array() ).matrix() * m_ElogV.exp() ).array();
    }

    // local parameters: phi (factor U)
    void gamPoisFactorStandard::localParam() {
        m_phi1cur = m_alpha1.array() + m_EZ_j.array();
        m_phi2cur = m_alpha2.rowwise() + m_EV.colwise().sum();
    }

    // global parameters: theta (factor V)
    void gamPoisFactorStandard::globalParam() {
        m_theta1cur = m_beta1.array() + m_EZ_i.array();
        m_theta2cur = m_beta2.rowwise() + m_EU.colwise().sum();
    }


    //-------------------//
    //     algorithm     //
    //-------------------//
    void gamPoisFactorStandard::algorithm() {
        // Iteration
        int nstab = 0; // number of successive iteration where the normalized gap betwwen two iteration is close to zero (convergence when nstab > rstab)
        int iter = 0;

        while( (iter < m_iterMax) && (m_converged==false)) {

            if(m_verbose==true) {
                Rcpp::Rcout << "iter " << iter << std::endl;
            }

            // Multinomial parameters
            this->multinomParam();

            // local parameters
            // U : param phi
            this->localParam();

            // expectation and log-expectation
            Egam(m_phi1cur, m_phi2cur, m_EU);
            Elgam(m_phi1cur, m_phi2cur, m_ElogU);

            // global parameters
            // V : param theta
            this->globalParam();

            // expectation and log-expectation
            Egam(m_theta1cur, m_theta2cur, m_EV);
            Elgam(m_theta1cur, m_theta2cur, m_ElogV);

            // Poisson rate
            this->poissonRate();

            // log-likelihood
            this->computeLogLike(iter);
            // ELBO
            this->computeELBO(iter);
            // deviance
            this->computeDeviance(iter);
            // explained variance
            this->computeExpVar(iter);
            // convergence
            this-> assessConvergence(iter, nstab);
            // increment values of parameters
            this->updates();
            // increment iteration
            iter++;
        }
    }

    //-------------------//
    //       return      //
    //-------------------//

    /*!
     * \brief create list with results to be return
     *
     * @param[out] list containing output
     */
    Rcpp::List gamPoisFactorStandard::returnObject(Rcpp::List &results) {
        Rcpp::List logLikelihood = Rcpp:List::create(Rcpp::Named("margLogLike") = m_margLogLike.head(m_nbIter),
                                                     Rcpp::Named("condLogLike") = m_condLogLike.head(m_nbIter),
                                                     Rcpp::Named("priorLogLike") = m_priorLogLike.head(m_nbIter),
                                                     Rcpp::Named("postLogLike") = m_postLogLike.head(m_nbIter),
                                                     Rcpp::Named("compLogLike") = m_compLogLike.head(m_nbIter),
                                                     Rcpp::Named("elbo") = m_elbo.head(nbIter));

        Rcpp::List expVariance = Rcpp:List::create(Rcpp::Named("expVar0") = m_expVar0.head(m_nbIter),
                                                   Rcpp::Named("expVarU") = m_expVarU.head(m_nbIter),
                                                   Rcpp::Named("expVarV") = m_expVarV.head(m_nbIter));

        Rcpp::List params = Rcpp:List::create(Rcpp::Named("phi1") = m_phi1cur,
                                              Rcpp::Named("phi2") = m_phi2cur,
                                              Rcpp::Named("theta1") = m_theta1cur,
                                              Rcpp::Named("theta2") = m_theta2cur);

        Rcpp::List stats = Rcpp:List::create(Rcpp::Named("EU") = m_EU,
                                             Rcpp::Named("EV") = m_EV,
                                             Rcpp::Named("ElogU") = m_ElogU,
                                             Rcpp::Named("ElogV") = m_ElogV);

        Rcpp::List returnObj = Rcpp::List::create(Rcpp::Named("U") = m_EU,
                                                  Rcpp::Named("V") = m_EV,
                                                  Rcpp::Named("logLikelihood") = logLikelihood,
                                                  Rcpp::Named("expVariance") = expVariance,
                                                  Rcpp::Named("params") = params,
                                                  Rcpp::Named("stats") = stats,
                                                  Rcpp::Named("normGap") = m_normGap.head(m_nbIter),
                                                  Rcpp::Named("deviance") = m_deviance.head(m_nbIter),
                                                  Rcpp::Named("converged") = m_converged,
                                                  Rcpp::Named("nbIter") = m_nbIter);

        results = returnObj;
    }
}