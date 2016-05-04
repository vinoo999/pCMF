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
     */
    void gamPoisFactorStandard::ELBO(int iter) {}

    /*!
     * \brief compute deviance between estimated and saturated model
     */
    void gamPoisFactorStandard::deviance(int iter) {}

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

            // deviance
            this->deviance(iter);


            //// explained variance
            explainedVar0(iter) = expVar0(X, EU, EV);
            explainedVar1(iter) = expVar1(X, EU);
            explainedVar2(iter) = expVar2(X, EV);


            //// breaking condition: convergence or not
            // Rcout << "convergence" << std::endl;
            double paramNorm = parameterNorm(phi1old, phi2old, theta1old, theta2old);
            double diffNorm = differenceNorm(phi1old, phi2old, theta1old, theta2old, phi1cur, phi2cur, theta1cur, theta2cur);

            normGap(iter) = diffNorm / paramNorm;

            // derivative order to consider
            double condition = convCondition(order, normGap, iter, 0);

            if(std::abs(condition) < epsilon) {
                nstab++;
            } else {
                nstab=0;
            }

            if(nstab > rstab) {
                converged=true;
                nbIter=iter;
            }


            //// increment values of parameters
            phi1old = phi1cur;
            phi2old = phi2cur;

            theta1old = theta1cur;
            theta2old = theta2cur;

            //// increment iteration
            iter++;

        }

    }

}