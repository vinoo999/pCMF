// Copyright 2016-06 Ghislain Durif
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
* \file gamPoisFactorEM.cpp
* \brief class definition for EM algorithm in Gamma Poisson Factor Model
* \author Ghislain Durif
* \version 0.1
* \date 02/06/2016
*/

#include <Rcpp.h>
#include <RcppEigen.h>
#include "gamPoisFactorEM.h"
#include "intermediate.h"

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision

namespace countMatrixFactor {

    // CONSTRUCTOR
    gamPoisFactorEM::gamPoisFactorEM(int n, int p, int K, int iterMax,
                    int iterMax_Estep, int iterMax_Mstep, int order,
                    int stabRange, double epsilon, bool verbose,
                    const MatrixXi &X,
                    const MatrixXd &phi1, const MatrixXd &phi2,
                    const MatrixXd &theta1, const MatrixXd &theta2,
                    const MatrixXd &alpha1, const MatrixXd &alpha2,
                    const MatrixXd &beta1, const MatrixXd &beta2)
        : gamPoisFactorStandard::gamPoisFactorStandard(n, p, K, iterMax, order,
                                                       stabRange, epsilon, verbose,
                                                       X, phi1, phi2, theta1, theta2,
                                                       alpha1, alpha2, beta1, beta2)
    {
        m_iterMax_Estep = iterMax_Estep;
        m_iterMax_Vstep = iterMax_Mstep;

        m_converged_Estep = false;
        m_converged_Mstep = false;

        m_normGap_Estep = VectorXd::Zero(iterMax * iterMax_Estep);
        m_normGap_Mstep = VectorXd::Zero(iterMax * iterMax_Mstep);

        // prior parameters
        m_alpha1cur = MatrixXd(m_alpha1);
        m_alpha2cur = MatrixXd(m_alpha2);

        m_alpha1old = MatrixXd(m_alpha1);
        m_alpha2old = MatrixXd(m_alpha2);

        m_beta1cur = MatrixXd(m_beta1);
        m_beta2cur = MatrixXd(m_beta2);

        m_beta1old = MatrixXd(m_beta1);
        m_beta2old = MatrixXd(m_beta2);
    }

    // DESTRUCTOR
    gamPoisFactorEM::~gamPoisFactorEM() {}

    // member functions: documented in src

    /*!
     * \brief compute E-step (variational) in EM algo
     */
    void gamPoisFactorEM::Estep() {
        m_converged_Estep = false;
        // Iteration
        int nstab = 0; // number of successive iteration where the normalized gap betwwen two iteration is close to zero (convergence when nstab > rstab)
        int iter = 0;

        while( (iter < m_iterMax_Estep) && (m_converged_Estep==false)) {

            if(m_verbose==true) {
                Rcpp::Rcout << "iter " << iter << std::endl;
            }

            // Multinomial parameters
            //Rcpp::Rcout << "algorithm: Multinomial parameters" << std::endl;
            this->multinomParam();

            // local parameters
            // U : param phi
            //Rcpp::Rcout << "algorithm: local parameters" << std::endl;
            this->localParam();

            // expectation and log-expectation
            Egam(m_phi1cur, m_phi2cur, m_EU);
            Elgam(m_phi1cur, m_phi2cur, m_ElogU);

            // global parameters
            // V : param theta
            //Rcpp::Rcout << "algorithm: global parameters" << std::endl;
            this->globalParam();

            // expectation and log-expectation
            Egam(m_theta1cur, m_theta2cur, m_EV);
            Elgam(m_theta1cur, m_theta2cur, m_ElogV);

            // Poisson rate
            //Rcpp::Rcout << "algorithm: Poisson rate" << std::endl;
            this->poissonRate();

            // log-likelihood
            //Rcpp::Rcout << "algorithm: loglikelihood" << std::endl;
            this->computeLogLike(iter);
            // ELBO
            //Rcpp::Rcout << "algorithm: ELBO" << std::endl;
            this->computeELBO(iter);
            // deviance
            //Rcpp::Rcout << "algorithm: deviance" << std::endl;
            this->computeDeviance(iter);
            // explained variance
            //Rcpp::Rcout << "algorithm: explained variance" << std::endl;
            this->computeExpVar(iter);
            // convergence
            //Rcpp::Rcout << "algorithm: convergence ?" << std::endl;
            this->assessConvergenceEstep(iter, nstab);
            // increment values of parameters
            //Rcpp::Rcout << "algorithm: next iteration" << std::endl;
            this->nextIterate();
            // increment iteration
            iter++;
        }
    }

    /*!
     * \brief compute M-step (prior parameter updates) in EM algo
     */
    void gamPoisFactorEM::Mstep() {}

    /*!
     * \brief compute EM algorithm in gamma Poisson factor model
     */
    void gamPoisFactorEM::EMalgorithm() {}

    // create list with results to be return
    void gamPoisFactorEM::returnObject(Rcpp::List &results) {}

    //-------------------//
    // parameter updates //
    //-------------------//

    // local parameters: alpha (factor U)
    void gamPoisFactorEM::localPriorParam() {}

    // global parameters: beta (factor V)
    void gamPoisFactorEM::globalPriorParam() {}

    //-------------------//
    //     algorithm     //
    //-------------------//

    /*!
     * \brief assess convergence in E-step (variational)
     */
    void gamPoisFactorEM::assessConvergenceEstep(int iter, int &nstab) {
        // breaking condition: convergence or not
        double paramNorm = sqrt(parameterNorm2(m_phi1old, m_phi2old) + parameterNorm2(m_theta1old, m_theta2old));
        double diffNorm = sqrt(differenceNorm2(m_phi1old, m_phi2old, m_phi1cur, m_phi2cur) + differenceNorm2(m_theta1old, m_theta2old, m_theta1cur, m_theta2cur));

        m_normGap(iter) = diffNorm / paramNorm;

        // derivative order to consider
        double condition = convCondition(m_order, m_normGap_Estep, iter, 0);

        if(std::abs(condition) < m_epsilon) {
            nstab++;
        } else {
            nstab=0;
        }

        if(nstab > m_stabRange) {
            m_converged=true;
            m_nbIter=iter;
        } else {
            if(iter == m_iterMax - 1) {
                m_nbIter = iter;
            }
        }
    }

    void gamPoisFactorEM::assessConvergenceMstep(int iter, int &nstab) {}
    void gamPoisFactorEM::assessConvergence(int iter, int &nstab) {}


}
