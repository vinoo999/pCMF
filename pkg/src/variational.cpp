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
* \file variational.cpp
* \brief tamplete class definition for variational algorithm
* \author Ghislain Durif
* \version 0.1
* \date 06/06/2016
*/

#include <Rcpp.h>
#include <RcppEigen.h>
#include "variational.h"

#define msquare() unaryExpr(std::bind2nd(std::pointer_to_binary_function<double,double,double>(std::pow),2))

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision

namespace countMatrixFactor {

    // CONSTRUCTOR
    template <typename model>
    variational<model>::variational(int iterMax, int order,
                int stabRange, double epsilon, bool verbose,
                int n, int p, int K,
                const MatrixXi &X,
                const MatrixXd &phi1, const MatrixXd &phi2,
                const MatrixXd &theta1, const MatrixXd &theta2,
                const MatrixXd &alpha1, const MatrixXd &alpha2,
                const MatrixXd &beta1, const MatrixXd &beta2)
    : loglikelihood(iterMax),
      explainedVariance(iterMax),
      m_model(n, p, K,
              X, phi1, phi2,
              theta1, theta2,
              alpha1, alpha2,
              beta1, beta2)
    {

        // parameters
        m_iterMax = iterMax;
        m_iter = 0;
        m_order = order;
        m_stabRange = stabRange;
        m_epsilon = epsilon;
        m_verbose = verbose;

        m_converged = false;
        m_nbIter = 0;

        // criterion
        m_normGap = VectorXd::Zero(iterMax);
    }

    // DESTRUCTOR
    template <typename model>
    variational<model>::~variational() {}

    /*!
     * \brief Initialization of sufficient statistics
     */
    template <typename model>
    void variational<model>::Init() {}

    /*!
     * \brief run algorithm
     */
    template <typename model>
    void variational<model>::algorithm() {

        // Iteration
        int nstab = 0; // number of successive iteration where the normalized gap betwwen two iteration is close to zero (convergence when nstab > rstab)

        while( (m_iter < m_iterMax) && (m_converged==false)) {

            if(m_verbose==true) {
                Rcpp::Rcout << "iter " << m_iter << std::endl;
            }

            // Multinomial parameters
            //Rcpp::Rcout << "algorithm: Multinomial parameters" << std::endl;
            m_model.multinomParam();

            // local parameters
            // U : param phi
            //Rcpp::Rcout << "algorithm: local parameters" << std::endl;
            m_model.localParam();

            // global parameters
            // V : param theta
            //Rcpp::Rcout << "algorithm: global parameters" << std::endl;
            m_model.globalParam();

            // Poisson rate
            //Rcpp::Rcout << "algorithm: Poisson rate" << std::endl;
            m_model.poissonRate();

            // log-likelihood
            //Rcpp::Rcout << "algorithm: loglikelihood" << std::endl;
            this->computeLogLike(m_iter);
            // ELBO
            //Rcpp::Rcout << "algorithm: ELBO" << std::endl;
            this->computeELBO(m_iter);
            // deviance
            //Rcpp::Rcout << "algorithm: deviance" << std::endl;
            this->computeDeviance(m_iter);
            // explained variance
            //Rcpp::Rcout << "algorithm: explained variance" << std::endl;
            this->computeExpVar(m_iter);
            // convergence
            //Rcpp::Rcout << "algorithm: convergence ?" << std::endl;
            this->assessConvergence(nstab);
            // increment values of parameters
            //Rcpp::Rcout << "algorithm: next iteration" << std::endl;
            m_model.nextIterate();
            // increment iteration
            m_iter++;
        }

    }

    /*!
     * \brief create list of object to return
     */
    template <typename model>
    void variational<model>::returnObject(Rcpp::List &results) {}

    /*!
     * \brief assess convergence
     *
     * @param[in] iter current iteration
     * @param[in,out] nstab number of successive iteration respecting the breaking condition
     */
    template <typename model>
    void variational<model>::assessConvergence(int &nstab) {
        // breaking condition: convergence or not
        double res = m_model.normGap();
        m_normGap(m_iter) = res;

        // derivative order to consider
        double condition = convCondition(m_order, m_normGap, m_iter, 0);

        if(std::abs(condition) < m_epsilon) {
            nstab++;
        } else {
            nstab=0;
        }

        if(nstab > m_stabRange) {
            m_converged=true;
            m_nbIter=m_iter;
        } else {
            if(m_iter == m_iterMax - 1) {
                m_nbIter = m_iter;
            }
        }
    }

    //-------------------//
    //      criteria     //
    //-------------------//

    /*!
    * \brief compute log-likelihood
    *
    * @param[in] iter current iteration
    */
    template <typename model>
    void variational<model>::computeLogLike(int iter) {
        m_condLogLike(iter) = m_model.computeCondLogLike();
        m_priorLogLike(iter) = m_model.computePriorLogLike();
        m_postLogLike(iter) = m_model.computePostLogLike();
        m_compLogLike(iter) = m_condLogLike(iter) + m_postLogLike(iter);
        m_margLogLike(iter) = m_condLogLike(iter) + m_priorLogLike(iter) - m_postLogLike(iter);
    }

    /*!
    * \brief compute evidence lower bound
    *
    * @param[in] iter current iteration
    */
    template <typename model>
    void variational<model>::computeELBO(int iter) {
        m_elbo(iter) = m_model.computeELBO();
    }

    /*!
     * \brief compute deviance between estimated and saturated model
     *
     * @param[in] iter current iteration
     */
    template <typename model>
    void variational<model>::computeDeviance(int iter) {
        m_deviance(iter) = m_model.computeDeviance();
    }

    /*!
     * \brief compute explained variance
     *
     * @param[in] iter current iteration
     */
    template <typename model>
    void variational<model>::computeExpVar(int iter) {
        m_expVar0(iter) = m_model.computeExpVar0();
        m_expVarU(iter) = m_model.computeExpVarU();
        m_expVarV(iter) = m_model.computeExpVarV();
    }


    //-------------------//
    //   convergence     //
    //-------------------//

    /*!
    * \brief l2 squared norm of all parameters
    *
    * Computation of sum_{ij} param1_{ij}^2 + param2_{ij}^2
    *
    * @param[in] param1 rows x cols, matrix of first parameters
    * @param[in] param2 rows x cols, matrix of second parameters
    * @return res l2 squared norm
    */
    double parameterNorm2(const MatrixXd &param1, const MatrixXd &param2) {
        double res = param1.msquare().sum() + param2.msquare().sum();
        return res;
    }

    /*!
    * \fn difference of squared euclidean norm (on parameters)
    *
    * @param[in] param1a matrix of parameters 1 (state a)
    * @param[in] param2a matrix of parameters 2 (state a)
    * @param[in] param1b matrix of parameters 1 (state b)
    * @param[in] param2b matrix of parameters 2 (state b)
    *
    * @return sum((param1a - param1b)^2) + sum((param2a-param2b)^2)
    */
    double differenceNorm2(const MatrixXd &param1a, const MatrixXd &param2a, const MatrixXd &param1b, const MatrixXd &param2b) {
        double res = parameterNorm2(param1a - param1b, param2a - param2b);
        return(res);
    }

    /*!
     * \fn assess convergence condition
     *
     * convergence assessed on the normalized gap between two iterates
     * order 0: value of normalized gap
     * order 1: first empirical derivative of normalized gap (speed)
     * order 2: second empirical derivative of normalized gap (acceleration)
     *
     * @return value of the condition
     */
    double convCondition(int order, const VectorXd &normGap, int iter, int drift) {
        double condition = 1;
        switch(order) {
            case 0 : {
                condition = normGap(iter-drift);
            } break;

            case 1 : {
                if(iter>1) {
                    condition = normGap(iter-drift) - normGap(iter-drift-1);
                }
            } break;

            case 2 : {
                if(iter>2) {
                    condition = normGap(iter-drift) - 2*normGap(iter-drift-1) + normGap(iter-drift-2);
                }
            } break;
        }
        return(condition);
    }

}
