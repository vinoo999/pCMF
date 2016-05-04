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
 * \file gamPoisFactor.cpp
 * \brief class definition for Gamma Poisson Factor Model
 * \author Ghislain Durif
 * \version 0.1
 * \date 22/04/2016
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/math/special_functions/digamma.hpp>
#include "gamPoisFactor.h"

using namespace Rcpp;

#define coeffExp() unaryExpr(std::ptr_fun<double,double>(exp))
#define coeffSum(X) unaryExpr(std::bind2nd(std::pointer_to_binary_function<double,double,double>(scalsum),X))
#define square() unaryExpr(std::bind2nd(std::pointer_to_binary_function<double,double,double>(pow),2))
#define digamma() unaryExpr(std::ptr_fun<double,double>(digamma))
#define log() unaryExpr(std::ptr_fun<double,double>(log))

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
    gamPoisFactor::gamPoisFactor(int n, int p, int K, int iterMax, int order,
                                 int stabRange, double epsilon, bool verbose,
                                 const MatrixXi &X,
                                 const MatrixXd &phi1, const MatrixXd &phi2,
                                 const MatrixXd &theta1, const MatrixXd &theta2,
                                 const MatrixXd &alpha1, const MatrixXd &alpha2,
                                 const MatrixXd &beta1, const MatrixXd &beta2) :
                loglikelihood(iterMax), explainedVariance(iterMax) {

        // dimensions
        m_N = n;
        m_P = p;
        m_K = K;

        // parameters
        m_iterMax = iterMax;
        m_order = order;
        m_stabRange = stabRange;
        m_epsilon = epsilon;
        m_verbose = verbose;

        m_converged = false;
        m_nbIter = 0;

        // data
        m_X = MatrixXi(X);

        // variational parameters
        m_phi1cur = MatrixXd(phi1);
        m_phi2cur = MatrixXd(phi2);

        m_phi1old = MatrixXd(phi1);
        m_phi2old = MatrixXd(phi2);

        m_theta1cur = MatrixXd(theta1);
        m_theta2cur = MatrixXd(theta2);

        m_theta1old = MatrixXd(theta1);
        m_theta2old = MatrixXd(theta2);

        // sufficient statistics
        m_EU = MatrixXd::Zero(n,K);
        Egam(phi1, phi2, m_EU);
        m_ElogU = MatrixXd::Zero(n,K);
        Elgam(phi1, phi2, m_ElogU);

        m_EV = MatrixXd::Zero(p,K);
        Egam(theta1, theta2, m_EV);
        m_ElogV = MatrixXd::Zero(p,K);
        Elgam(theta1, theta2, m_ElogV);

        m_EZ_i = MatrixXd::Zero(p,K);
        m_EZ_j = MatrixXd::Zero(n,K);

        // prior parameter
        m_alpha1 = MatrixXd(alpha1);
        m_alpha2 = MatrixXd(alpha1);
        m_beta1 = MatrixXd(beta1);
        m_beta2 = MatrixXd(beta2);

        // criterion
        m_normGap = VectorXd::Zero(iterMax);

    }

    // DESTRUCTOR
    gamPoisFactor::~gamPoisFactor() {}

    //-------------------//
    //   convergence     //
    //-------------------//

    /*!
     * \fn squared norm of parameters
     *
     * @param[in] param1 matrix of parameters 1
     * @param[in] param2 matrix of parameters 2
     *
     * @return sum(param1^2) + sum(param2^2)
     */
    double parameterNorm2(const MatrixXd &param1, const MatrixXd &param2) {
        double res = param1.square().sum() + param2.square().sum();
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