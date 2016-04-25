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
* \file loglikelihood.cpp
* \brief class definition for log-likelihood
* \author Ghislain Durif
* \version 0.1
* \date 22/04/2016
*/

#include <Rcpp.h>
#include <RcppEigen.h>
#include "loglikelihood.h"
#include "intermediate.h"

#define digamma() unaryExpr(std::ptr_fun<double,double>(digamma))
#define log() unaryExpr(std::ptr_fun<double,double>(log))
#define lgamma() unaryExpr(std::ptr_fun<double,double>(lgamma))
#define dirac() unaryExpr(std::ptr_fun<double,double>(intermediate::dirac))

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;              // variable size matrix, double precision
using Eigen::MatrixXi;              // variable size matrix, integer
using Eigen::VectorXd;              // variable size matrix, double precision

namespace countMatrixFactor {

    //------------------------------------------------------------------------//
    // CONSTRUCTOR
    loglikelihood::loglikelihood(int size) {
        m_margLogLike = VectorXd::Zero(size);
        m_condLogLike = VectorXd::Zero(size);
        m_priorLogLike = VectorXd::Zero(size);
        m_postLogLike = VectorXd::Zero(size);
        m_compLogLike = VectorXd::Zero(size);
    }

    // DESTRUCTOR
    loglikelihood::~loglikelihood() {}

    // getter
    /*!
     * \brief getter for marginal log-likelihood
     *
     * @param[out] res vector of marginal log-likelihood
     */
    void loglikelihood::getMarginal(VectorXd &res) {
        res = m_margLogLike;
    }
    /*!
     * \brief getter for conditional log-likelihood
     *
     * @param[out] res vector of conditional log-likelihood
     */
    void loglikelihood::getConditional(VectorXd &res) {
        res = m_condLogLike;
    }
    /*!
     * \brief getter for prior log-likelihood
     *
     * @param[out] res vector of prior log-likelihood
     */
    void loglikelihood::getPrior(VectorXd &res) {
        res = m_priorLogLike;
    }
    /*!
     * \brief getter for posterior log-likelihood
     *
     * @param[out] res vector of posterior log-likelihood
     */
    void loglikelihood::getPosterior(VectorXd &res) {
        res = m_postLogLike;
    }

    /*!
     * \brief getter for complete log-likelihood
     *
     * @param[out] res vector of complete log-likelihood
     */
    void loglikelihood::getComplete(VectorXd &res) {
        res = m_compLogLike;
    }


    // FUNCTIONS
    // local log-likelihood function

    double gammaLogLike(const MatrixXd &X, const MatrixXd &alpha, const MatrixXd &beta) {
        double res;
        // sum( (alpha.ij - 1)*log(X.ij)  + alpha.ij*log(beta.ij) - beta.ij * X.ij - lgamma(alpha.ij) )
        res = ( alpha.array() * beta.log().array() + (alpha.array() - 1) * X.cast<double>().log().array() - beta.array() * X.array() - alpha.lgamma().array() ).sum();
        return res;
    }

    double poisLogLike(const MatrixXi &X, const MatrixXd &lambda) {
        double res;
        // sum(X.ij * log(lambda.ij) - lambda.ij - lfactorial(X.ij))
        res = ( X.cast<double>().array() * lambda.log().array() - lambda.array() - (X.array() + 1).lgamma() ).sum();
        return(res);
    }

    double ZIpoisLogLike(const MatrixXi &X, const MatrixXd &lambda, const MatrixXd &pi) {
        double res1(0);
        double res2(0);
        int n(X.rows());
        int p(X.cols());
        // sum( log( (1-pi.ij) * delta_0(X.ij) + pi.ij * Poisson(X.ij, lambda.ij) ))
        // sum( (1 - delta_0(X.ij)) * (X.ij * log(lambda.ij) - lambda.ij - lfactorial(X.ij)) + delta_0(X.ij) * log( (1-pi.ij) + pi.ij * Poisson(X.ij, lambda.ij) ) )
        res1 = ((1 - X.cast<double>().dirac().array()) * ( X.cast<double>().array() * lambda.log().array() - lambda.array() - (X.array() + 1).lgamma() )).sum();
        res2 = (( X.cast<double>().dirac().array()) * ( 1 - pi.array() + pi.array() * (X.cast<double>().array() * lambda.log().array() - lambda.array() - (X.array() + 1).lgamma().array()).exp() )).log().sum();
        //std::cout << X.cast<double>().unaryExpr(std::pointer_to_unary_function<double,double>(dirac)) << std::endl; //.array().sum;
        return res1+res2;
    }
}










