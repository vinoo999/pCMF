// Copyright 2016-04 Ghislain Durif
//
// This file is part of the `cmf' library for R and related languages.
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
* \brief functions for log-likelihood computation
* \author Ghislain Durif
* \version 0.1
* \date 22/04/2016
*/

#include <Rcpp.h>
#include <RcppEigen.h>
#include "loglikelihood.h"
#include "intermediate.h"

#define mlog() unaryExpr(std::ptr_fun<double,double>(log))
#define mlgamma() unaryExpr(std::ptr_fun<double,double>(lgamma))
#define mdirac() unaryExpr(std::ptr_fun<double,double>(intermediate::dirac))

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;              // variable size matrix, double precision
using Eigen::MatrixXi;              // variable size matrix, integer
using Eigen::VectorXd;              // variable size matrix, double precision

namespace countMatrixFactor {

    //------------------------------------------------------------------------//
    // FUNCTIONS
    /*!
     * \fn log-likelihood for Gamma distribution
     *
     * @param[in] X matrix of observations
     * @param[in] param1 matrix of corresponding first parameters from Gamma distribution
     * @param[in] param2 matrix of corresponding second parameters from Gamma distribution
     *
     * @return effective value
     */
    double gammaLogLike(const MatrixXd &X, const MatrixXd &param1, const MatrixXd &param2) {
        double res;
        // sum( (param1.ij - 1)*log(X.ij)  + param1.ij*log(param2.ij) - param2.ij * X.ij - lgamma(param1.ij) )
        res = ( param1.array() * param2.mlog().array() + (param1.array() - 1) * X.cast<double>().mlog().array() - param2.array() * X.cast<double>().array() - param1.mlgamma().array() ).sum();
        return res;
    }

    /*!
     * \fn log-likelihood for Poisson distribution
     *
     * @param[in] X matrix of observations
     * @param[in] lambda matrix of corresponding Poisson rates
     *
     * @return effective value
     */
    double poisLogLike(const MatrixXi &X, const MatrixXd &lambda) {
        double res;
        // sum(X.ij * log(lambda.ij) - lambda.ij - lfactorial(X.ij))
        res = ( X.cast<double>().array() * lambda.mlog().array() - lambda.array() - (X.cast<double>().array() + 1).mlgamma() ).sum();
        return(res);
    }

    /*!
     * \fn log-likelihood for Zero-Inflated Poisson distribution
     *
     * Mixture between a Bernoulli and Poisson distribution (with proba pi)
     *
     * @param[in] X matrix of observations
     * @param[in] lambda matrix of corresponding Poisson rates
     * @param[in] pi matrix of zero-inflation probabilities
     *
     * @return effective value
     */
    double ZIpoisLogLike(const MatrixXi &X, const MatrixXd &lambda, const MatrixXd &pi) {
        double res1(0);
        double res2(0);
        int n(X.rows());
        int p(X.cols());
        // sum( log( (1-pi.ij) * delta_0(X.ij) + pi.ij * Poisson(X.ij, lambda.ij) ))
        // sum( (1 - delta_0(X.ij)) * (X.ij * log(lambda.ij) - lambda.ij - lfactorial(X.ij)) + delta_0(X.ij) * log( (1-pi.ij) + pi.ij * Poisson(X.ij, lambda.ij) ) )
        res1 = ((1 - X.cast<double>().mdirac().array()) * ( X.cast<double>().array() * lambda.mlog().array() - lambda.array() - (X.array() + 1).mlgamma() )).sum();
        res2 = (( X.cast<double>().mdirac().array()) * ( 1 - pi.array() + pi.array() * (X.cast<double>().array() * lambda.mlog().array() - lambda.array() - (X.cast<double>().array() + 1).mlgamma().array()).exp() )).mlog().sum();
        //std::cout << X.cast<double>().unaryExpr(std::pointer_to_unary_function<double,double>(dirac)) << std::endl; //.array().sum;
        return res1+res2;
    }

    /*!
     * \fn log-likelihood for Poisson saturated model, i.e. with lambda = X
     *
     * @param[in] X matrix of observations
     * @param[in] lambda matrix of corresponding Poisson rates,
     * supposed to be X in that case (where null values are replaced by 1,
     * for compatibility with log, convention 0 * log 0 = 0)
     *
     * @return effective value
     */
    double poisLogLikeSaturated(const MatrixXi &X, const MatrixXd &lambda0) {
        double res;
        // sum(X.ij * log(lambda0.ij) - X.ij - lfactorial(X.ij))
        // lambda0 = X (where 0 are replaced by 1 for convention 0 * log 0 = 0)
        res = ( X.cast<double>().array() * lambda0.mlog().array() - X.cast<double>().array() - (X.cast<double>().array() + 1).mlgamma() ).sum();
        return(res);
    }

    /*!
     * \fn deviance between estimated Poisson and saturated Poisson models
     *
     * @param[in] X matrix of observations
     * @param[in] lambda matrix of corresponding Poisson rates,
     * @param[in] lambda0 matrix of observations, i.e. X (where null values are replaced by 1,
     * for compatibility with log)
     *
     * @return effective value
     */
    double poisDeviance(const MatrixXi &X, const MatrixXd &lambda, const MatrixXd &lambda0) {
        MatrixXd lambda_tmp = lambda;
        intermediate::eraseZero(lambda_tmp);
        double res1 = ( X.cast<double>().array() * lambda_tmp.mlog().array() - lambda.array() ).sum();
        double res2 = ( X.cast<double>().array() * lambda0.mlog().array() - X.cast<double>().array() ).sum();
        double res = -2*(res1 - res2);
        return(res);
    }
}










