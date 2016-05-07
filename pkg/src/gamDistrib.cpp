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
 * \file gamDistrib.cpp
 * \fn class definition for Gamma Prior
 * \author Ghislain Durif
 * \version 0.1
 * \date 22/04/2016
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <boost/math/special_functions/digamma.hpp>
#include "gamDistrib.h"

#define digamma() unaryExpr(std::ptr_fun<double,double>(digamma))
#define lgamma() unaryExpr(std::ptr_fun<double,double>(lgamma))
#define log() unaryExpr(std::ptr_fun<double,double>(std::log))
#define square() unaryExpr(std::bind2nd(std::pointer_to_binary_function<double,double,double>(std::pow),2))

// [[Rcpp::depends(BH)]]
using boost::math::digamma;

// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision

// SRC
namespace countMatrixFactor {

    /*!
     * \fn Compute expection of Gamma distribution
     *
     * E[U] = alpha/beta when alpha=param1 (shape) and beta=param2 (rate)
     *
     * @param[in] param1 rows x cols, matrix of first parameters
     * @param[in] param2 rows x cols, matrix of second parameters
     * @param[out] res rows x cols, matrix of expectation for each couple
     * of parameters
     */
    void Egam(const MatrixXd &param1, const MatrixXd &param2, MatrixXd &res) {
        res = param1.array() / param2.array();
    }

    /*!
    * \fn Compute expection of log of Gamma distribution
    *
    * E[log U] = digamma(alpha) - log(beta)
    * when alpha=param1 (shape) and beta=param2 (rate)
    *
    * @param[in] param1 rows x cols, matrix of first parameters
    * @param[in] param2 rows x cols, matrix of second parameters
    * @param[out] res rows x cols, matrix of expectation of log for each couple
    * of parameters
    */
    void Elgam(const MatrixXd &param1, const MatrixXd &param2, MatrixXd &res) {
        res = param1.digamma().array() - param2.log().array();
    }

    /*!
     * \fn Compute entropy of Gamma distribution
     *
     * E[-log p(U)] = (1-alpha)*digamma(alpha) + alpha - log(beta)
     *               + log(gamma(alpha))
     * when alpha=param1 (shape) and beta=param2 (rate)
     *
     * @param[in] param1 rows x cols, matrix of first parameters
     * @param[in] param2 rows x cols, matrix of second parameters
     * @param[out] res rows x cols, matrix of entropy for each couple
     * of parameters
     */
    void entropyGam(const MatrixXd &param1, const MatrixXd &param2, MatrixXd &res) {
        res = (1-param1.array()) * param1.digamma().array() + param1.array()
        + param1.lgamma().array() - param2.log().array();
    }


    /*!
     * \fn Estimate shape and rate parameters
     *
     * Simulate n repetition of Gamma(param1, param2) and compute estimation of param1 and param2
     *
     * Note: param1 = shape, param2 = rate
     *
     * Estimations: param1 = mean^2 / variance and param2 = mean/variance
     *
     * @param[in] n size of the sample
     * @param[in] param1 shape parameter
     * @param[in] param2 rate parameter
     * @param[out] param1e estimation of shape parameter
     * @param[out] param2e estimation of rate parameter
     */
    void estimParam(double n, double param1, double param2, double &param1e, double &param2e) {
        Rcpp::NumericVector sample = Rcpp::rgamma(n, param1, param2);
        VectorXd sample2 = Rcpp::as<Map<VectorXd> >(sample);
        double mean = sample2.mean();
        double var = sample2.square().mean() - std::pow(mean,2);
        param1e = std::pow(mean,2) / var;
        param2e = mean / var;
    }

}