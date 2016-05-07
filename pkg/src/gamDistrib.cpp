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
 * \brief class definition for Gamma Prior
 * \author Ghislain Durif
 * \version 0.1
 * \date 22/04/2016
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/math/special_functions/digamma.hpp>
#include "gamDistrib.h"

#define digamma() unaryExpr(std::ptr_fun<double,double>(digamma))
#define lgamma() unaryExpr(std::ptr_fun<double,double>(lgamma))
#define log() unaryExpr(std::ptr_fun<double,double>(log))

// [[Rcpp::depends(BH)]]
using boost::math::digamma;

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision

// SRC
namespace countMatrixFactor {

    /*!
     * \brief Compute expection of Gamma distribution
     *
     * E[U] = alpha/beta when alpha=param1 and beta=param2
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
    * \brief Compute expection of log of Gamma distribution
    *
    * E[log U] = digamma(alpha) - log(beta) when alpha=param1 and beta=param2
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
     * \brief Compute entropy of Gamma distribution
     *
     * E[-log p(U)] = (1-alpha)*digamma(alpha) + alpha - log(beta)
     *               + log(gamma(alpha))
     * when alpha=param1 and beta=param2
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

}