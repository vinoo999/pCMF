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
* \file gamParam.cpp
* \brief class definition for Gamma Prior
* \author Ghislain Durif
* \version 0.1
* \date 22/04/2016
*/

#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <boost/math/special_functions/digamma.hpp>
#include "gamParam.h"

#define square() unaryExpr(std::bind2nd(std::pointer_to_binary_function<double,double,double>(pow),2))
#define digamma() unaryExpr(std::ptr_fun<double,double>(digamma))
#define log() unaryExpr(std::ptr_fun<double,double>(log))
#define lgamma() unaryExpr(std::ptr_fun<double,double>(lgamma))

// [[Rcpp::depends(BH)]]
using boost::math::digamma;

// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                   // 'maps' rather than copies
using Eigen::MatrixXd;              // variable size matrix, double precision
using Eigen::MatrixXi;              // variable size matrix, integer
using Eigen::VectorXd;              // variable size vector, double precision

// SRC
namespace countMatrixFactor {

    //------------------------------------------------------------------------//
    // CONSTRUCTOR
    gamParam::gamParam(int rows, int cols) {
        // dimensions
        m_rows = rows;
        m_cols = cols;

        // hyper-parameters
        m_param1 = MatrixXd::Zero(rows,cols);
        m_param2 = MatrixXd::Zero(rows,cols);
    }

    gamParam::gamParam(int rows, int cols,
                       const MatrixXd &param1, const MatrixXd &param2) {
        // dimensions
        m_rows = rows;
        m_cols = cols;

        // hyper-parameters
        m_param1 = MatrixXd(param1);
        m_param2 = MatrixXd(param2);
    }

    // DESTRUCTOR
    gamParam::~gamParam() {}

    //------------------------------------------------------------------------//
    // GETTER
    /*!
    * \brief getter for param1
    *
    * @param[out] rows x cols, matrix of param1
    */
    void gamParam::getParam1(MatrixXd &param1) {
        param1 = m_param1;
    }
    /*!
    * \brief getter for param2
    *
    * @param[out] rows x cols, matrix of param2
    */
    void gamParam::getParam2(MatrixXd &param2) {
        param2 = m_param2;
    }

    // SETTER
    /*!
    * \brief setter for param1
    *
    * @param[in] rows x cols, matrix of param1
    */
    void gamParam::setParam1(MatrixXd &param1) {
        m_param1 = param1;
    }
    /*!
    * \brief setter for param2
    *
    * @param[in] rows x cols, matrix of param2
    */
    void gamParam::setParam2(MatrixXd &param2) {
        m_param2 = param2;
    }

    //------------------------------------------------------------------------//
    // MEMBER FUNCTIONS

    // expectation
    /*!
     * \brief Compute expection of Gamma distribution
     *
     * E[U] = alpha/beta when alpha=param1 and beta=param2
     *
     * @param[out] res rows x cols, matrix of expectation for each couple
     * of parameters
     */
    void gamParam::E(MatrixXd &res) {
        res = m_param1.array() / m_param2.array();
    }

    // log-expectation
    /*!
     * \brief Compute expection of log of Gamma distribution
     *
     * E[log U] = digamma(alpha) - log(beta) when alpha=param1 and beta=param2
     *
     * @param[out] res rows x cols, matrix of expectation of log
     * for each couple of parameters
     */
    void gamParam::Elog(MatrixXd &res) {
        res = m_param1.digamma().array() - m_param2.log().array();
    }

    // entropy
    /*!
     * \brief Compute entropy of Gamma distribution
     *
     * E[-log p(U)] = (1-alpha)*digamma(alpha) + alpha - log(beta)
     *                  + log(gamma(alpha))
     * when alpha=param1 and beta=param2
     *
     * @param[out] res rows x cols, matrix of entropy for each couple
     * of parameters
     */
    void gamParam::entropy(MatrixXd &res) {
        res = (1-m_param1.array()) * m_param1.digamma().array()
                + m_param1.array()
                + m_param1.lgamma().array() - m_param2.log().array();
    }

    // parameter norm
    /*!
     * \brief l2 squared norm of all parameters
     *
     * Computation of sum_{ij} param1_{ij}^2 + param2_{ij}^2
     *
     * @return res l2 squared norm
     */
    double gamParam::parameterNorm2() {
        double res = m_param1.square().sum() + m_param2.square().sum();
        return res;
    }

    //------------------------------------------------------------------------//
    // FUNCTIONS

    // expectation
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

    // log-expectation
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

    // entropy
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

    // parameter norm
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
        double res = param1.square().sum() + param2.square().sum();
        return res;
    }

}