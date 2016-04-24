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
#define log() unaryExpr(std::ptr_fun<double,double>(log))
#define lgamma() unaryExpr(std::ptr_fun<double,double>(lgamma))

// [[Rcpp::depends(BH)]]
using boost::math::digamma;

// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision

// SRC
namespace countMatrixFactor {

    //------------------------------------------------------------------------//
    // CONSTRUCTOR
    gamDistrib::gamDistrib(int rows, int cols) : gamParam(rows, cols) {
        m_Egam = MatrixXd::Zero(rows,cols);
        m_Elgam = MatrixXd::Zero(rows,cols);
        m_entropy = MatrixXd::Zero(rows,cols);
    }

    gamDistrib::gamDistrib(int rows, int cols, const MatrixXd &param1, const MatrixXd &param2)
            : gamParam(rows, cols, param1, param2) {
        m_Egam = MatrixXd::Zero(rows,cols);
        m_Elgam = MatrixXd::Zero(rows,cols);
        m_entropy = MatrixXd::Zero(rows,cols);
    }

    // DESTRUCTOR
    gamDistrib::~gamDistrib() {}

    //------------------------------------------------------------------------//
    // GETTER
    /*!
     * \brief getter for expectation
     *
     * @param[out] Egam rows x cols, matrix of expectation
     */
    void gamDistrib::getE(MatrixXd &Egam) {
        Egam = m_Egam;
    }
    /*!
     * \brief getter for expectation of log
     *
     * @param[out] Elgam rows x cols, matrix of log-expectation
     */
    void gamDistrib::getElog(MatrixXd &Elgam) {
        Elgam = m_Elgam;
    }
    /*!
     * \brief getter for entropy
     *
     * @param[out] entropy rows x cols, matrix of entropy
     */
    void gamDistrib::getEntropy(MatrixXd &entropy) {
        entropy = m_entropy;
    }

    //------------------------------------------------------------------------//
    // member functions: documented in src

    // expectation
    /*!
     * \brief Compute expection of Gamma distribution
     *
     * E[U] = alpha/beta when alpha=param1 and beta=param2
     */
    void gamDistrib::E() {
        m_Egam = m_param1.array() / m_param2.array();
    }

    // log-expectation
    /*!
     * \brief Compute expection of log of Gamma distribution
     *
     * E[log U] = digamma(alpha) - log(beta) when alpha=param1 and beta=param2
     */
    void gamDistrib::Elog() {
        m_Elgam = m_param1.digamma().array() - m_param2.log().array();
    }

    // entropy
    /*!
     * \brief Compute entropy of Gamma distribution
     *
     * E[-log p(U)] = (1-alpha)*digamma(alpha) + alpha - log(beta)
     *                  + log(gamma(alpha))
     * when alpha=param1 and beta=param2
     */
    void gamDistrib::entropy() {
        m_entropy = (1-m_param1.array()) * m_param1.digamma().array()
                        + m_param1.array()
                        + m_param1.lgamma().array() - m_param2.log().array();
    }

}