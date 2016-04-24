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

#ifndef gamDistrib_H
#define gamDistrib_H

/*!
 * \file gamDistrib.h
 * \brief class definition for Gamma Prior
 * \author Ghislain Durif
 * \version 0.1
 * \date 22/04/2016
 */

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;                  // variable size matrix, double precision

namespace countMatrixFactor {
    /*!
     * \class gamPoisFactor
     * \brief class to define Gamma prior for latent factors
     */

    class gamDistrib {
    protected:
        // dimensions
        int m_rows;         /*!< number of rows (dimension of observation or variable space) */
        int m_cols;         /*!< number of factors */

        // hyper-parameters
        MatrixXd m_param1;  /*!< rows x cols, shape parameter (called alpha) */
        MatrixXd m_param2;  /*!< rows x cols, rate parameter (called beta) */

        // sufficient statistics
        MatrixXd m_Egam;    /*!< rows x cols, expectation, alpha/beta */
        MatrixXd m_Elgam;   /*!< rows x cols, log-expectation, digamma(alpha) - log(beta) */

        MatrixXd m_entropy; /*!< rows x cols, entropy of the Gamma distribution,
                                (1-alpha)*digamma(alpha) + alpha - log(beta) + log(gamma(alpha))*/

    public:
        /*!
         * \brief Constructor
         *
         * Constructor of the class gamDistrib without initialization
         */
        gamDistrib(int rows, int cols);

        /*!
         * \brief Constructor
         *
         * Constructor of the class gamDistrib with initialization
         */
        gamDistrib(int rows, int cols, const MatrixXd &param1, const MatrixXd &param2);

        /*!
         * \brief Destructor
         *
         * Destructor of the class gamDistrib
         */
        ~gamDistrib();

    public:
        // member functions: documented in src

        void expectation();
        void logexpectation();
        void entropy();
    };

}
