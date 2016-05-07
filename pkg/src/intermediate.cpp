// Copyright 2016-04 Ghislain Durif
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
 * \file intermediate.cpp
 * \brief intermediate functions
 * \author Ghislain Durif
 * \version 0.1
 * \date 25/04/2016
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <cstdio>
#include "intermediate.h"

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;                  // variable size matrix, double precision

namespace intermediate {

    /*!
     * \fn Dirac function in zero
     *
     * Indicate if the considered value is null
     *
     * @param x value tested
     * @return 0 or 1 indicating if x==0 or not
     *
     * @TODO verify the 0-1 return
     */
    double dirac(double x) {

        if(x==0) {
            return 1;
        } else {
            return 0;
        }
    }

    /*!
     * \fn erase zeros in matrix
     *
     * @param[in,out] X matrix of data (count)
     */
    void eraseZero(MatrixXd &X) {
        int n = X.rows();
        int p = X.cols();
        for(int i=0; i<n; i++) {
            for(int j=0; j<n; j++) {
                if(X(i,j) == 0) {
                    X(i,j) = X(i,j) + 1;
                }
            }
        }
    }

    /*!
     * \fn exponential that throw error when under or overflowing
     *
     * @param[in] x real
     */
    double myExp(double x) {
        if(std::abs(x) > 700) {
            Rcpp::stop("under or overflow in exponential");
            return(-1);
        } else {
            return exp(x);
        }
    }

    /*!
     * \fn check if exponential can be applied to a matrix
     *
     * throw error when potential under or overflowing
     *
     * @param[in] A matrix to be checked
     */
    void checkExp(const MatrixXd &A) {
        // check for 300 and -300 so a product of such elements does not exceed the limit
        if( A.maxCoeff() > 300) {
            Rcpp::stop("overflow in element-wise exponentiation of matrix");
        }
        if( A.minCoeff() < - 300) {
            Rcpp::stop("underflow in element-wise exponentiation of matrix");
        }
    }

    /*!
     * \fn check matrix (dim, min and max element)
     *
     * @param[in] A matrix to be checked
     */
    void checkMat(const MatrixXd &A) {
        Rcpp::Rcout << "dimension: " << A.rows() << " x " << A.cols() << " / min = " << A.minCoeff() << " / max = " << A.maxCoeff() << std::endl;
    }
}
