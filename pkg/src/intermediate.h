// Copyright 2016-04 Ghislain Durif
//
// This file is part of the `pCMF' library for R and related languages.
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

#ifndef intermediate_H
#define intermediate_H

/*!
 * \file intermediate.h
 * \brief intermediate functions
 * \author Ghislain Durif
 * \version 0.1
 * \date 25/04/2016
 */

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;                  // variable size matrix, double precision


/*!
 * \namespace intermediate
 *
 * A specific namespace for internal functions
 */
namespace intermediate {

    // Dirac function in zero
    double dirac(double x);

    // erase zero
    void eraseZero(MatrixXd &X);

    // threshold at t
    double threshold(double x, double t);

    // exponential that throw error when under or overflowing
    double myExp(double x);

    // exponential that throw error when under or overflowing
    double safeExp(double x);

    // logit function
    double logit(double x);

    // logit inverse function
    double expit(double x);

    // compute E[log(Z!)] where Z follows a binomial distribution
    double lgamBinom(int N, double p);

    // check if exponential can be applied to a matrix
    void checkExp(const MatrixXd &A);

    // check for overflow in exponential
    int checkMaxExp(const MatrixXd &A);

    // check for underflow in exponential
    int checkMinExp(const MatrixXd &A);

    // check matrix (dim, min and max element)
    void checkMat(const MatrixXd &A);

    // inverse psi (digamma function)
    double psiInv(double y, int nbIter);

    // parameter squared euclidean norm
    double parameterNorm2(const MatrixXd &param1, const MatrixXd &param2);

    // difference of squared euclidean norm (on parameters)
    double differenceNorm2(const MatrixXd &param1a, const MatrixXd &param2a, const MatrixXd &param1b, const MatrixXd &param2b);

}

#endif
