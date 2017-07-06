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

/*!
* \file wrapper_variational_ZI_GaP.cpp
* \brief definition of wrapper (to call in R) for zero-inflated Gamma Poisson Factor Model
* \author Ghislain Durif
* \version 0.1
* \date 09/06/2016
*/

#include <Rcpp.h>
#include <RcppEigen.h>
#include <cstdio>
#include "variational.h"
#include "gamPoisFactorZI.h"
#include "wrapper_template.h"

#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace countMatrixFactor;


// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision

//' @title wrapper_variational_ZI_GaP
//' @keywords internal
//'
//' @description
//' Description
//'
//' @details
//' Wrapper for Cpp function
//'
//' @author
//' Ghislain Durif, \email{gd.dev@libertymail.net}
//'
//'
//' @seealso aaa
//'
//' @import Rcpp
//' @import RcppEigen
//' @useDynLib pCMF
//'
//' @param X n x p, count data matrix
//' @param K dimension of the latent subspace
//' @param ZI boolean, indicating if the model is zero-inflated
//' @param phi1 n x K, initial values of first parameter of Gamma distribution on U
//' @param phi2 n x K, initial values of second parameter of Gamma distribution on U
//' @param theta1 n x K, initial values of first parameter of Gamma distribution on V
//' @param theta2 n x K, initial values of second parameter of Gamma distribution on V
//' @param alpha1 n x K, initial values of first parameter of Gamma prior on U
//' @param alpha2 n x K, initial values of second parameter of Gamma prior on U
//' @param beta1 n x K, initial values of first parameter of Gamma prior on V
//' @param beta2 n x K, initial values of second parameter of Gamma prior on V
//' @param iterMax maximum number of iterations
//' @param iterMin minimum number of iterations
//' @param epsilon precision for comparison when assessing convergence
//' @param order derivative order on normalized gap to assess convergence
//' (0 is the current value, 1 the first order empirical derivative, 2 the second order empirical derivative)
//' @param stabRange range of stability (number of iterations where parameter values are stable to confirm convergence)
//' @param verbose boolean indicating verbosity in the output
//' @param ncores number of cores to use for multi-threading
//' @param nbInit number of initialization to try
//' @param iterMaxInit number of iteration to for each initialization in multi-initialization case (nbInit > 1)
//' @param noise nosie level (between 0 and 1) for multi-init case
//' @param seed for RNG (-1 means no seed)
//'
//' @return return
//' \item{Y}{Y}
//'
//' @export
// [[Rcpp::export]]
SEXP wrapper_variational_ZI_GaP(SEXP Xin, int K, bool ZI,
                                SEXP phi01in, SEXP phi02in,
                                SEXP theta01in, SEXP theta02in,
                                SEXP alpha1in, SEXP alpha2in,
                                SEXP beta1in, SEXP beta2in,
                                int iterMax, int iterMin, double epsilon,
                                int order, int stabRange, bool verbose, int ncores,
                                int nbInit, int iterMaxInit, double noise, int seed) {

    MatrixXi X = Rcpp::as< Map<MatrixXi> >(Xin);
    MatrixXd phi01 = Rcpp::as< Map<MatrixXd> >(phi01in);
    MatrixXd phi02 = Rcpp::as< Map<MatrixXd> >(phi02in);
    MatrixXd theta01 = Rcpp::as< Map<MatrixXd> >(theta01in);
    MatrixXd theta02 = Rcpp::as< Map<MatrixXd> >(theta02in);
    MatrixXd alpha1 = Rcpp::as< Map<MatrixXd> >(alpha1in);
    MatrixXd alpha2 = Rcpp::as< Map<MatrixXd> >(alpha2in);
    MatrixXd beta1 = Rcpp::as< Map<MatrixXd> >(beta1in);
    MatrixXd beta2 = Rcpp::as< Map<MatrixXd> >(beta2in);

    return wrapper_template<gamPoisFactorZI, variational>(X, K, ZI,
                                                          phi01, phi02,
                                                          theta01, theta02,
                                                          alpha1, alpha2,
                                                          beta1, beta2,
                                                          iterMax, iterMin, epsilon,
                                                          order, stabRange, verbose, ncores,
                                                          nbInit, iterMaxInit, noise, seed);

}
