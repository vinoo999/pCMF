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
* \file gamPoisFactor_wrapper.cpp
* \brief definition of wrapper (to call in R) for standard Gamma Poisson Factor Model
* \author Ghislain Durif
* \version 0.1
* \date 06/05/2016
*/

#include <Rcpp.h>
#include <RcppEigen.h>
#include <cstdio>
#include "internal/variational.h"
#include "internal/model/gamPoisFactorPen.h"

using namespace countMatrixFactor;


// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision


// [[Rcpp::export]]
SEXP gamPoisFactor_wrapper(SEXP Xin, int K,
                           SEXP phi01in, SEXP phi02in,
                           SEXP theta01in, SEXP theta02in,
                           SEXP alpha1in, SEXP alpha2in,
                           SEXP beta1in, SEXP beta2in,
                           SEXP lambda_in, SEXP mu_in,
                           int iterMax, double epsilon,
                           int order, int stabRange, bool verbose) {

    MatrixXi X = Rcpp::as< Map<MatrixXi> >(Xin);
    MatrixXd phi01 = Rcpp::as< Map<MatrixXd> >(phi01in);
    MatrixXd phi02 = Rcpp::as< Map<MatrixXd> >(phi02in);
    MatrixXd theta01 = Rcpp::as< Map<MatrixXd> >(theta01in);
    MatrixXd theta02 = Rcpp::as< Map<MatrixXd> >(theta02in);
    MatrixXd alpha1 = Rcpp::as< Map<MatrixXd> >(alpha1in);
    MatrixXd alpha2 = Rcpp::as< Map<MatrixXd> >(alpha2in);
    MatrixXd beta1 = Rcpp::as< Map<MatrixXd> >(beta1in);
    MatrixXd beta2 = Rcpp::as< Map<MatrixXd> >(beta2in);

    VectorXd lambda = Rcpp::as< Map<VectorXd> >(lambda_in);
    VectorXd mu = Rcpp::as< Map<VectorXd> >(mu_in);

    int n = X.rows();
    int p = X.cols();

    // declaration of object gamPoisFactorStandard
    Rcpp::Rcout << "Declaration" << std::endl;
    variational<gamPoisFactorPen> myModel(iterMax, order,
                                          stabRange, epsilon, verbose,
                                          n, p, K, X,
                                          phi01, phi02, theta01, theta02,
                                          alpha1, alpha2, beta1, beta2,
                                          lambda, mu);

    // initialization
    Rcpp::Rcout << "Initialization" << std::endl;
    myModel.Init();

    // computations
    Rcpp::Rcout << "Algorithm" << std::endl;
    myModel.algorithm();

    // factor order
    Rcpp::Rcout << "factor order" << std::endl;
    myModel.computeOrder();

    // returns
    Rcpp::Rcout << "Output" << std::endl;
    Rcpp::List results;
    myModel.returnObject(results);

    return results;

}