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

#ifndef loglikelihood_H
#define loglikelihood_H

/*!
* \file loglikelihood.h
* \brief functions for log-likelihood computation
* \author Ghislain Durif
* \version 0.1
* \date 22/04/2016
*/

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;              // variable size matrix, double precision
using Eigen::MatrixXi;              // variable size matrix, integer
using Eigen::VectorXd;              // variable size matrix, double precision

namespace countMatrixFactor {

    // FUNCTIONS
    // local log-likelihood function
    double gammaLogLike(const MatrixXd &X, const MatrixXd &alpha, const MatrixXd &beta);

    double poisLogLike(const MatrixXi &X, const MatrixXd &lambda);

    double ZIpoisLogLike(const MatrixXi &X, const MatrixXd &lambda, const MatrixXd &pi);

    // Saturated Poisson log-likelihood for deviance
    double poisLoglikeSaturated(const MatrixXi &X, const MatrixXd &lambda);

    // deviance between estimated Poisson and saturated Poisson models
    double poisDeviance(const MatrixXi &X, const MatrixXd &lambda, const MatrixXd &lambda0);

}

#endif
