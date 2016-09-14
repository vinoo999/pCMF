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

#ifndef explainedVariance_H
#define explainedVariance_H

/*!
* \file explainedVariance.h
* \brief functions for explained variance computation
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

    // functions, documented in src

    // explained variance as sum_{i,j} (X_{i,j} - (UtV)_{i,j})^2 / sum_{i,j} X_{i,j}^2
    double expVar0(const MatrixXi &X, const MatrixXd &U, const MatrixXd &V);

    // explained variance as t(X.proj) %*% X.proj / t(X) %*% X, where X.proj is the projection of X onto the new components
    // Xproj = U %*% ginv(t(U) %*% U) %*% t(U) %*% X
    double expVarU(const MatrixXi &X, const MatrixXd &U);

    // explained variance as t(X.proj) %*% X.proj / t(X) %*% X, where X.proj is the projection of X onto the new components
    // Xproj = X %*% V %*% ginv(t(V) %*% V) %*% t(V)
    double expVarV(const MatrixXi &X, const MatrixXd &V);
}

#endif
