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

#ifndef gamDistrib_H
#define gamDistrib_H

/*!
 * \file gamDistrib.h
 * \brief function definition for Gamma distribution
 * \author Ghislain Durif
 * \version 0.1
 * \date 22/04/2016
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include "random.h"

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;                  // variable size matrix, double precision

namespace countMatrixFactor {

    // expectation
    void Egam(const MatrixXd &param1, const MatrixXd &param2, MatrixXd &res);

    // log-expectation
    void Elgam(const MatrixXd &param1, const MatrixXd &param2, MatrixXd &res);

    // entropy
    void entropyGam(const MatrixXd &param1, const MatrixXd &param2, MatrixXd &res);

    // estimate shape and rate parameters
    void estimParam(double n, double param1, double param2, double &param1e, double &param2e, myRandom::RNGType &rng);

}

#endif
