// Copyright 2017-06 Ghislain Durif
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

#ifndef myRANDOM_H
#define myRANDOM_H

/*!
* \file random.h
* \brief header file for random number generation functions
* \author Ghislain Durif
* \version 0.1
* \date 05/06/2017
*/

#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/random/mersenne_twister.hpp>

// [[Rcpp::depends(BH)]]
using boost::random::mt19937;

// [[Rcpp::depends(RcppEigen)]]
using Eigen::VectorXd;


/*!
* \namespace myRandom
*
* A specific namespace for random related functions
*/
namespace myRandom {

    /*!
     * \typedef RNG type (Mersenne-Twister)
     */
    typedef mt19937 RNGType;

    // initialize rng
    RNGType rngInit();

    // initialize rng (with a user specified seed)
    RNGType rngInit(uint32_t seed);


    // generate random sample from Gamma distribution
    void rGamma(VectorXd &vec, int n, double param1, double param2, RNGType rng);


    // generate random sample of unsigned int32 integers
    void rInt32(uint32_t* vec, int n, RNGType rng);

}

#endif // myRANDOM_H
