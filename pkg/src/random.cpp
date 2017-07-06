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

/*!
* \file random.cpp
* \brief source file for random number generation functions
* \author Ghislain Durif
* \version 0.1
* \date 05/06/2017
*/

#include <Rcpp.h>
#include <RcppEigen.h>
#include <ctime>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include "random.h"

// [[Rcpp::depends(BH)]]
using boost::random::mt19937;
using boost::random::gamma_distribution;
using boost::random::uniform_int_distribution;
using boost::random::uniform_real_distribution;
using boost::random::variate_generator;

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;


namespace myRandom {

    /*!
     * \function initialize rng
     *
     * \return a Boost rrandom generator
     */
    RNGType rngInit() {
        RNGType rng(static_cast<unsigned int>(std::time(0)));
        return rng;
    }

    /*!
    * \function initialize rng (with a user specified seed)
    *
    * \param[in] seed uint32_t integer (between 0 and 2^32-1)
    *
    * \return a Boost rrandom generator
    */
    RNGType rngInit(uint32_t seed) {
        RNGType rng(seed);
        return rng;
    }


    /*!
     * \function rGamma
     * \brief generate random sample from Gamma distribution
     *
     * Simulate n repetition of Gamma(param1, param2)
     *
     * Note: param1 = shape, param2 = rate
     *
     * \param[out] vec vector to store the simulated values
     * \param[in] n number of values to generate
     * \param[in] param1 shape Gamma parameter
     * \param[in] param1 rate Gamma parameter
     * \param[in] rng Random Number Generator from boost
     */
    void rGamma(VectorXd &vec, int n, double param1, double param2, RNGType &rng) {

        // Gamma generator
        gamma_distribution<> myGamma(param1, 1/param2);
        variate_generator< RNGType &, gamma_distribution<> >generator(rng, myGamma);

        // n is assumed to be the length of vector vec
        for(int ind=0; ind<n; ind++) {
            vec(ind) = generator();
        }
    }


    /*!
     * \function rInt32
     * \brief generate random sample of unsigned int32 integers
     *
     * \param[out] vec vector to store the simulated values
     * \param[in] n number of values to generate
     * \param[in] rng Random Number Generator from boost
     */
    void rInt32(uint32_t* vec, int n, RNGType &rng) {

        // integer generator
        uniform_int_distribution<> candidate(0, std::numeric_limits<uint32_t>::max());
        variate_generator< RNGType &, uniform_int_distribution<> >generator(rng, candidate);

        // n is assumed to be the length of vector vec
        for(int ind=0; ind<n; ind++) {
            vec[ind] = generator();
        }
    }

    /*!
     * \function rUnif
     * \brief generate random sample from uniform distribution
     *
     * Simulate n repetition of Uniform(param1, param2)
     *
     * Note: param1 = min, param2 = max
     *
     * \param[out] vec vector to store the simulated values
     * \param[in] n number of values to generate
     * \param[in] param1 min uniform parameter
     * \param[in] param1 max uniform parameter
     * \param[in] rng Random Number Generator from boost
     */
    void rUnif(VectorXd &vec, int n, double param1, double param2, RNGType &rng) {

        // Uniform generator
        uniform_real_distribution<> myUnif(param1, param2);
        variate_generator< RNGType &, uniform_real_distribution<> >generator(rng, myUnif);

        // n is assumed to be the length of vector vec
        for(int ind=0; ind<n; ind++) {
            vec(ind) = generator();
        }
    }


    /*!
     * \function rUnif
     * \brief generate random sample from uniform distribution
     *
     * Simulate n repetition of Uniform(param1, param2)
     *
     * Note: param1 = min, param2 = max with parameter pair (param1, param2) for each
     * drawning
     *
     * \param[out] mat matrix to store the simulated values
     * \param[in] nrow number of rows in the matrix mat
     * \param[in] ncol number of cols in the matrix mat
     * \param[in] param1 min uniform parameter
     * \param[in] param1 max uniform parameter
     * \param[in] rng Random Number Generator from boost
     */
    void rUnif(MatrixXd &mat, int nrow, int ncol,
               const MatrixXd &param1, const MatrixXd &param2, RNGType &rng) {

        // nrow and ncol are assumed to be consistent with the dimension of the different input matrices
        for(int rowInd=0; rowInd<nrow; rowInd++) {
            for(int colInd=0; colInd<ncol; colInd++) {
                // Uniform generator
                uniform_real_distribution<> myUnif(param1(rowInd, colInd),
                                                   param2(rowInd, colInd));
                variate_generator< RNGType &, uniform_real_distribution<> >generator(rng, myUnif);
                mat(rowInd, colInd) = generator();
            }
        }
    }


}
