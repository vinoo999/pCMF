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

#ifndef WRAPPER_TEMPLATE_H
#define WRAPPER_TEMPLATE_H

/*!
* \file wrapper_template.cpp
* \brief definition of a wrapper template for varEM algo in Gamma Poisson Factor Model
* \author Ghislain Durif
* \version 0.1
* \date 05/06/2017
*/

#include <Rcpp.h>
#include <RcppEigen.h>
#include <cstdio>
#include <ctime>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include "variationalEM.h"
#include "random.h"

#include <omp.h>
// [[Rcpp::plugins(openmp)]]

// [[Rcpp::depends(BH)]]
using boost::random::mt19937;
using boost::random::uniform_int_distribution;
using boost::random::variate_generator;

// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision

/*!
 * \namespace countMatrixFactor
 *
 * A common namespace for the entire package
 */
namespace countMatrixFactor {

    /*!
     * \function wrapper_template
     * \brief caller for different inference algo on different GaP factor model
     *
     * @tparam algo an inference algo (variational or variational EM)
     * @tparam model a gamma Poisson factor model
     *
     * \param X n x p, count data matrix (const reference)
     * \param K dimension of the latent subspace
     * \param ZI boolean, indicating if the model is zero-inflated
     * \param phi1 n x K, initial values of first parameter of Gamma distribution on U (const reference)
     * \param phi2 n x K, initial values of second parameter of Gamma distribution on U (const reference)
     * \param theta1 n x K, initial values of first parameter of Gamma distribution on V (const reference)
     * \param theta2 n x K, initial values of second parameter of Gamma distribution on V (const reference)
     * \param alpha1 n x K, initial values of first parameter of Gamma prior on U (const reference)
     * \param alpha2 n x K, initial values of second parameter of Gamma prior on U (const reference)
     * \param beta1 n x K, initial values of first parameter of Gamma prior on V (const reference)
     * \param beta2 n x K, initial values of second parameter of Gamma prior on V (const reference)
     * \param iterMax maximum number of iterations
     * \param iterMin minimum number of iterations
     * \param epsilon precision for comparison when assessing convergence
     * \param order derivative order on normalized gap to assess convergence
     * (0 is the current value, 1 the first order empirical derivative, 2 the second order empirical derivative)
     * \param stabRange range of stability (number of iterations where parameter values are stable to confirm convergence)
     * \param verbose boolean indicating verbosity in the output
     * \param ncores number of cores to use for multi-threading
     * \param nbInit number of initialization to try
     * \param iterMaxInit number of iteration to for each initialization in multi-initialization case (nbInit > 1)
     * \param noise nosie level (between 0 and 1) for multi-init case
     * \param seed for RNG (-1 means no seed)
     *
     */
    template<typename model, template<typename> class algo>
    SEXP wrapper_template(const MatrixXi &X, int K, bool ZI,
                          const MatrixXd &phi01, const MatrixXd &phi02,
                          const MatrixXd &theta01, const MatrixXd &theta02,
                          const MatrixXd &alpha1, const MatrixXd &alpha2,
                          const MatrixXd &beta1, const MatrixXd &beta2,
                          int iterMax, int iterMin, double epsilon,
                          int order, int stabRange, bool verbose, int ncores,
                          int nbInit, int iterMaxInit, double noise, int seed);


    //------------------------------------------------------------------------//
    //                        IMPLEMENTATION                                  //
    //------------------------------------------------------------------------//


    template<typename model, template<typename> class algo>
    SEXP wrapper_template(const MatrixXi &X, int K, bool ZI,
                          const MatrixXd &phi01, const MatrixXd &phi02,
                          const MatrixXd &theta01, const MatrixXd &theta02,
                          const MatrixXd &alpha1, const MatrixXd &alpha2,
                          const MatrixXd &beta1, const MatrixXd &beta2,
                          int iterMax, int iterMin, double epsilon,
                          int order, int stabRange, bool verbose, int ncores,
                          int nbInit, int iterMaxInit, double noise, int seed) {


        int n = X.rows();
        int p = X.cols();

        // parallelizing
#if defined(_OPENMP)
        omp_set_num_threads(ncores);
        Eigen::initParallel();
#endif

        // random seed state
        myRandom::RNGType rng(static_cast<unsigned int>(std::time(0)));
        if(seed >= 0) {
            rng.seed(seed);
        }

        // multi-init case
        if(nbInit > 1) {

            Rcpp::Rcout << "Mutli-initialization" << std::endl;

            // generation of random seed for multiple initialization
            uint32_t* candidateSeeds = new uint32_t[nbInit];
            myRandom::rInt32(candidateSeeds, nbInit, rng);

            // store best SEED id and score
            uint32_t bestSeed = 0;
            int bestSeedId = 0;
            double bestELBO = - std::numeric_limits<double>::max();

            // LOOP OVER MULTI-INITIALIZATION
            for(int initId=0; initId<nbInit; initId++) {

                if(verbose) {
                    Rcpp::Rcout << "Seed " << initId+1 << std::endl;
                }

                // set seed
                rng.seed(candidateSeeds[initId]);

                // generate alpha and beta
                MatrixXd tmpAlpha1(alpha1.rows(), alpha1.cols());
                myRandom::rUnif(tmpAlpha1, alpha1.rows(), alpha1.cols(),
                                (-noise)*alpha1, noise * alpha1, rng);
                tmpAlpha1 = alpha1.array() + tmpAlpha1.array();

                MatrixXd tmpAlpha2(alpha2.rows(), alpha2.cols());
                myRandom::rUnif(tmpAlpha2, alpha2.rows(), alpha2.cols(),
                                (-noise)*alpha2, noise*alpha2, rng);
                tmpAlpha2 = alpha2.array() + tmpAlpha2.array();

                MatrixXd tmpBeta1(beta1.rows(), beta1.cols());
                myRandom::rUnif(tmpBeta1, beta1.rows(), beta1.cols(),
                                (-noise)*beta1, noise * beta1, rng);
                tmpBeta1 = beta1.array() + tmpBeta1.array();

                MatrixXd tmpBeta2(beta2.rows(), beta2.cols());
                myRandom::rUnif(tmpBeta2, beta2.rows(), beta2.cols(),
                                (-noise)*beta2, noise*beta2, rng);
                tmpBeta2 = beta2.array() + tmpBeta2.array();

                // declaration of object model
                // Rcpp::Rcout << "1) Declaration / ";
                algo<model> myModel(iterMaxInit, 1, order,
                                    stabRange, epsilon, 0,
                                    n, p, K, X,
                                    phi01, phi02, theta01, theta02,
                                    tmpAlpha1, tmpAlpha2,
                                    tmpBeta1, tmpBeta2);

                // initialization
                // Rcpp::Rcout << "2) Init / ";
                myModel.Init(rng);

                // computations
                // Rcpp::Rcout << "3) Algo" << std::endl;
                myModel.algorithm();

                // get ELBO
                double currentELBO = myModel.getELBO();

                // comparison
                if(currentELBO > bestELBO) {
                    bestELBO = currentELBO;
                    bestSeed = candidateSeeds[initId];
                    bestSeedId = initId;
                }

            }

            delete[] candidateSeeds;

            // set seed
            rng.seed(bestSeed);

        }


        // declaration of object gamPoisFactorStandard
        Rcpp::Rcout << "Declaration" << std::endl;
        algo<model> myModel(iterMax, iterMin, order,
                            stabRange, epsilon, verbose,
                            n, p, K, X,
                            phi01, phi02, theta01, theta02,
                            alpha1, alpha2, beta1, beta2);

        // initialization
        Rcpp::Rcout << "Initialization" << std::endl;
        myModel.Init(rng);

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

}

#endif // WRAPPER_TEMPLATE_H
