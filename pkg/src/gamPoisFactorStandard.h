// Copyright 2016-04 Ghislain Durif
//
// This file is part of the `hClustering' library for R and related languages.
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

#ifndef gamPoisFactorStandard_H
#define gamPoisFactorStandard_H

/*!
 * \file gamPoisFactorStandard.h
 * \brief class definition for standard Gamma Poisson Factor Model
 * \author Ghislain Durif
 * \version 0.1
 * \date 04/05/2016
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include "gamPoisFactor.h"
#include "intermediate.h"

// [[Rcpp::depends(RcppEigen)]]
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
     * \class gamPoisFactor
     * \brief class to process variational inference in Gamma Poisson Factor Model
     *
     * MODEL
     *  X_{ij} = \sum_k Z_{ijk}
     *  X | U,V ~ Poisson(U t(V))
     *  U ~ Gamma(alpha)
     *  V ~ Gamma(beta)
     *  Z_{ij.} | X,U,V ~ Multinom((omega_{ijk})_k)
     * Variational distribution
     *  Z_{ijk} | X,U,V ~ Multinom(omega_{ijk})
     *  U_{ik} ~ Gamma(phi_{ik})
     *  V_{jk} ~ Gamma(theta_{jk})
     */

    class gamPoisFactorStandard : public gamPoisFactor {

    public:
        /*!
         * \brief Constructor
         *
         * Constructor of the class gamPoisFactor
         */
        gamPoisFactorStandard(int n, int p, int K, int iterMax, int order,
                              int stabRange, double epsilon, bool verbose,
                              const MatrixXi &X,
                              const MatrixXd &phi1, const MatrixXd &phi2,
                              const MatrixXd &theta1, const MatrixXd &theta2,
                              const MatrixXd &alpha1, const MatrixXd &alpha2,
                              const MatrixXd &beta1, const MatrixXd &beta2);

        /*!
         * \brief Destructor
         *
         * Destructor of the class gamPoisFactor
         */
        ~gamPoisFactorStandard();

    public:

        // initialization
        void Init();

        // run algorithm
        void algorithm();

        // create list with results to be return
        void returnObject(Rcpp::List &results);

    protected:

        //-------------------//
        //      criteria     //
        //-------------------//

        // compute log-likelihood
        void computeLogLike(int iter);

        // compute evidence lower bound
        void computeELBO(int iter);

        // evidence lower bound for the specific gamma Poisson factor model
        double ELBO();

        // compute deviance between estimated and saturated model
        void computeDeviance(int iter);

        // deviance between estimated and saturated model for Poisson model
        double deviance();

        // compute explained variance
        void computeExpVar(int iter);

        //-------------------//
        // parameter updates //
        //-------------------//

        // poisson rate
        void poissonRate();

        // multinomial parameters
        void multinomParam();

        // local parameters: phi (factor U)
        void localParam();

        // global parameters: theta (factor V)
        void globalParam();

        // update parameters between iterations
        void nextIterate();

        //-------------------//
        //     algorithm     //
        //-------------------//

        // assess convergence
        void assessConvergence(int iter, int &nstab);

    };

}

#endif