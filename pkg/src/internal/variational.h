// Copyright 2016-06 Ghislain Durif
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

#ifndef VARIATIONAL_H
#define VARIATIONAL_H

/*!
* \file variational.h
* \brief tamplete class definition for variational algorithm
* \author Ghislain Durif
* \version 0.1
* \date 06/06/2016
*/

#include <Rcpp.h>
#include <RcppEigen.h>
#include "internal/explainedVariance.h"
#include "internal/loglikelihood.h"

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
    * \class variational
    * \brief tamplate class to process variational inference on a given model
    */

    template <typename model>
    class variational : public loglikelihood, public explainedVariance {
    protected:

        // MODEL
        model m_model;

        // parameters
        int m_iterMax;          /*!< maximum number of iterations */
        int m_iter;             /*!< current iteration */
        int m_order;            /*!< derivative order on normalized gap to assess convergence
         (0 is the current value, 1 the first order empirical derivative, 2 the second order empirical derivative) */
        int m_stabRange;        /*!< range of stability (number of iterations where parameter values are stable to confirm convergence) */
        double m_epsilon;       /*!< precision for comparison when assessing convergence */
        bool m_verbose;         /*!< boolean indicating verbosity in the output */

        bool m_converged;       /*!< status of convergence */
        int m_nbIter;           /*!< number of effective iterations */

        // criterion
        VectorXd m_normGap;         /*!< normalized gap between two iterates (to assess convergence) */

    public:
        /*!
        * \brief Constructor
        *
        * Constructor of the class variational
        *
        * \param n number of observations (rows)
        * \param p number of variables (columns)
        * \param K dimension of the latent subspace
        * \param iterMax maomegamum number of iterations
        * \param order derivative order on normalized gap to assess convergence
        * (0 is the current value, 1 the first order empirical derivative, 2 the second order empirical derivative)
        * \param stabRange range of stability (number of iterations where parameter values are stable to confirm convergence)
        * \param epsilon precision for comparison when assessing convergence
        * \param verbose boolean indicating verbosity in the output
        * \param X n x p, count data matrix (const reference)
        * \param phi1 n x K, initial values of first parameter of Gamma distribution on U (const reference)
        * \param phi2 n x K, initial values of second parameter of Gamma distribution on U (const reference)
        * \param theta1 n x K, initial values of first parameter of Gamma distribution on V (const reference)
        * \param theta2 n x K, initial values of second parameter of Gamma distribution on V (const reference)
        * \param alpha1 n x K, initial values of first parameter of Gamma prior on U (const reference)
        * \param alpha2 n x K, initial values of second parameter of Gamma prior on U (const reference)
        * \param beta1 n x K, initial values of first parameter of Gamma prior on V (const reference)
        * \param beta2 n x K, initial values of second parameter of Gamma prior on V (const reference)
        */
        variational(int iterMax, int order,
                    int stabRange, double epsilon, bool verbose,
                    int n, int p, int K,
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
        ~variational();

    public:

        // Initialization of sufficient statistics
        void Init();

        // run algorithm
        void algorithm();

        // create list of object to return
        void returnObject(Rcpp::List &results);

        // assess convergence
        void assessConvergence(int &nstab);

        // compute factor order
        void computeOrder();

    protected :

        //-------------------//
        //      criteria     //
        //-------------------//

        // compute log-likelihood
        void computeLogLike(int iter);

        // compute evidence lower bound
        void computeELBO(int iter);

        // compute deviance between estimated and saturated model
        void computeDeviance(int iter);

        // compute explained variance
        void computeExpVar(int iter);

    };

    //-------------------//
    //   convergence     //
    //-------------------//

    // convergence condition
    double convCondition(int order, const VectorXd &normGap, int iter, int drift);

}

#endif