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

#ifndef gamPoisFactor_H
#define gamPoisFactor_H

/*!
 * \file gamPoisFactor.h
 * \brief class definition for Gamma Poisson Factor Model (abstract class)
 * \author Ghislain Durif
 * \version 0.1
 * \date 22/04/2016
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include "gamDistrib.h"
#include "gamParam.h"
#include "loglikelihood.h"
#include "explainedVariance.h"

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::VectorXi;                  // variable size vector, integer

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
     *
     *  Abstract class: skeleton for the different implementations
     */

    class gamPoisFactor : public loglikelihood, public explainedVariance {
    protected:
        // dimensions
        int m_N;      /*!< number of observations (rows) */
        int m_P;      /*!< number of variables (columns) */
        int m_K;      /*!< dimension of the latent subspace */

        // parameters
        int m_iterMax;          /*!< maomegamum number of iterations */
        int m_order;            /*!< derivative order on normalized gap to assess convergence
                                (0 is the current value, 1 the first order empirical derivative, 2 the second order empirical derivative) */
        int m_stabRange;        /*!< range of stability (number of iterations where parameter values are stable to confirm convergence) */
        double m_epsilon;       /*!< precision for comparison when assessing convergence */
        bool m_verbose;         /*!< boolean indicating verbosity in the output */

        bool m_converged;       /*!< status of convergence */
        int m_nbIter;           /*!< number of effective iterations */

        // data
        MatrixXi m_X;           /*!< n x p, count data matrix */
        MatrixXd m_lambda0;     /*!< n x p, Poisson rate for saturated model (X without zeros) */
        MatrixXd m_lambda;      /*!< n x p, Poisson rate */

        // variational parameters
        MatrixXd m_phi1cur;       /*!< n x K, current values of first parameter of Gamma distribution on U */
        MatrixXd m_phi2cur;       /*!< n x K, current values of second parameter of Gamma distribution on U */

        MatrixXd m_phi1old;       /*!< n x K, previous values of first parameter of Gamma distribution on U */
        MatrixXd m_phi2old;       /*!< n x K, previous values of second parameter of Gamma distribution on U */

        MatrixXd m_theta1cur;     /*!< p x K, current values of first parameter of Gamma distribution on V */
        MatrixXd m_theta2cur;     /*!< p x K, current values of second parameter of Gamma distribution on V */

        MatrixXd m_theta1old;     /*!< p x K, previous values of first parameter of Gamma distribution on V */
        MatrixXd m_theta2old;     /*!< p x K, previous values of second parameter of Gamma distribution on V */

        // sufficient statistics
        MatrixXd m_EU;            /*!< n x K, Expectation of U */
        MatrixXd m_ElogU;         /*!< n x K, Expectation of log U */

        MatrixXd m_EV;            /*!< p x K, Expectation of V */
        MatrixXd m_ElogV;         /*!< p x K, Expectation of log V */

        MatrixXd m_EZ_i;          /*!< p x k, \sum_i X_{ij} xi_{ijk} = \sum_i E[Z_{ijk}] */
        MatrixXd m_EZ_j;          /*!< n x k, \sum_j X_{ij} xi_{ijk} = \sum_j E[Z_{ijk}] */

        // prior parameter
        MatriXd m_alpha1;         /*!< n x K, values of first parameter of Gamma prior on U */
        MatrixXd m_alpha2;        /*!< n x K, values of second parameter of prior Gamma prior on U */
        MatrixXd m_beta1;         /*!< p x K, values of first parameter of prior Gamma prior on V */
        MatrixXd m_beta2;         /*!< p x K, values of second parameter of prior Gamma prior on V */

        // criterion
        VectorXd m_normGap;         /*!< normalized gap between two iterates (to assess convergence) */

        // order of factors
        VectorXi m_orderDeviance;   /*!< order of factors according to increasing deviance */
        VectorXi m_orderExpVar0;    /*!< order of factors according to increasing expVar0 */
        VectorXi m_orderExpVarU;    /*!< order of factors according to increasing expVarU */
        VectorXi m_orderExpVarV;    /*!< order of factors according to increasing expVarV */

    public:
        /*!
         * \brief Constructor
         *
         * Constructor of the class gamPoisFactor
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
        gamPoisFactor(int n, int p, int K, int iterMax, int order,
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
        ~gamPoisFactor();

    protected:

        /*!
         * \brief Initialization of sufficient statistics
         */
        virtual void Init() = 0;

        //-------------------//
        //      criteria     //
        //-------------------//

        /*!
        * \brief compute all different log-likelihood
        *
        * Pure virtual member function, to be implemented, depending on the model
        */
        virtual void computeLogLike() = 0;

        /*!
         * \brief compute evidence lower bound
         *
         * Pure virtual member function, to be implemented, depending on the model
         */
        virtual void computeELBO() = 0;

        /*!
         * \brief evidence lower bound for a specific model
         *
         * Pure virtual member function, to be implemented, depending on the model
         */
        virtual double ELBO() = 0;

        /*!
         * \brief compute deviance between estimated and saturated model
         *
         * Pure virtual member function, to be implemented, depending on the model
         */
        virtual void computeDeviance() = 0;

        /*!
         * \brief deviance between estimated and saturated model for Poisson model
         *
         * Pure virtual member function, to be implemented, depending on the model
         */
        virtual double deviance() = 0;

        // compute explained variance
        double computeExpVar(int iter);

        //-------------------//
        // parameter updates //
        //-------------------//

        // poisson rate
        virtual void poissonRate() = 0;

        // multinomial parameters
        virtual void multinomParam() = 0;

        // local parameters: phi (factor U)
        virtual void localParam() = 0;

        // global parameters: theta (factor V)
        virtual void globalParam() = 0;

        // update parameters between iterations
        void update();

        //-------------------//
        //     algorithm     //
        //-------------------//

        // compute algorithm
        virtual void algorithm() = 0;

        // assess convergence
        void assessConvergence(int iter, int &nstab);

        //-------------------//
        //       return      //
        //-------------------//

        // create list of object to return
        virtual Rcpp::List returnObject() = 0;

    };

    //-------------------//
    //   convergence     //
    //-------------------//

    // parameter squared euclidean norm
    double parameterNorm2(const MatrixXd &param1, const MatrixXd &param2);

    // difference of squared euclidean norm (on parameters)
    double differenceNorm2(const MatrixXd &param1a, const MatrixXd &param2a, const MatrixXd &param1b, const MatrixXd &param2b);

    // convergence condition
    double convCondition(int order, const VectorXd &normGap, int iter, int drift);

}

#endif