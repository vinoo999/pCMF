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
 * \brief class definition for Gamma Poisson Factor Model
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
     *  Z_{ijk} | X,U,V ~ Multinom(xi_{ijk})
     *  U_{ik} ~ Gamma(phi_{ik})
     *  V_{jk} ~ Gamma(theta_{jk})
     */

    class gamPoisFactor : public loglikelihood, public explainedVariance {
    protected:
        // dimensions
        int m_N;      /*!< number of observations (rows) */
        int m_P;      /*!< number of variables (columns) */
        int m_K;      /*!< dimension of the latent subspace */

        // parameters
        int m_iterMax;          /*!< maximum number of iterations */
        int m_order;            /*!< derivative order on normalized gap to assess convergence
                                (0 is the current value, 1 the first order empirical derivative, 2 the second order empirical derivative) */
        int m_stabRange;        /*!< range of stability (number of iterations where parameter values are stable to confirm convergence) */
        double m_epsilon;       /*!< precision for comparison when assessing convergence */
        bool m_verbose;         /*!< boolean indicating verbosity in the output */

        bool m_converged;       /*!< status of convergence */
        int m_nbIter;           /*!< number of effective iterations */

        // data
        MatrixXi m_X;         /*!< n x p, count data matrix */

        // variational parameters
        gamDistrib m_UphiCur;    /*!< n x K x 2, current values
                                of variational parameters (phi1, phi2)
                                of Gamma distribution on U */

        gamParam m_UphiOld;      /*!< n x K x 2, previous values
                                of variational parameters (phi1, phi2)
                                of Gamma distribution on U */

        gamDistrib m_VthetaCur;  /*!< p x K x 2, current values
                                of variational parameters (theta1, theta2)
                                of Gamma distribution on V */

        gamParam m_VthetaOld;    /*!< p x K, previous values
                                of variational parameters (theta1, theta2)
                                of Gamma distribution on V */

        // sufficient statistics
        MatrixXd m_EZ_i;          /*!< p x k, \sum_i X_{ij} xi_{ijk} = \sum_i E[Z_{ijk}] */
        MatrixXd m_EZ_j;          /*!< n x k, \sum_j X_{ij} xi_{ijk} = \sum_j E[Z_{ijk}] */

        // prior parameter
        gamParam m_alpha;       /*!< n x K, values of first parameter of Gamma prior on U */
        gamParam m_beta;        /*!< p x K, values of first parameter of prior Gamma prior on V */

        // criterion
        VectorXd m_normGap;           /*!< normalized gap between two iterates (to assess convergence) */

    public:
        /*!
         * \brief Constructor
         *
         * Constructor of the class gamPoisFactor
         *
         * \param n number of observations (rows)
         * \param p number of variables (columns)
         * \param K dimension of the latent subspace
         * \param iterMax maximum number of iterations
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
        // member functions: documented in src

        // initialization
        void Init();

        //-------------------//
        //   convergence     //
        //-------------------//

        // difference norm (on parameters)
        double differenceNorm();

        // convergence condition
        double convCondition(int iter, int drift);

        //-------------------//
        // parameter updates //
        //-------------------//

        // Poisson intensity
        void poisRate();

        // local parameters: phi (factor U)
        void localParam();

        // global parameters: theta (factor V)
        void globalParam();


    };





}

#endif