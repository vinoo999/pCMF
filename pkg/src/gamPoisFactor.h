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

    class gamPoisFactor {
    protected:
        // dimensions
        int m_N;      /*!< number of observations (rows) */
        int m_P;      /*!< number of variables (columns) */
        int m_K;      /*!< dimension of the latent subspace */

        // parameters
        int m_iterMax;        /*!< maximum number of iterations */
        int m_order;          /*!< derivative order on normalized gap to assess convergence (0 is the current value, 1 the first order empirical derivative, 2 the second order empirical derivative) */
        int m_stabRange;      /*!< range of stability (number of iterations where parameter values are stable to confirm convergence) */
        double m_epsilon;     /*!< precision for comparison when assessing convergence */
        bool m_verbose;       /*!< boolean indicating verbosity in the output */

        // data
        MatrixXi m_X;         /*!< n x p, count data matrix */

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
        MAtrixXd m_ElogV;         /*!< p x K, Expectation of log V */

        MatrixXd m_EZ_i;          /*!< p x k, \sum_i X_{ij} xi_{ijk} = \sum_i E[Z_{ijk}] */
        MatrixXd m_EZ_j;          /*!< n x k, \sum_j X_{ij} xi_{ijk} = \sum_j E[Z_{ijk}] */

        // prior parameter
        MatriXd m_alpha1;         /*!< n x K, values of first parameter of Gamma prior on U */
        MatrixXd m_alpha2;        /*!< n x K, values of second parameter of prior Gamma prior on U */
        MatrixXd m_beta1;         /*!< p x K, values of first parameter of prior Gamma prior on V */
        MatrixXd m_beta2;         /*!< p x K, values of second parameter of prior Gamma prior on V */

        // criterion
        VectorXd m_margLogLike;       /*!< marginal log-likelihood of the data */
        VectorXd m_condLogLike;       /*!< conditional log-likelihood of the data */
        VectorXd m_priorLogLike;      /*!< log-likelihood of factor priors */
        VectorXd m_postLogLike;       /*!< log-likelihood of factor posterior */
        VectorXd m_compLogLike;       /*!< complete log-likelihood of the model */

        VectorXd m_deviance;          /*!< deviance between estimated and saturated model */
        VectorXd m_normGap;           /*!< normalized gap between two iterates (to assess convergence) */
        VectorXd m_expVar0;           /*!< proportion of explained variance as residuals sum of squares */
        VectorXd m_expVarU;           /*!< proportion of variance explained by columns of U (as the ratio between variance of the projection over total variance) */
        VectorXd m_expVarV;           /*!< proportion of variance explained by columns of V (as the ratio between variance of the projection over total variance) */

    public:
        /*!
         * \brief Constructor
         *
         * Constructor of the class gamPoisFactor
         *
         * \param
         */
        gamPoisFactor(int n, in)


    };





}

#endif