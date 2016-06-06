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

#ifndef gamPoisFactor_H
#define gamPoisFactor_H

/*!
 * \file gamPoisFactor.h
 * \brief class definition for standard Gamma Poisson Factor Model
 * \author Ghislain Durif
 * \version 0.2
 * \date 04/05/2016
 */

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::VectorXi;                  // variable size vector, double precision

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

    class gamPoisFactor {

    protected:
        // dimensions
        int m_N;      /*!< number of observations (rows) */
        int m_P;      /*!< number of variables (columns) */
        int m_K;      /*!< dimension of the latent subspace */

        // data
        MatrixXi m_X;           /*!< n x p, count data matrix */
        MatrixXd m_lambda0;     /*!< n x p, Poisson rate for saturated model (X without zeros) */
        MatrixXd m_lambda;      /*!< n x p, Poisson rate */

        // variational parameters
        MatrixXd m_phi1cur;       /*!< n x K, current values of first parameter of Gamma distribution on U */
        MatrixXd m_phi2cur;       /*!< n x K, current values of second parameter of Gamma distribution on U */
        MatrixXd m_theta1cur;     /*!< p x K, current values of first parameter of Gamma distribution on V */
        MatrixXd m_theta2cur;     /*!< p x K, current values of second parameter of Gamma distribution on V */

        MatrixXd m_phi1old;       /*!< n x K, previous values of first parameter of Gamma distribution on U */
        MatrixXd m_phi2old;       /*!< n x K, previous values of second parameter of Gamma distribution on U */
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
        MatrixXd m_alpha1cur;      /*!< n x K, current values of first parameter of Gamma prior on U */
        MatrixXd m_alpha2cur;      /*!< n x K, current values of second parameter of prior Gamma prior on U */
        MatrixXd m_beta1cur;       /*!< p x K, current values of first parameter of prior Gamma prior on V */
        MatrixXd m_beta2cur;       /*!< p x K, current values of second parameter of prior Gamma prior on V */

        MatrixXd m_alpha1old;      /*!< n x K, previous values of first parameter of Gamma prior on U */
        MatrixXd m_alpha2old;      /*!< n x K, previous values of second parameter of prior Gamma prior on U */
        MatrixXd m_beta1old;       /*!< p x K, previous values of first parameter of prior Gamma prior on V */
        MatrixXd m_beta2old;       /*!< p x K, previous values of second parameter of prior Gamma prior on V */

        // order of factors
        VectorXi m_orderDeviance;   /*!< order of factors according to increasing deviance */
        VectorXi m_orderExpVar0;    /*!< order of factors according to increasing expVar0 */
        VectorXi m_orderExpVarU;    /*!< order of factors according to increasing expVarU */
        VectorXi m_orderExpVarV;    /*!< order of factors according to increasing expVarV */

        // criterion following order of factors
        VectorXd m_kDeviance;       /*!< deviance depending on k */
        VectorXd m_kExpVar0;        /*!< expVar0 depending on k */
        VectorXd m_kExpVarU;        /*!< expVarU depending on k */
        VectorXd m_kExpVarV;        /*!< expVarV depending on k */

    public:
        /*!
         * \brief Constructor
         *
         * Constructor of the class gamPoisFactor
         */
        gamPoisFactor(int n, int p, int K,
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

    public:

        // initialization
        void Init();

        // create list with results to be return
        void returnObject(Rcpp::List &results);

        // compute factor order
        void computeOrder();

        //-------------------//
        //      criteria     //
        //-------------------//

        // compute conditonal log-likelihood
        double computeCondLogLike();

        // compute Prior log-likelihood
        double computePriorLogLike();

        // compute Posterior log-likelihood
        double computePostLogLike();

        // compute complete log-likelihood
        double computeCompLogLike();

        // compute complete log-likelihood
        double computeMargLogLike();

        // compute evidence lower bound
        double computeELBO();

        // compute deviance between estimated and saturated model
        double computeDeviance();

        // compute explained variance regarding residuals sum of squares
        double computeExpVar0();

        // compute explained variance regarding U
        double computeExpVarU();

        // compute explained variance regarding V
        double computeExpVarV();

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

        // compute normalized gap between two iterates
        double normGap();

        //-------------------//
        //   order factors   //
        //-------------------//

        // order factors according to expVar0
        void orderExpVar0(VectorXi &order);

        // order factors according to expVarU
        void orderExpVarU(VectorXi &order);

        // order factors according to expVarV
        void orderExpVarV(VectorXi &order);

        // order factors according to deviance
        void orderDeviance(VectorXi &order);

    };

}

#endif