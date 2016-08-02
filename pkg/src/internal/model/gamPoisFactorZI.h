// Copyright 2016-05 Ghislain Durif
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

#ifndef gamPoisFactorZI_H
#define gamPoisFactorZI_H

/*!
* \file gamPoisFactorZI.h
* \brief class definition for Zero-Inflated Gamma Poisson Factor Model
* \author Ghislain Durif
* \version 0.1
* \date 08/06/2016
*/

#include <Rcpp.h>
#include <RcppEigen.h>
#include "gamPoisFactor.h"

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision

namespace countMatrixFactor {
    /*!
     * \class gamPoisFactorZI
     * \brief class to process variational inference in Zero-Inflated Gamma Poisson Factor Model
     *
     * MODEL
     *  X_{ij} = \sum_k Z_{ijk}
     *  X | U,V ~ Poisson(U t(V))
     *  U ~ Gamma(alpha)
     *  V ~ Gamma(beta)
     *  Z_{ij.} | X,U,V ~ Multinom((omega_{ijk})_k)
     *  D_{ij} ~ Bernoulli(prob0_j)
     * Variational distribution
     *  Z_{ijk} | X,U,V ~ Multinom(omega_{ijk})
     *  U_{ik} ~ Gamma(phi_{ik})
     *  V_{jk} ~ Gamma(theta_{jk})
     *  Y_{ij} ~ Bernoulli(prob_j)
     */

    class gamPoisFactorZI : public gamPoisFactor {

    protected:

        // ZI probabilities and frequencies
        VectorXd m_prob;        /*!< vector of probability for variational distribution of D */
        VectorXd m_prob0;       /*!< vector of probability for prior distribution of D */
        VectorXd m_freq;        /*!< vector of frequence of non null values in each column of D */

        // sufficient stats
        MatrixXd m_EZ_logU_k;       /*!< n x p, \sum_k E[Z_{ijk}] * E[log U_{ik}] */
        MatrixXd m_EZ_logV_k;       /*!< n x p, \sum_k E[Z_{ijk}] * E[log V_{jk}] */
        MatrixXd m_EU_EV_k;         /*!< n x p, \sum_k E[U_{ik}] * E[V_{jk}] */
        MatrixXd m_ElogU_ElogV_k;   /*!< n x p, \sum_k exp(E[log(U_{ik})]) * exp(E[log(V_{jk})]) */
        MatrixXd m_ElgamZ_k;        /*!< n x p, \sum_k E[log(Z_{ijk}!)] */

    public:
        /*!
        * \brief Constructor
        *
        * Constructor of the class gamPoisFactorZI
        */
        gamPoisFactorZI(int n, int p, int K, const MatrixXi &X,
                        const MatrixXd &phi1, const MatrixXd &phi2,
                        const MatrixXd &theta1, const MatrixXd &theta2,
                        const MatrixXd &alpha1, const MatrixXd &alpha2,
                        const MatrixXd &beta1, const MatrixXd &beta2);

        /*!
        * \brief Destructor
        *
        * Destructor of the class gamPoisFactorZI
        */
        ~gamPoisFactorZI();

    public:

        // initialization
        void Init();

        // create list with results to be return
        void returnObject(Rcpp::List &results);

    public:

        //-------------------//
        // parameter updates //
        //-------------------//

        // multinomial parameters
        void multinomParam();

        // local parameters: phi (factor U)
        void localParam();

        // global parameters: theta (factor V)
        void globalParam();

        // zi proba
        void ZIproba();

        // zi proba in prior
        void priorZIproba();

        // parameter update variational standard
        void updateVarational();

        // parameter update variational EM (E-step)
        void updateEstep();

        // parameter update variational EM (M-step)
        void updateMstep();

    };

}

#endif