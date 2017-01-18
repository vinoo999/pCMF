// Copyright 2016-08 Ghislain Durif
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

#ifndef gamPoisFactorSparse2_H
#define gamPoisFactorSparse2_H

/*!
* \file gamPoisFactorSparse2.h
* \brief class definition for Sparse Gamma Poisson Factor Model
* \author Ghislain Durif
* \version 0.1
* \date 02/08/2016
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
* \class gamPoisFactorSparse2
* \brief class to process variational inference in Sparse Gamma Poisson Factor Model
*
* MODEL
*  X_{ij} = \sum_k Z_{ijk}
*  X | U,V ~ Poisson(U t(V))
*  U ~ Gamma(alpha)
*  V ~ (1-p_{s,k}) delta_0(v_{jk}) + p_{s,k} Gamma(beta)
*  Z_{ij.} | X,U,V ~ Multinom((omega_{ijk})_k)
*  S_{jk} ~ Bernoulli(p_{s,k})
* Variational distribution
*  Z_{ijk} | X,U,V ~ Multinom(omega_{ijk})
*  U_{ik} ~ Gamma(phi_{ik})
*  \tilde{V}_{jk} ~ Gamma(theta_{jk})
*  Y_{ij} ~ Bernoulli(prob_j)
*/

class gamPoisFactorSparse22 : public gamPoisFactor {

protected:

    // Sparse indicator and probabilities
    MatrixXd m_probSparse;           /*!< p x K, matrix of probabilities for variational distribution of S */
    VectorXd m_probSparsePrior;      /*!< vector (size K) of probabilities for prior distribution of S */
    MatrixXd m_S;               /*!< p x K, current states od indicator, depending on the current probabilities */

    // sufficient stats (sum on i)
    MatrixXd m_EZ_logU_i;       /*!< p x K, \sum_i E[Z_{ijk}] * E[log U_{ik}] */
    MatrixXd m_EZ_logV_i;       /*!< p x K, \sum_i E[Z_{ijk}] * E[log V_{jk}] */
    MatrixXd m_EU_EV_i;         /*!< p x K, \sum_i E[U_{ik}] * E[V_{jk}] */
    MatrixXd m_ElgamZ_i;        /*!< p x K, \sum_i E[log(Z_{ijk}!)] */

public:
    /*!
    * \brief Constructor
    *
    * Constructor of the class gamPoisFactorSparse2
    */
    gamPoisFactorSparse2(int n, int p, int K, const MatrixXi &X,
                        const MatrixXd &phi1, const MatrixXd &phi2,
                        const MatrixXd &theta1, const MatrixXd &theta2,
                        const MatrixXd &alpha1, const MatrixXd &alpha2,
                        const MatrixXd &beta1, const MatrixXd &beta2);

    /*!
    * \brief Destructor
    *
    * Destructor of the class gamPoisFactorSparse2
    */
    ~gamPoisFactorSparse2();

public:

    // initialization
    void Init();

    // create list with results to be return
    void returnObject(Rcpp::List &results);

public:
    //-------------------//
    //      criteria     //
    //-------------------//

    // compute evidence lower bound
    double computeELBO();


    //-------------------//
    // parameter updates //
    //-------------------//

    // multinomial parameters
    void multinomParam();

    // local parameters: phi (factor U)
    void localParam();

    // global parameters: theta (factor V)
    void globalParam();

    // sparse proba
    void Sproba();

    // sparse proba in prior
    void priorSproba();

    // parameter update variational standard
    void updateVarational();

    //--------------------------------------//
    // parameter updates for variational EM //
    //--------------------------------------//

    // local parameters: alpha (factor U)
    void localPriorParam();

    // global parameters: beta (factor V)
    void globalPriorParam();

    // parameter update variational EM (E-step)
    void updateEstep();

    // parameter update variational EM (M-step)
    void updateMstep(int iter);

};

}

#endif