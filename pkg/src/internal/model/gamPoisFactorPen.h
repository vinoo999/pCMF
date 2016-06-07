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

#ifndef gamPoisFactorPen_H
#define gamPoisFactorPen_H

/*!
* \file gamPoisFactorPen.h
* \brief class definition for penalized Gamma Poisson Factor Model
* \author Ghislain Durif
* \version 0.1
* \date 09/05/2016
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
    * \class gamPoisFactorPen
    * \brief class to process variational inference in penalized Gamma Poisson Factor Model
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
    * Penalization
    *   l2 penalty on 1/theta_{jk,2}
    *   l2 penalty on 1/phi_{ik,2}
    */

    class gamPoisFactorPen : public gamPoisFactor {

    protected:
        // penalty constant
        VectorXd m_lambda_k;        /*!< penalty constant for l2 constraint on theta_{jk,2} */
        VectorXd m_mu_k;            /*!< penalty constant for l2 constraint on phi_{ik,2} */

        // variational parameters
        MatrixXd m_phi2inter;       /*!< n x K, intermediate values of second parameter of Gamma distribution on U */
        MatrixXd m_theta2inter;     /*!< p x K, intermediate values of second parameter of Gamma distribution on V */

    public:
        /*!
        * \brief Constructor
        *
        * Constructor of the class gamPoisFactorPen
        */
        gamPoisFactorPen(int n, int p, int K, const MatrixXi &X,
                         const MatrixXd &phi1, const MatrixXd &phi2,
                         const MatrixXd &theta1, const MatrixXd &theta2,
                         const MatrixXd &alpha1, const MatrixXd &alpha2,
                         const MatrixXd &beta1, const MatrixXd &beta2,
                         const VectorXd &lambda_k, const VectorXd &mu_k);

        /*!
        * \brief Destructor
        *
        * Destructor of the class gamPoisFactorPen
        */
        ~gamPoisFactorPen();

    public:

        // initialization
        void Init();

        // create list with results to be return
        void returnObject(Rcpp::List &results);

    protected:

        //-------------------//
        // parameter updates //
        //-------------------//

        // local parameters: phi (factor U)
        void localParam();

        // penalized local parameters: phi (factor U)
        void penLocalParam();

        // global parameters: theta (factor V)
        void globalParam();

        // penalized global parameters: theta (factor V)
        void penGlobalParam();

    };

}

#endif