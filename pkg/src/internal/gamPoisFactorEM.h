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

#ifndef gamPoisFactorEM_H
#define gamPoisFactorEM_H

/*!
* \file gamPoisFactorEM.h
* \brief class definition for EM algorithm in Gamma Poisson Factor Model
* \author Ghislain Durif
* \version 0.1
* \date 02/06/2016
*/

#include <Rcpp.h>
#include <RcppEigen.h>
#include "gamPoisFactorStandard.h"
#include "intermediate.h"

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision

namespace countMatrixFactor {
    /*!
    * \class gamPoisFactorEM
    * \brief class to process EM algorithm with variational inference in Gamma Poisson Factor Model
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
    * EM algo:
    *  E-step = variational
    *  M-step = updates of prior parameters
    */

    class gamPoisFactorEM : public gamPoisFactorStandard {

    protected:
        int m_iterMax_Estep;        /*!< max nb of iterations in E-step (variational) */
        int m_iterMax_Mstep;        /*!< max nb of iterations in M-step */

        VectorXd m_nbIter_Estep;    /*!< nb of iterations in each E-step */
        VectorXd m_nbIter_Mstep;    /*!< nb of iterations in each M-step */

        int m_iter;             /*!< current iterations in EM algo */
        int m_globalIter;       /*!< current inner iterations (counting all E-steps and M-steps) */
        int m_nbGlobalIter;     /*!< number of inner effective iterations */

        bool m_converged_Estep;       /*!< status of convergence in E-step */
        bool m_converged_Mstep;       /*!< status of convergence in M-step */

        VectorXd m_normGap_Estep;         /*!< normalized gap between two iterates (to assess convergence) in E-step */
        VectorXd m_normGap_Mstep;         /*!< normalized gap between two iterates (to assess convergence) in M-step */

        // prior parameters
        MatrixXd m_alpha1old;       /*!< n x K, previous values of first parameter of Gamma prior distribution on U */
        MatrixXd m_alpha2old;       /*!< n x K, previous values of second parameter of Gamma prior distribution on U */

        MatrixXd m_beta1old;     /*!< p x K, previous values of first parameter of Gamma prior distribution on V */
        MatrixXd m_beta2old;     /*!< p x K, previous values of second parameter of Gamma prior distribution on V */

    public:
        /*!
        * \brief Constructor
        *
        * Constructor of the class gamPoisFactorEM
        */
        gamPoisFactorEM(int n, int p, int K, int iterMax,
                        int iterMax_Estep, int iterMax_Mstep, int order,
                        int stabRange, double epsilon, bool verbose,
                        const MatrixXi &X,
                        const MatrixXd &phi1, const MatrixXd &phi2,
                        const MatrixXd &theta1, const MatrixXd &theta2,
                        const MatrixXd &alpha1, const MatrixXd &alpha2,
                        const MatrixXd &beta1, const MatrixXd &beta2);

        /*!
        * \brief Destructor
        *
        * Destructor of the class gamPoisFactorEM
        */
        ~gamPoisFactorEM();

    public:

        // run algorithm
        void Estep();
        void Mstep();
        void EMalgorithm();

        // create list with results to be return
        void returnObject(Rcpp::List &results);

    protected:

        //-------------------//
        // parameter updates //
        //-------------------//

        // local parameters: alpha (factor U)
        void localPriorParam();

        // global parameters: beta (factor V)
        void globalPriorParam();

        // update parameters between iterations
        void nextIterateEstep();

        // update parameters between iterations
        void nextIterateMstep();

        //-------------------//
        //     algorithm     //
        //-------------------//

        // assess convergence
        void assessConvergenceEstep(int iter, int &nstab);
        void assessConvergenceMstep(int iter, int &nstab);
        void assessConvergence(int iter, int &nstab);

    };

}

#endif