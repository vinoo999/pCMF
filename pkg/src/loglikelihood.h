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

#ifndef loglikelihood_H
#define loglikelihood_H

/*!
* \file loglikelihood.h
* \brief class definition for log-likelihood
* \author Ghislain Durif
* \version 0.1
* \date 22/04/2016
*/

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using Eigen::VectorXd;              // variable size matrix, double precision

namespace countMatrixFactor {
    /*!
    * \class gamParam
    * \brief class to define Gamma distribution parameters
    */
    class loglikelihood {
    protected:
        // dimensions
        int m_iterMax;              /*!< maximum number of iterations */

        VectorXd m_margLogLike;     /*!< marginal log-likelihood of the data */
        VectorXd m_condLogLike;     /*!< conditional log-likelihood of the data */
        VectorXd m_priorLogLike;    /*!< log-likelihood of factor priors */
        VectorXd m_postLogLike;     /*!< log-likelihood of factor posterior */
        VectorXd m_compLogLike;     /*!< complete log-likelihood of the model */

        public:
            /*!
            * \brief Constructor
            *
            * Constructor of the class loglikelihood
            */
            loglikelihood(int iterMax);

            /*!
            * \brief Destructor
            *
            * Destructor of the class gamDistrib
            */
            ~loglikelihood();

        public:
            // getter
            void getMarginal(VectorXd &res);
            void getConditional(VectorXd &res);
            void getPrior(VectorXd &res);
            void getPosterior(VectorXd &res);
            void getComplete(VectorXd &res);

            // member functions: doc in src
            void computeLogLike();
        };
}

#endif
