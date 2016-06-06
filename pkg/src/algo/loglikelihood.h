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
using Eigen::MatrixXd;              // variable size matrix, double precision
using Eigen::MatrixXi;              // variable size matrix, integer
using Eigen::VectorXd;              // variable size matrix, double precision

namespace countMatrixFactor {
    /*!
    * \class loglikelihood
    * \brief class to define all attributes and member functions regarding log-likelihood
    */
    class loglikelihood {
    protected:

        VectorXd m_margLogLike;     /*!< marginal log-likelihood of the data */
        VectorXd m_condLogLike;     /*!< conditional log-likelihood of the data */
        VectorXd m_priorLogLike;    /*!< log-likelihood of factor priors */
        VectorXd m_postLogLike;     /*!< log-likelihood of factor posterior */
        VectorXd m_compLogLike;     /*!< complete log-likelihood of the model */
        VectorXd m_elbo;            /*!< Evidence lower bound of the model */
        VectorXd m_deviance;        /*!< deviance between estimated and saturated model */

        public:
            /*!
            * \brief Constructor
            *
            * Constructor of the class loglikelihood
            */
            loglikelihood(int size);

            /*!
            * \brief Destructor
            *
            * Destructor of the class gamDistrib
            */
            ~loglikelihood();

        public:
            // getter
            void getMarginal(VectorXd &res, int size);
            void getConditional(VectorXd &res, int size);
            void getPrior(VectorXd &res, int size);
            void getPosterior(VectorXd &res, int size);
            void getComplete(VectorXd &res, int size);
            void getELBO(VectorXd &res, int size);
            void getDeviance(VectorXd &res, int size);

            /*!
             * \brief compute all different log-likelihood
             *
             * Pure virtual member function, to be implemented, depending on the model
             *
             * @param[in] iter current iteration
             */
            virtual void computeLogLike(int iter) = 0;

            /*!
             * \brief compute evidence lower bound
             *
             * Pure virtual member function, to be implemented, depending on the model
             *
             * @param[in] iter current iteration
             */
            virtual void computeELBO(int iter) = 0;

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
             *
             * @param[in] iter current iteration
             */
            virtual void computeDeviance(int iter) = 0;

            /*!
             * \brief deviance between estimated and saturated model for Poisson model
             *
             * Pure virtual member function, to be implemented, depending on the model
             */
            virtual double deviance() = 0;
        };

    // FUNCTIONS
    // local log-likelihood function
    double gammaLogLike(const MatrixXd &X, const MatrixXd &alpha, const MatrixXd &beta);

    double poisLogLike(const MatrixXi &X, const MatrixXd &lambda);

    double ZIpoisLogLike(const MatrixXi &X, const MatrixXd &lambda, const MatrixXd &pi);

    // Saturated Poisson log-likelihood for deviance
    double poisLoglikeSaturated(const MatrixXi &X, const MatrixXd &lambda);

    // deviance between estimated Poisson and saturated Poisson models
    double poisDeviance(const MatrixXi &X, const MatrixXd &lambda, const MatrixXd &lambda0);

}

#endif
