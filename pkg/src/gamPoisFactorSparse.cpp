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

/*!
* \file gamPoisFactorSparse.cpp
* \brief class definition for penalized Gamma Poisson Factor Model
* \author Ghislain Durif
* \version 0.1
* \date 09/05/2016
*/

#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <cstdio>
#include <vector>
#include "gamPoisFactorSparse.h"

#define sqrt() unaryExpr(std::ptr_fun<double,double>(std::sqrt))
#define square() unaryExpr(std::bind2nd(std::pointer_to_binary_function<double,double,double>(std::pow),2))

// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision

using std::vector;

// SRC
namespace countMatrixFactor {

    // CONSTRUCTOR
    gamPoisFactorSparse::gamPoisFactorSparse(int n, int p, int K, int iterMax, int order,
                                       int stabRange, double epsilon, bool verbose,
                                       const MatrixXi &X,
                                       const MatrixXd &phi1, const MatrixXd &phi2,
                                       const MatrixXd &theta1, const MatrixXd &theta2,
                                       const MatrixXd &alpha1, const MatrixXd &alpha2,
                                       const MatrixXd &beta1, const MatrixXd &beta2,
                                       const VectorXd &lambda_k, const VectorXd &mu_k)
        : gamPoisFactorPen::gamPoisFactorPen(n, p, K, iterMax, order,
                                                       stabRange, epsilon, verbose,
                                                       X, phi1, phi2, theta1, theta2,
                                                       alpha1, alpha2, beta1, beta2,
                                                       lambda_k, mu_k)
    {}

    // DESTRUCTOR
    gamPoisFactorSparse::~gamPoisFactorSparse() {}

    //-------------------//
    // parameter updates //
    //-------------------//

    /*!
     * \brief update rule for penalized global parameters theta (factor V) in variational inference
     */
    void gamPoisFactorSparse::sparseGlobalParam() {
        m_theta2cur = m_theta2inter.array() + (m_theta1cur.array().inverse().rowwise() * m_lambda_k.transpose().array());
    }

    //-------------------//
    //     algorithm     //
    //-------------------//

    /*!
    * \brief compute algorithm for variational inference in gamma Poisson factor model
    */
    void gamPoisFactorSparse::algorithm() {

        // Iteration
        int nstab = 0; // number of successive iteration where the normalized gap betwwen two iteration is close to zero (convergence when nstab > rstab)
        int iter = 0;

        while( (iter < m_iterMax) && (m_converged==false)) {

            if(m_verbose==true) {
                Rcpp::Rcout << "iter " << iter << std::endl;
            }

            // Multinomial parameters
            //Rcpp::Rcout << "algorithm: Multinomial parameters" << std::endl;
            this->multinomParam();

            // local parameters
            // U : param phi
            //Rcpp::Rcout << "algorithm: local parameters" << std::endl;
            this->localParam();
            this->penLocalParam();

            // expectation and log-expectation
            Egam(m_phi1cur, m_phi2cur, m_EU);
            Elgam(m_phi1cur, m_phi2cur, m_ElogU);

            // global parameters
            // V : param theta
            //Rcpp::Rcout << "algorithm: global parameters" << std::endl;
            this->globalParam();
            this->sparseGlobalParam();

            // expectation and log-expectation
            Egam(m_theta1cur, m_theta2cur, m_EV);
            Elgam(m_theta1cur, m_theta2cur, m_ElogV);

            // Poisson rate
            //Rcpp::Rcout << "algorithm: Poisson rate" << std::endl;
            this->poissonRate();

            // log-likelihood
            //Rcpp::Rcout << "algorithm: loglikelihood" << std::endl;
            this->computeLogLike(iter);
            // ELBO
            //Rcpp::Rcout << "algorithm: ELBO" << std::endl;
            this->computeELBO(iter);
            // deviance
            //Rcpp::Rcout << "algorithm: deviance" << std::endl;
            this->computeDeviance(iter);
            // explained variance
            //Rcpp::Rcout << "algorithm: explained variance" << std::endl;
            this->computeExpVar(iter);
            // convergence
            //Rcpp::Rcout << "algorithm: convergence ?" << std::endl;
            this-> assessConvergence(iter, nstab);
            // increment values of parameters
            //Rcpp::Rcout << "algorithm: next iteration" << std::endl;
            this->nextIterate();
            // increment iteration
            iter++;
        }
    }

}