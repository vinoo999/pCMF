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
* \file gamPoisFactorPen.cpp
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
#include "gamPoisFactorPen.h"

#define msqrt() unaryExpr(std::ptr_fun<double,double>(std::sqrt))
#define msquare() unaryExpr(std::bind2nd(std::pointer_to_binary_function<double,double,double>(std::pow),2))

// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision

using std::vector;

// SRC
namespace countMatrixFactor {

    // CONSTRUCTOR
    gamPoisFactorPen::gamPoisFactorPen(int n, int p, int K, int iterMax, int order,
                                       int stabRange, double epsilon, bool verbose,
                                       const MatrixXi &X,
                                       const MatrixXd &phi1, const MatrixXd &phi2,
                                       const MatrixXd &theta1, const MatrixXd &theta2,
                                       const MatrixXd &alpha1, const MatrixXd &alpha2,
                                       const MatrixXd &beta1, const MatrixXd &beta2,
                                       const VectorXd &lambda_k, const VectorXd &mu_k)
        : gamPoisFactorStandard::gamPoisFactorStandard(n, p, K, iterMax, order,
                                                       stabRange, epsilon, verbose,
                                                       X, phi1, phi2, theta1, theta2,
                                                       alpha1, alpha2, beta1, beta2)
        {
            // penalty constant
            m_lambda_k = VectorXd(lambda_k);
            m_mu_k = VectorXd(mu_k);

            // variational parameters
            m_phi2inter = MatrixXd::Zero(n,K);
            m_theta2inter = MatrixXd::Zero(p,K);

        }

    // DESTRUCTOR
    gamPoisFactorPen::~gamPoisFactorPen() {}

    // member functions: documented in src

    /*!
    * \brief Initialization of sufficient statistics
    */
    void gamPoisFactorPen::Init() {

        // Gamma variational parameter
        Rcpp::Rcout << "Init: Gamma variational parameter" << std::endl;
        for(int k=0; k<m_K; k++) {
            // local parameters
            for(int i=0; i<m_N; i++) {
                double param1 = 0;
                double param2 = 0;
                estimParam(1000, m_alpha1(i,k), m_alpha2(i,k), param1, param2);
                m_phi1old(i,k) = param1;
                m_phi1old(i,k) = param1;
                m_phi2inter(i,k) = param2;
            }
            // global parameters
            for(int j=0; j<m_P; j++) {
                double param1 = 0;
                double param2 = 0;
                estimParam(1000, m_beta1(j,k), m_beta2(j,k), param1, param2);
                m_theta1cur(j,k) = param1;
                m_theta1old(j,k) = param1;
                m_theta2inter(j,k) = param2;
            }
        }

        // penalized parameters
        this->penLocalParam();
        this->penGlobalParam();

        m_phi2old = m_phi2cur;
        m_theta2old = m_theta2cur;


        // sufficient statistics
        Rcpp::Rcout << "Init: sufficient statistics" << std::endl;
        Egam(m_phi1cur, m_phi2cur, m_EU);
        Elgam(m_phi1cur, m_phi2cur, m_ElogU);
        Egam(m_theta1cur, m_theta2cur, m_EV);
        Elgam(m_theta1cur, m_theta2cur, m_ElogV);

        // Multinomial Z parameters
        this->multinomParam();

    }

    //-------------------//
    // parameter updates //
    //-------------------//

    /*!
    * \brief update rule for local parameters phi (factor U) in variational inference
    */
    void gamPoisFactorPen::localParam() {
        m_phi1cur = m_alpha1.array() + m_EZ_j.array();
        m_phi2inter = m_alpha2.rowwise() + m_EV.colwise().sum();
    }

    /*!
    * \brief update rule for penalized local parameters phi (factor U) in variational inference
    */
    void gamPoisFactorPen::penLocalParam() {
        m_phi2cur = (m_phi2inter.array() * m_phi1cur.array() + ( (m_phi2inter.array() * m_phi1cur.array()).msquare()
                    + 8 * (m_phi1cur.array().rowwise() * m_mu_k.array())).msqrt()).array() / (2 * m_phi1cur.array());

        // test
        for(int i=0; i<m_N; i++) {
            for(int k = 0; k<m_K; k++) {
                double test = ( m_phi2inter(i,k) * m_phi1cur(i,k) + sqrt( std::pow(m_phi2inter(i,k) * m_phi1cur(i,k),2)
                                + 8 * m_phi1cur(i,k) * m_mu_k(k)) )/(2*m_phi1cur(i,k));
                if(test != m_phi2cur(i,k)) {
                    Rcpp::Rcout << "error in penalized updates" << std::endl;
                }
            }
        }
    }

    /*!
    * \brief update rule for global parameters theta (factor V) in variational inference
    */
    void gamPoisFactorPen::globalParam() {
        m_theta1cur = m_beta1.array() + m_EZ_i.array();
        m_theta2inter = m_beta2.rowwise() + m_EU.colwise().sum();
    }

    /*!
     * \brief update rule for penalized global parameters theta (factor V) in variational inference
     */
    void gamPoisFactorPen::penGlobalParam() {
        m_theta2cur = (m_theta2inter.array() * m_theta1cur.array() + ( (m_theta2inter.array() * m_theta1cur.array()).msquare()
                    + 8 * (m_theta1cur.array().rowwise() * m_mu_k.array())).msqrt()).array() / (2 * m_theta1cur.array());
    }

    //-------------------//
    //     algorithm     //
    //-------------------//

    /*!
    * \brief compute algorithm for variational inference in gamma Poisson factor model
    */
    void gamPoisFactorPen::algorithm() {

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
            this->penGlobalParam();

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

    //       return      //
    //-------------------//

    /*!
    * \brief create list with results to be return
    *
    * @param[out] list containing output
    */
    void gamPoisFactorPen::returnObject(Rcpp::List &results) {
        Rcpp::List logLikelihood = Rcpp::List::create(Rcpp::Named("margLogLike") = m_margLogLike.head(m_nbIter),
                                                      Rcpp::Named("condLogLike") = m_condLogLike.head(m_nbIter),
                                                      Rcpp::Named("priorLogLike") = m_priorLogLike.head(m_nbIter),
                                                      Rcpp::Named("postLogLike") = m_postLogLike.head(m_nbIter),
                                                      Rcpp::Named("compLogLike") = m_compLogLike.head(m_nbIter),
                                                      Rcpp::Named("elbo") = m_elbo.head(m_nbIter));

        Rcpp::List expVariance = Rcpp::List::create(Rcpp::Named("expVar0") = m_expVar0.head(m_nbIter),
                                                    Rcpp::Named("expVarU") = m_expVarU.head(m_nbIter),
                                                    Rcpp::Named("expVarV") = m_expVarV.head(m_nbIter));

        Rcpp::List params = Rcpp::List::create(Rcpp::Named("phi1") = m_phi1cur,
                                               Rcpp::Named("phi2") = m_phi2cur,
                                               Rcpp::Named("theta1") = m_theta1cur,
                                               Rcpp::Named("theta2") = m_theta2cur);

        Rcpp::List stats = Rcpp::List::create(Rcpp::Named("EU") = m_EU,
                                              Rcpp::Named("EV") = m_EV,
                                              Rcpp::Named("ElogU") = m_ElogU,
                                              Rcpp::Named("ElogV") = m_ElogV);

        Rcpp::List order = Rcpp::List::create(Rcpp::Named("orderDeviance") = m_orderDeviance,
                                              Rcpp::Named("orderExpVar0") = m_orderExpVar0,
                                              Rcpp::Named("orderExpVarU") = m_orderExpVarU,
                                              Rcpp::Named("orderExpVarV") = m_orderExpVarV,
                                              Rcpp::Named("kDeviance") = m_kDeviance,
                                              Rcpp::Named("kExpVar0") = m_kExpVar0,
                                              Rcpp::Named("kExpVarU") = m_kExpVarU,
                                              Rcpp::Named("kExpVarV") = m_kExpVarV);

        Rcpp::List returnObj = Rcpp::List::create(Rcpp::Named("U") = m_EU,
                                                  Rcpp::Named("V") = m_EV,
                                                  Rcpp::Named("logLikelihood") = logLikelihood,
                                                  Rcpp::Named("expVariance") = expVariance,
                                                  Rcpp::Named("params") = params,
                                                  Rcpp::Named("stats") = stats,
                                                  Rcpp::Named("order") = order,
                                                  Rcpp::Named("normGap") = m_normGap.head(m_nbIter),
                                                  Rcpp::Named("deviance") = m_deviance.head(m_nbIter),
                                                  Rcpp::Named("converged") = m_converged,
                                                  Rcpp::Named("nbIter") = m_nbIter);

        Rcpp::List pen = Rcpp::List::create(Rcpp::Named("lambda_k") = m_lambda_k,
                                              Rcpp::Named("mu_k") = m_mu_k,
                                              Rcpp::Named("ElogU") = m_ElogU,
                                              Rcpp::Named("ElogV") = m_ElogV);

        SEXP tmp = Rcpp::Language("c", results, returnObj).eval();

        results = tmp;
    }

}