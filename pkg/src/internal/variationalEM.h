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

#ifndef variationalEM_H
#define variationalEM_H

/*!
* \file variationalEM.h
* \brief template class definition for variationalEM algorithm
* \author Ghislain Durif
* \version 0.1
* \date 06/06/2016
*/

#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include "variational.h"

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision

using std::vector;                      // vector template

/*!
* \namespace countMatrixFactor
*
* A common namespace for the entire package
*/
namespace countMatrixFactor {

    /*!
     * \class variationalEM
     * \brief class to process EM algorithm with variational inference in Gamma Poisson Factor Model
     * @tparam model a gamma Poisson factor model
     *
     * EM algo:
     *  E-step = variational
     *  M-step = updates of prior parameters
     */
    template <typename model>
    class variationalEM : public variational<model> {
    protected:

        int m_iterMax_Estep;        /*!< max nb of iterations in E-step (variational) */
        int m_iterMax_Mstep;        /*!< max nb of iterations in M-step */

        VectorXd m_nbIter_Estep;    /*!< nb of iterations in each E-step */
        VectorXd m_nbIter_Mstep;    /*!< nb of iterations in each M-step */

        int m_globalIter;       /*!< current inner iterations (counting all E-steps and M-steps) */
        int m_nbGlobalIter;     /*!< number of inner effective iterations */

        vector<bool> m_converged_Estep;       /*!< status of convergence of each E-step */
        vector<bool> m_converged_Mstep;       /*!< status of convergence of each M-step */

        VectorXd m_normGap_Estep;         /*!< normalized gap between two iterates (to assess convergence) in E-step */
        VectorXd m_normGap_Mstep;         /*!< normalized gap between two iterates (to assess convergence) in M-step */


    public:
        /*!
        * \brief Constructor of the template class variationalEM for gamma Poisson Factor Model
        * @tparam model a gamma Poisson factor model
        *
        * \param iterMax maximum number of iterations
        * \param iterMax_Estep maximum number of iterations
        * \param order derivative order on normalized gap to assess convergence
        * (0 is the current value, 1 the first order empirical derivative, 2 the second order empirical derivative)
        * \param stabRange range of stability (number of iterations where parameter values are stable to confirm convergence)
        * \param epsilon precision for comparison when assessing convergence
        * \param verbose boolean indicating verbosity in the output
        * \param n number of observations (rows)
        * \param p number of variables (columns)
        * \param K dimension of the latent subspace
        * \param X n x p, count data matrix (const reference)
        * \param phi1 n x K, initial values of first parameter of Gamma distribution on U (const reference)
        * \param phi2 n x K, initial values of second parameter of Gamma distribution on U (const reference)
        * \param theta1 n x K, initial values of first parameter of Gamma distribution on V (const reference)
        * \param theta2 n x K, initial values of second parameter of Gamma distribution on V (const reference)
        * \param alpha1 n x K, initial values of first parameter of Gamma prior on U (const reference)
        * \param alpha2 n x K, initial values of second parameter of Gamma prior on U (const reference)
        * \param beta1 n x K, initial values of first parameter of Gamma prior on V (const reference)
        * \param beta2 n x K, initial values of second parameter of Gamma prior on V (const reference)
        */
        variationalEM(int iterMax, int iterMax_Estep, int iterMax_Mstep, int order,
                      int stabRange, double epsilon, bool verbose,
                      int n, int p, int K, bool ZI, const MatrixXi &X,
                      const MatrixXd &phi1, const MatrixXd &phi2,
                      const MatrixXd &theta1, const MatrixXd &theta2,
                      const MatrixXd &alpha1, const MatrixXd &alpha2,
                      const MatrixXd &beta1, const MatrixXd &beta2);

        /*!
        * \brief Constructor of the template class variationalEM for Penalized (l2) gamma Poisson Factor Model
        * @tparam model a gamma Poisson factor model
        *
        * \param n number of observations (rows)
        * \param p number of variables (columns)
        * \param K dimension of the latent subspace
        * \param iterMax maomegamum number of iterations
        * \param order derivative order on normalized gap to assess convergence
        * (0 is the current value, 1 the first order empirical derivative, 2 the second order empirical derivative)
        * \param stabRange range of stability (number of iterations where parameter values are stable to confirm convergence)
        * \param epsilon precision for comparison when assessing convergence
        * \param verbose boolean indicating verbosity in the output
        * \param X n x p, count data matrix (const reference)
        * \param phi1 n x K, initial values of first parameter of Gamma distribution on U (const reference)
        * \param phi2 n x K, initial values of second parameter of Gamma distribution on U (const reference)
        * \param theta1 n x K, initial values of first parameter of Gamma distribution on V (const reference)
        * \param theta2 n x K, initial values of second parameter of Gamma distribution on V (const reference)
        * \param alpha1 n x K, initial values of first parameter of Gamma prior on U (const reference)
        * \param alpha2 n x K, initial values of second parameter of Gamma prior on U (const reference)
        * \param beta1 n x K, initial values of first parameter of Gamma prior on V (const reference)
        * \param beta2 n x K, initial values of second parameter of Gamma prior on V (const reference)
        * \param r_theta2 vector of l2 penalty constraint on theta2
        * \param r_phi2 vector of l2 penalty constraint on phi2
        */
        variationalEM(int iterMax, int iterMax_Estep, int iterMax_Mstep, int order,
                      int stabRange, double epsilon, bool verbose,
                      int n, int p, int K, bool ZI, const MatrixXi &X,
                      const MatrixXd &phi1, const MatrixXd &phi2,
                      const MatrixXd &theta1, const MatrixXd &theta2,
                      const MatrixXd &alpha1, const MatrixXd &alpha2,
                      const MatrixXd &beta1, const MatrixXd &beta2,
                      const VectorXd &lambda_k, const VectorXd &mu_k);

        /*!
        * \brief Destructor
        *
        * Destructor of the template class
        */
        ~variationalEM();

    public:

        // run algorithm
        void Estep();
        void Mstep();
        void algorithm();

        // create list of object to return
        void returnObject(Rcpp::List &results);

        // assess convergence
        void assessConvergence(int &nstab);
        void assessConvergenceEstep(int iter, int &nstab);
        void assessConvergenceMstep(int iter, int &nstab);

    };

// ############################################################################################## //

    // Class implementation

    // CONSTRUCTOR 1
    template <typename model>
    variationalEM<model>::variationalEM(int iterMax, int iterMax_Estep, int iterMax_Mstep, int order,
                                        int stabRange, double epsilon, bool verbose,
                                        int n, int p, int K, bool ZI, const MatrixXi &X,
                                        const MatrixXd &phi1, const MatrixXd &phi2,
                                        const MatrixXd &theta1, const MatrixXd &theta2,
                                        const MatrixXd &alpha1, const MatrixXd &alpha2,
                                        const MatrixXd &beta1, const MatrixXd &beta2)
    : variational<model>(iterMax, order, stabRange, epsilon, verbose,
                          n, p, K, ZI, X,
                          phi1, phi2, theta1, theta2,
                          alpha1, alpha2, beta1, beta2),
      m_converged_Estep(iterMax, false),
      m_converged_Mstep(iterMax, false)
    {
        m_iterMax_Estep = iterMax_Estep;
        m_iterMax_Mstep = iterMax_Mstep;

        m_nbIter_Estep = VectorXd::Zero(iterMax);
        m_nbIter_Mstep = VectorXd::Zero(iterMax);

        m_globalIter = 0;

        m_normGap_Estep = VectorXd::Zero(iterMax * iterMax_Estep);
        m_normGap_Mstep = VectorXd::Zero(iterMax * iterMax_Mstep);

        // criteria
        this->m_expVar0 = VectorXd::Zero(iterMax * (iterMax_Estep + iterMax_Mstep));
        this->m_expVarU = VectorXd::Zero(iterMax * (iterMax_Estep + iterMax_Mstep));
        this->m_expVarV = VectorXd::Zero(iterMax * (iterMax_Estep + iterMax_Mstep));

        this->m_margLogLike = VectorXd::Zero(iterMax * (iterMax_Estep + iterMax_Mstep));
        this->m_condLogLike = VectorXd::Zero(iterMax * (iterMax_Estep + iterMax_Mstep));
        this->m_priorLogLike = VectorXd::Zero(iterMax * (iterMax_Estep + iterMax_Mstep));
        this->m_postLogLike = VectorXd::Zero(iterMax * (iterMax_Estep + iterMax_Mstep));
        this->m_compLogLike = VectorXd::Zero(iterMax * (iterMax_Estep + iterMax_Mstep));
        this->m_elbo = VectorXd::Zero(iterMax * (iterMax_Estep + iterMax_Mstep));
        this->m_deviance = VectorXd::Zero(iterMax * (iterMax_Estep + iterMax_Mstep));
    }

    // CONSTRUCTOR 2
    template <typename model>
    variationalEM<model>::variationalEM(int iterMax, int iterMax_Estep, int iterMax_Mstep, int order,
                                    int stabRange, double epsilon, bool verbose,
                                    int n, int p, int K, bool ZI, const MatrixXi &X,
                                    const MatrixXd &phi1, const MatrixXd &phi2,
                                    const MatrixXd &theta1, const MatrixXd &theta2,
                                    const MatrixXd &alpha1, const MatrixXd &alpha2,
                                    const MatrixXd &beta1, const MatrixXd &beta2,
                                    const VectorXd &lambda_k, const VectorXd &mu_k)
    : variational<model>(iterMax, order, stabRange, epsilon, verbose,
                          n, p, K, ZI, X,
                          phi1, phi2, theta1, theta2,
                          alpha1, alpha2, beta1, beta2,
                          lambda_k, mu_k),
      m_converged_Estep(iterMax, false),
      m_converged_Mstep(iterMax, false)
    {
        m_iterMax_Estep = iterMax_Estep;
        m_iterMax_Mstep = iterMax_Mstep;

        m_nbIter_Estep = VectorXd::Zero(iterMax);
        m_nbIter_Mstep = VectorXd::Zero(iterMax);

        m_globalIter = 0;

        m_normGap_Estep = VectorXd::Zero(iterMax * iterMax_Estep);
        m_normGap_Mstep = VectorXd::Zero(iterMax * iterMax_Mstep);

        // criteria
        this->m_expVar0 = VectorXd::Zero(iterMax * (iterMax_Estep + iterMax_Mstep));
        this->m_expVarU = VectorXd::Zero(iterMax * (iterMax_Estep + iterMax_Mstep));
        this->m_expVarV = VectorXd::Zero(iterMax * (iterMax_Estep + iterMax_Mstep));

        this->m_margLogLike = VectorXd::Zero(iterMax * (iterMax_Estep + iterMax_Mstep));
        this->m_condLogLike = VectorXd::Zero(iterMax * (iterMax_Estep + iterMax_Mstep));
        this->m_priorLogLike = VectorXd::Zero(iterMax * (iterMax_Estep + iterMax_Mstep));
        this->m_postLogLike = VectorXd::Zero(iterMax * (iterMax_Estep + iterMax_Mstep));
        this->m_compLogLike = VectorXd::Zero(iterMax * (iterMax_Estep + iterMax_Mstep));
        this->m_elbo = VectorXd::Zero(iterMax * (iterMax_Estep + iterMax_Mstep));
        this->m_deviance = VectorXd::Zero(iterMax * (iterMax_Estep + iterMax_Mstep));
    }

    // DESTRUCTOR
    template <typename model>
    variationalEM<model>::~variationalEM() {}


    /*!
     * \brief compute E-step (variational) in EM algo
     * @tparam model a gamma Poisson factor model
     */
    template <typename model>
    void variationalEM<model>::Estep() {

        // Iteration
        int nstab = 0; // number of successive iteration where the normalized gap betwwen two iteration is close to zero (convergence when nstab > rstab)
        int iter = 0;

        while( (iter < m_iterMax_Estep) && (m_converged_Estep[this->m_iter]==false)) {

            if(this->m_verbose==true) {
                // Rcpp::Rcout << "E-step : iter " << iter << std::endl;
            }

            // Multinomial parameters
            // Rcpp::Rcout << "algorithm: Multinomial parameters" << std::endl;
            this->m_model.multinomParam();

            // local parameters
            // U : param phi
            // Rcpp::Rcout << "algorithm: local parameters" << std::endl;
            this->m_model.localParam();

            // global parameters
            // V : param theta
            // Rcpp::Rcout << "algorithm: global parameters" << std::endl;
            this->m_model.globalParam();

            // Poisson rate
            // Rcpp::Rcout << "algorithm: Poisson rate" << std::endl;
            this->m_model.poissonRate();

            // log-likelihood
            // Rcpp::Rcout << "algorithm: loglikelihood" << std::endl;
            this->computeLogLike(m_globalIter);
            // ELBO
            // Rcpp::Rcout << "algorithm: ELBO" << std::endl;
            this->computeELBO(m_globalIter);
            // deviance
            // Rcpp::Rcout << "algorithm: deviance" << std::endl;
            this->computeDeviance(m_globalIter);
            // explained variance
            // Rcpp::Rcout << "algorithm: explained variance" << std::endl;
            this->computeExpVar(m_globalIter);
            // convergence
            // Rcpp::Rcout << "algorithm: convergence ?" << std::endl;
            this->assessConvergenceEstep(iter, nstab);
            // update values of parameters
            // Rcpp::Rcout << "algorithm: next iteration" << std::endl;
            if(m_converged_Estep[this->m_iter] == false && iter < m_iterMax_Estep - 1) {
                this->m_model.nextIterateEstep();
            }
            // increment iteration
            iter++;
            m_globalIter++;
        }
    }

    /*!
    * \brief compute M-step (prior parameter updates) in EM algo
    * @tparam model a gamma Poisson factor model
    */
    template <typename model>
    void variationalEM<model>::Mstep() {

        // Iteration
        int nstab = 0; // number of successive iteration where the normalized gap betwwen two iteration is close to zero (convergence when nstab > rstab)
        int iter = 0;

        while( (iter < m_iterMax_Mstep) && (m_converged_Mstep[this->m_iter]==false)) {

            if(this->m_verbose==true) {
                // Rcpp::Rcout << "M-step : iter " << iter << std::endl;
            }

            // local parameters
            // U : param phi
            // Rcpp::Rcout << "algorithm: local parameters" << std::endl;
            this->m_model.localPriorParam();

            // global parameters
            // V : param theta
            // Rcpp::Rcout << "algorithm: global parameters" << std::endl;
            this->m_model.globalPriorParam();

            /// log-likelihood
            // Rcpp::Rcout << "algorithm: loglikelihood" << std::endl;
            this->computeLogLike(m_globalIter);
            // ELBO
            // Rcpp::Rcout << "algorithm: ELBO" << std::endl;
            this->computeELBO(m_globalIter);
            // deviance
            // Rcpp::Rcout << "algorithm: deviance" << std::endl;
            this->computeDeviance(m_globalIter);
            // explained variance
            // Rcpp::Rcout << "algorithm: explained variance" << std::endl;
            this->computeExpVar(m_globalIter);
            // convergence
            // Rcpp::Rcout << "algorithm: convergence ?" << std::endl;
            this->assessConvergenceMstep(iter, nstab);
            // update values of parameters
            // Rcpp::Rcout << "algorithm: next iteration" << std::endl;
            if(m_converged_Mstep[this->m_iter] == false && iter < m_iterMax_Mstep - 1) {
                this->m_model.nextIterateMstep();
            }
            // increment iteration
            iter++;
            m_globalIter++;
        }
    }



    /*!
     * \brief run EM algorithm
     * @tparam model a gamma Poisson factor model
     */
    template <typename model>
    void variationalEM<model>::algorithm() {

        // Iteration
        int nstab = 0; // number of successive iteration where the normalized gap betwwen two iteration is close to zero (convergence when nstab > rstab)

        while( (this->m_iter < this->m_iterMax) && (this->m_converged==false)) {

            if(this->m_verbose==true) {
                Rcpp::Rcout << "################" << std::endl;
                Rcpp::Rcout << "iter " << this->m_iter << std::endl;
                Rcpp::Rcout << "globalIter " << m_globalIter << std::endl;
            }

            // E-step
            Rcpp::Rcout << "E-step" << std::endl;
            this->Estep();

            // M-step
            Rcpp::Rcout << "M-step" << std::endl;
            this->Mstep();

            // convergence
            Rcpp::Rcout << "algorithm: convergence ?" << std::endl;
            this->assessConvergence(nstab);

            // increment iteration
            this->m_iter++;
        }

    }


    /*!
     * \brief assess convergence in E-step (variational)
     * @tparam model a gamma Poisson factor model
     *
     * @param[in] iter current iteration
     * @param[in,out] nstab number of successive iteration respecting the breaking condition
     */
    template <typename model>
    void variationalEM<model>::assessConvergenceEstep(int iter, int &nstab) {
        // breaking condition: convergence or not
        // Rcpp::Rcout << "norm gap" << std::endl;
        double res = this->m_model.normGapEstep();
        // Rcpp::Rcout << "index = " << this->m_iter * m_iterMax_Estep + iter << std::endl;
        m_normGap_Estep(this->m_iter * m_iterMax_Estep + iter) = res;
        // Rcpp::Rcout << "norm gap = " << m_normGap_Estep(this->m_iter * m_iterMax_Estep + iter) << std::endl;

        // derivative order to consider
        double condition = variational<model>::convCondition(this->m_order, m_normGap_Estep, iter, this->m_iter * m_iterMax_Estep);

        if(std::abs(condition) < this->m_epsilon) {
            nstab++;
        } else {
            nstab=0;
        }

        if(nstab > this->m_stabRange) {
            m_converged_Estep[this->m_iter] = true;
            m_nbIter_Estep(this->m_iter) = iter;
        } else {
            if(iter == m_iterMax_Estep - 1) {
                // Rcpp::Rcout << "convergence ?" << std::endl;
                // Rcpp::Rcout << "iter = " << iter << std::endl;
                // Rcpp::Rcout << "m_iter = " << this->m_iter << std::endl;
                m_nbIter_Estep(this->m_iter) = iter;
            }
        }
    }

    /*!
     * \brief assess convergence in M-step
     * @tparam model a gamma Poisson factor model
     *
     * @param[in] iter current iteration
     * @param[in,out] nstab number of successive iteration respecting the breaking condition
     */
    template <typename model>
    void variationalEM<model>::assessConvergenceMstep(int iter, int &nstab) {
        // breaking condition: convergence or not
        // Rcpp::Rcout << "norm gap" << std::endl;
        double res = this->m_model.normGapMstep();
        // Rcpp::Rcout << "index = " << this->m_iter * m_iterMax_Mstep + iter << std::endl;
        m_normGap_Mstep(this->m_iter * m_iterMax_Mstep + iter) = res;
        // Rcpp::Rcout << "norm gap = " << m_normGap_Mstep(this->m_iter * m_iterMax_Mstep + iter) << std::endl;

        // derivative order to consider
        double condition = variational<model>::convCondition(this->m_order, m_normGap_Mstep, iter, this->m_iter * m_iterMax_Estep);

        if(std::abs(condition) < this->m_epsilon) {
            nstab++;
        } else {
            nstab=0;
        }

        if(nstab > this->m_stabRange) {
            m_converged_Mstep[this->m_iter] = true;
            m_nbIter_Mstep(this->m_iter) = iter;
        } else {
            if(iter == m_iterMax_Mstep - 1) {
                // Rcpp::Rcout << "convergence ?" << std::endl;
                // Rcpp::Rcout << "iter = " << iter << std::endl;
                // Rcpp::Rcout << "m_iter = " << this->m_iter << std::endl;
                m_nbIter_Mstep(this->m_iter) = iter;
            }
        }
    }

    /*!
     * \brief assess convergence of the EM algo
     * @tparam model a gamma Poisson factor model
     *
     * @param[in,out] nstab number of successive iteration respecting the breaking condition
     */
    template <typename model>
    void variationalEM<model>::assessConvergence(int &nstab) {
        // breaking condition: convergence or not
        double res = this->m_model.normGapEM();
        Rcpp::Rcout << "norm gap = " << res << std::endl;
        this->m_normGap(this->m_iter) = res;

        // derivative order to consider
        double condition = variational<model>::convCondition(this->m_order, this->m_normGap, this->m_iter, 0);

        if(std::abs(condition) < this->m_epsilon) {
            nstab++;
        } else {
            nstab=0;
        }

        if(nstab > this->m_stabRange) {
            this->m_converged=true;
            this->m_nbIter=this->m_iter;
            m_nbGlobalIter = m_globalIter;
        } else {
            if(this->m_iter == this->m_iterMax - 1) {
                this->m_nbIter = this->m_iter;
                m_nbGlobalIter = m_globalIter;
            }
        }
    }


    /*!
    * \brief create list of object to return
    * @tparam model a gamma Poisson factor model
    */
    template <typename model>
    void variationalEM<model>::returnObject(Rcpp::List &results) {

        Rcpp::List res;
        this->m_model.returnObject(res);

        Rcpp::List logLikelihood = Rcpp::List::create(Rcpp::Named("margLogLike") = this->m_margLogLike.head(m_nbGlobalIter),
                                                      Rcpp::Named("condLogLike") = this->m_condLogLike.head(m_nbGlobalIter),
                                                      Rcpp::Named("priorLogLike") = this->m_priorLogLike.head(m_nbGlobalIter),
                                                      Rcpp::Named("postLogLike") = this->m_postLogLike.head(m_nbGlobalIter),
                                                      Rcpp::Named("compLogLike") = this->m_compLogLike.head(m_nbGlobalIter),
                                                      Rcpp::Named("elbo") = this->m_elbo.head(m_nbGlobalIter));

        Rcpp::List expVariance = Rcpp::List::create(Rcpp::Named("expVar0") = this->m_expVar0.head(m_nbGlobalIter),
                                                    Rcpp::Named("expVarU") = this->m_expVarU.head(m_nbGlobalIter),
                                                    Rcpp::Named("expVarV") = this->m_expVarV.head(m_nbGlobalIter));

        Rcpp::List EM = Rcpp::List::create(Rcpp::Named("normGap_Estep") = m_normGap_Estep.head(m_nbIter_Estep.sum()+this->m_nbIter),
                                           Rcpp::Named("normGap_Mstep") = m_normGap_Mstep.head(m_nbIter_Mstep.sum()+this->m_nbIter),
                                           Rcpp::Named("nbIter_Estep") = m_nbIter_Estep.head(this->m_nbIter),
                                           Rcpp::Named("nbIter_Mstep") = m_nbIter_Estep.head(this->m_nbIter));

        Rcpp::List returnObj = Rcpp::List::create(Rcpp::Named("logLikelihood") = logLikelihood,
                                                  Rcpp::Named("expVariance") = expVariance,
                                                  Rcpp::Named("EM") = EM,
                                                  Rcpp::Named("normGap") = this->m_normGap.head(this->m_nbIter),
                                                  Rcpp::Named("deviance") = this->m_deviance.head(m_nbGlobalIter),
                                                  Rcpp::Named("converged") = this->m_converged,
                                                  Rcpp::Named("nbIter") = this->m_nbIter);

        SEXP tmp1 = Rcpp::Language("c", res, returnObj).eval();

        SEXP tmp2 = Rcpp::Language("c", results, tmp1).eval();

        results = tmp2;

    }

}

#endif