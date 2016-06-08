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

#ifndef VARIATIONAL_H
#define VARIATIONAL_H

/*!
* \file variational.h
* \brief template class definition for variational algorithm
* \author Ghislain Durif
* \version 0.1
* \date 06/06/2016
*/

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision

/*!
* \namespace countMatrixFactor
*
* A common namespace for the entire package
*/
namespace countMatrixFactor {

    /*!
    * \class variational
    * \brief template class to process variational inference on a given model
    * @tparam model a gamma Poisson factor model
    */
    template <typename model>
    class variational {
    protected:

        // MODEL
        model m_model;
        bool m_ZI;

        // parameters
        int m_iterMax;          /*!< maximum number of iterations */
        int m_iter;             /*!< current iteration */
        int m_order;            /*!< derivative order on normalized gap to assess convergence
        (0 is the current value, 1 the first order empirical derivative, 2 the second order empirical derivative) */
        int m_stabRange;        /*!< range of stability (number of iterations where parameter values are stable to confirm convergence) */
        double m_epsilon;       /*!< precision for comparison when assessing convergence */
        bool m_verbose;         /*!< boolean indicating verbosity in the output */

        bool m_converged;       /*!< status of convergence */
        int m_nbIter;           /*!< number of effective iterations */

        // criterion
        VectorXd m_normGap;         /*!< normalized gap between two iterates (to assess convergence) */

        VectorXd m_expVar0;         /*!< proportion of explained variance as residual sum of squares */
        VectorXd m_expVarU;         /*!< proportion of variance explained by columns of U (as the ratio between variance of the projection over total variance) */
        VectorXd m_expVarV;         /*!< proportion of variance explained by columns of V (as the ratio between variance of the projection over total variance) */

        VectorXd m_margLogLike;     /*!< marginal log-likelihood of the data */
        VectorXd m_condLogLike;     /*!< conditional log-likelihood of the data */
        VectorXd m_priorLogLike;    /*!< log-likelihood of factor priors */
        VectorXd m_postLogLike;     /*!< log-likelihood of factor posterior */
        VectorXd m_compLogLike;     /*!< complete log-likelihood of the model */
        VectorXd m_elbo;            /*!< Evidence lower bound of the model */
        VectorXd m_deviance;        /*!< deviance between estimated and saturated model */

    public:
        /*!
         * \brief Constructor of the template class variational for gamma Poisson Factor Model
         * @tparam model a gamma Poisson factor model
         *
         * \param iterMax maximum number of iterations
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
        variational(int iterMax, int order,
                    int stabRange, double epsilon, bool verbose,
                    int n, int p, int K, bool ZI,
                    const MatrixXi &X,
                    const MatrixXd &phi1, const MatrixXd &phi2,
                    const MatrixXd &theta1, const MatrixXd &theta2,
                    const MatrixXd &alpha1, const MatrixXd &alpha2,
                    const MatrixXd &beta1, const MatrixXd &beta2);

        /*!
         * \brief Constructor of the template class variational for Penalized (l2) gamma Poisson Factor Model
         * @tparam model a gamma Poisson factor model
         *
         * \param iterMax maximum number of iterations
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
         * \param r_theta2 vector of l2 constant penalty over theta_{2}, one for each factor
         * \param r_phi2 vector of l2 constant penalty over pi_{2}, one for each factor
         * \param r_theta2 vector of l2 penalty constraint on theta2
         * \param r_phi2 vector of l2 penalty constraint on phi2
         */
        variational(int iterMax, int order,
                    int stabRange, double epsilon, bool verbose,
                    int n, int p, int K, bool ZI,
                    const MatrixXi &X,
                    const MatrixXd &phi1, const MatrixXd &phi2,
                    const MatrixXd &theta1, const MatrixXd &theta2,
                    const MatrixXd &alpha1, const MatrixXd &alpha2,
                    const MatrixXd &beta1, const MatrixXd &beta2,
                    const VectorXd &r_theta2, const VectorXd &r_phi2);

        /*!
         * \brief Destructor
         *
         * Destructor of the template class
         */
        ~variational();

    public:

        // Initialization of model
        void Init();

        // run algorithm
        void algorithm();

        // create list of object to return
        void returnObject(Rcpp::List &results);

        // assess convergence
        void assessConvergence(int &nstab);

        // compute factor order
        void computeOrder();

    protected :

        //-------------------//
        //      criteria     //
        //-------------------//

        // compute log-likelihood
        void computeLogLike(int iter);

        // compute evidence lower bound
        void computeELBO(int iter);

        // compute deviance between estimated and saturated model
        void computeDeviance(int iter);

        // compute explained variance
        void computeExpVar(int iter);

        //-------------------//
        //   convergence     //
        //-------------------//

        // convergence condition
        static double convCondition(int order, const VectorXd &normGap, int iter, int topArray);

    };

// ############################################################################################## //

    // Class implementation

    // CONSTRUCTOR 1
    template <typename model>
    variational<model>::variational(int iterMax, int order,
                                    int stabRange, double epsilon, bool verbose,
                                    int n, int p, int K, bool ZI,
                                    const MatrixXi &X,
                                    const MatrixXd &phi1, const MatrixXd &phi2,
                                    const MatrixXd &theta1, const MatrixXd &theta2,
                                    const MatrixXd &alpha1, const MatrixXd &alpha2,
                                    const MatrixXd &beta1, const MatrixXd &beta2)
    : m_model(n, p, K,
              X, phi1, phi2,
              theta1, theta2,
              alpha1, alpha2,
              beta1, beta2)
    {
        m_ZI = ZI;

        // parameters
        m_iterMax = iterMax;
        m_iter = 0;
        m_order = order;
        m_stabRange = stabRange;
        m_epsilon = epsilon;
        m_verbose = verbose;

        m_converged = false;
        m_nbIter = 0;

        // criterion
        m_normGap = VectorXd::Zero(iterMax);

        m_expVar0 = VectorXd::Zero(iterMax);
        m_expVarU = VectorXd::Zero(iterMax);
        m_expVarV = VectorXd::Zero(iterMax);

        m_margLogLike = VectorXd::Zero(iterMax);
        m_condLogLike = VectorXd::Zero(iterMax);
        m_priorLogLike = VectorXd::Zero(iterMax);
        m_postLogLike = VectorXd::Zero(iterMax);
        m_compLogLike = VectorXd::Zero(iterMax);
        m_elbo = VectorXd::Zero(iterMax);
        m_deviance = VectorXd::Zero(iterMax);
    }

    // CONSTRUCTOR 2
    template <typename model>
    variational<model>::variational(int iterMax, int order,
                                    int stabRange, double epsilon, bool verbose,
                                    int n, int p, int K, bool ZI,
                                    const MatrixXi &X,
                                    const MatrixXd &phi1, const MatrixXd &phi2,
                                    const MatrixXd &theta1, const MatrixXd &theta2,
                                    const MatrixXd &alpha1, const MatrixXd &alpha2,
                                    const MatrixXd &beta1, const MatrixXd &beta2,
                                    const VectorXd &r_theta2, const VectorXd &r_phi2)
    : m_model(n, p, K,
              X, phi1, phi2,
              theta1, theta2,
              alpha1, alpha2,
              beta1, beta2,
              r_theta2, r_phi2)
    {
        m_ZI=ZI;

        // parameters
        m_iterMax = iterMax;
        m_iter = 0;
        m_order = order;
        m_stabRange = stabRange;
        m_epsilon = epsilon;
        m_verbose = verbose;

        m_converged = false;
        m_nbIter = 0;

        // criterion
        m_normGap = VectorXd::Zero(iterMax);

        m_expVar0 = VectorXd::Zero(iterMax);
        m_expVarU = VectorXd::Zero(iterMax);
        m_expVarV = VectorXd::Zero(iterMax);

        m_margLogLike = VectorXd::Zero(iterMax);
        m_condLogLike = VectorXd::Zero(iterMax);
        m_priorLogLike = VectorXd::Zero(iterMax);
        m_postLogLike = VectorXd::Zero(iterMax);
        m_compLogLike = VectorXd::Zero(iterMax);
        m_elbo = VectorXd::Zero(iterMax);
        m_deviance = VectorXd::Zero(iterMax);
    }

    // DESTRUCTOR
    template <typename model>
    variational<model>::~variational() {}

    /*!
     * \brief Initialization of sufficient statistics
     * @tparam model a gamma Poisson factor model
     */
    template <typename model>
    void variational<model>::Init() {
        m_model.Init();
    }

    /*!
     * \brief run algorithm
     * @tparam model a gamma Poisson factor model
     */
    template <typename model>
    void variational<model>::algorithm() {

        // Iteration
        int nstab = 0; // number of successive iteration where the normalized gap betwwen two iteration is close to zero (convergence when nstab > rstab)

        while( (m_iter < m_iterMax) && (m_converged==false)) {

            if(m_verbose==true) {
                Rcpp::Rcout << "iter " << m_iter << std::endl;
            }

            // Multinomial parameters
            //Rcpp::Rcout << "algorithm: Multinomial parameters" << std::endl;
            m_model.multinomParam();

            // local parameters
            // U : param phi
            //Rcpp::Rcout << "algorithm: local parameters" << std::endl;
            m_model.localParam();

            // global parameters
            // V : param theta
            //Rcpp::Rcout << "algorithm: global parameters" << std::endl;
            m_model.globalParam();

            // Poisson rate
            //Rcpp::Rcout << "algorithm: Poisson rate" << std::endl;
            m_model.poissonRate();

            // log-likelihood
            //Rcpp::Rcout << "algorithm: loglikelihood" << std::endl;
            this->computeLogLike(m_iter);
            // ELBO
            //Rcpp::Rcout << "algorithm: ELBO" << std::endl;
            this->computeELBO(m_iter);
            // deviance
            //Rcpp::Rcout << "algorithm: deviance" << std::endl;
            this->computeDeviance(m_iter);
            // explained variance
            //Rcpp::Rcout << "algorithm: explained variance" << std::endl;
            this->computeExpVar(m_iter);
            // convergence
            //Rcpp::Rcout << "algorithm: convergence ?" << std::endl;
            this->assessConvergence(nstab);
            // update values of parameters
            //Rcpp::Rcout << "algorithm: next iteration" << std::endl;
            m_model.nextIterate();
            // increment iteration
            m_iter++;
        }

    }

    /*!
     * \brief assess convergence
     * @tparam model a gamma Poisson factor model
     *
     * @param[in] iter current iteration
     * @param[in,out] nstab number of successive iteration respecting the breaking condition
     */
    template <typename model>
    void variational<model>::assessConvergence(int &nstab) {
        // breaking condition: convergence or not
        double res = m_model.normGap();
        m_normGap(m_iter) = res;

        // derivative order to consider
        double condition = variational<model>::convCondition(m_order, m_normGap, m_iter, 0);

        if(std::abs(condition) < m_epsilon) {
            nstab++;
        } else {
            nstab=0;
        }

        if(nstab > m_stabRange) {
            m_converged=true;
            m_nbIter=m_iter;
        } else {
            if(m_iter == m_iterMax - 1) {
                m_nbIter = m_iter;
            }
        }
    }

    /*!
     * \brief compute factor order
     * @tparam model a gamma Poisson factor model
     */
    template <typename model>
    void variational<model>::computeOrder() {
        m_model.computeOrder();
    }

    /*!
     * \brief create list of object to return
     * @tparam model a gamma Poisson factor model
     */
    template <typename model>
    void variational<model>::returnObject(Rcpp::List &results) {

        Rcpp::List res;
        m_model.returnObject(res);

        Rcpp::List logLikelihood = Rcpp::List::create(Rcpp::Named("margLogLike") = m_margLogLike.head(m_nbIter),
                                                      Rcpp::Named("condLogLike") = m_condLogLike.head(m_nbIter),
                                                      Rcpp::Named("priorLogLike") = m_priorLogLike.head(m_nbIter),
                                                      Rcpp::Named("postLogLike") = m_postLogLike.head(m_nbIter),
                                                      Rcpp::Named("compLogLike") = m_compLogLike.head(m_nbIter),
                                                      Rcpp::Named("elbo") = m_elbo.head(m_nbIter));

        Rcpp::List expVariance = Rcpp::List::create(Rcpp::Named("expVar0") = m_expVar0.head(m_nbIter),
                                                    Rcpp::Named("expVarU") = m_expVarU.head(m_nbIter),
                                                    Rcpp::Named("expVarV") = m_expVarV.head(m_nbIter));

        Rcpp::List returnObj = Rcpp::List::create(Rcpp::Named("logLikelihood") = logLikelihood,
                                                  Rcpp::Named("expVariance") = expVariance,
                                                  Rcpp::Named("normGap") = m_normGap.head(m_nbIter),
                                                  Rcpp::Named("deviance") = m_deviance.head(m_nbIter),
                                                  Rcpp::Named("converged") = m_converged,
                                                  Rcpp::Named("nbIter") = m_nbIter);

        SEXP tmp1 = Rcpp::Language("c", res, returnObj).eval();

        SEXP tmp2 = Rcpp::Language("c", results, tmp1).eval();

        results = tmp2;

    }

    //-------------------//
    //      criteria     //
    //-------------------//

    /*!
     * \brief compute log-likelihood
     * @tparam model a gamma Poisson factor model
     *
     * @param[in] iter current iteration
     */
    template <typename model>
    void variational<model>::computeLogLike(int iter) {
        m_condLogLike(iter) = m_model.computeCondLogLike();
        m_priorLogLike(iter) = m_model.computePriorLogLike();
        m_postLogLike(iter) = m_model.computePostLogLike();
        m_compLogLike(iter) = m_condLogLike(iter) + m_postLogLike(iter);
        m_margLogLike(iter) = m_condLogLike(iter) + m_priorLogLike(iter) - m_postLogLike(iter);
    }

    /*!
     * \brief compute evidence lower bound
     * @tparam model a gamma Poisson factor model
     *
     * @param[in] iter current iteration
     */
    template <typename model>
    void variational<model>::computeELBO(int iter) {
        m_elbo(iter) = m_model.computeELBO();
    }

    /*!
     * \brief compute deviance between estimated and saturated model
     * @tparam model a gamma Poisson factor model
     *
     * @param[in] iter current iteration
     */
    template <typename model>
    void variational<model>::computeDeviance(int iter) {
        m_deviance(iter) = m_model.computeDeviance();
    }

    /*!
     * \brief compute explained variance
     * @tparam model a gamma Poisson factor model
     *
     * @param[in] iter current iteration
     */
    template <typename model>
    void variational<model>::computeExpVar(int iter) {
        m_expVar0(iter) = m_model.computeExpVar0();
        m_expVarU(iter) = m_model.computeExpVarU();
        m_expVarV(iter) = m_model.computeExpVarV();
    }

    /*!
    * \brief assess convergence condition
    *
    * convergence assessed on the normalized gap between two iterates
    * order 0: value of normalized gap
    * order 1: first empirical derivative of normalized gap (speed)
    * order 2: second empirical derivative of normalized gap (acceleration)
    *
    * @return value of the condition
    */
    template <typename model>
    double variational<model>::convCondition(int order, const VectorXd &normGap, int iter, int topArray) {
        double condition = 1;

        iter = iter + topArray;

        switch(order) {
            case 0 : {
                condition = normGap(iter);
            } break;

            case 1 : {
                if(iter>+1) {
                    condition = normGap(iter) - normGap(iter-1);
                }
            } break;

            case 2 : {
                if(iter>+2) {
                    condition = normGap(iter) - 2*normGap(iter-1) + normGap(iter-2);
                }
            } break;
        }
        return(condition);
    }

}

#endif