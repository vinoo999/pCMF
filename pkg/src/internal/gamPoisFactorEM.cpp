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

/*!
* \file gamPoisFactorEM.cpp
* \brief class definition for EM algorithm in Gamma Poisson Factor Model
* \author Ghislain Durif
* \version 0.1
* \date 02/06/2016
*/

#include <Rcpp.h>
#include <RcppEigen.h>
#include "gamPoisFactorEM.h"
#include "intermediate.h"

#define mlog() unaryExpr(std::ptr_fun<double,double>(std::log))
#define mpsiInv() unaryExpr(std::bind2nd(std::pointer_to_binary_function<double,int,double>(intermediate::psiInv),6))


// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision

namespace countMatrixFactor {

    // CONSTRUCTOR
    gamPoisFactorEM::gamPoisFactorEM(int n, int p, int K, int iterMax,
                    int iterMax_Estep, int iterMax_Mstep, int order,
                    int stabRange, double epsilon, bool verbose,
                    const MatrixXi &X,
                    const MatrixXd &phi1, const MatrixXd &phi2,
                    const MatrixXd &theta1, const MatrixXd &theta2,
                    const MatrixXd &alpha1, const MatrixXd &alpha2,
                    const MatrixXd &beta1, const MatrixXd &beta2)
        : gamPoisFactorStandard::gamPoisFactorStandard(n, p, K, iterMax, order,
                                                       stabRange, epsilon, verbose,
                                                       X, phi1, phi2, theta1, theta2,
                                                       alpha1, alpha2, beta1, beta2)
    {
        m_iterMax_Estep = iterMax_Estep;
        m_iterMax_Mstep = iterMax_Mstep;

        m_nbIter_Estep = VectorXd::Zero(m_iterMax);
        m_nbIter_Mstep = VectorXd::Zero(m_iterMax);

        m_iter = 0;
        m_globalIter = 0;

        m_converged_Estep = false;
        m_converged_Mstep = false;

        m_normGap_Estep = VectorXd::Zero(iterMax * iterMax_Estep);
        m_normGap_Mstep = VectorXd::Zero(iterMax * iterMax_Mstep);

        // prior parameters
        m_alpha1old = MatrixXd(m_alpha1);
        m_alpha2old = MatrixXd(m_alpha2);

        m_beta1old = MatrixXd(m_beta1);
        m_beta2old = MatrixXd(m_beta2);

        // criteria
        m_expVar0 = VectorXd::Zero(iterMax * iterMax_Estep * iterMax_Mstep);
        m_expVarU = VectorXd::Zero(iterMax * iterMax_Estep * iterMax_Mstep);
        m_expVarV = VectorXd::Zero(iterMax * iterMax_Estep * iterMax_Mstep);

        m_margLogLike = VectorXd::Zero(iterMax * iterMax_Estep * iterMax_Mstep);
        m_condLogLike = VectorXd::Zero(iterMax * iterMax_Estep * iterMax_Mstep);
        m_priorLogLike = VectorXd::Zero(iterMax * iterMax_Estep * iterMax_Mstep);
        m_postLogLike = VectorXd::Zero(iterMax * iterMax_Estep * iterMax_Mstep);
        m_compLogLike = VectorXd::Zero(iterMax * iterMax_Estep * iterMax_Mstep);
        m_elbo = VectorXd::Zero(iterMax * iterMax_Estep * iterMax_Mstep);
        m_deviance = VectorXd::Zero(iterMax * iterMax_Estep * iterMax_Mstep);
    }

    // DESTRUCTOR
    gamPoisFactorEM::~gamPoisFactorEM() {}

    // member functions: documented in src

    /*!
     * \brief compute E-step (variational) in EM algo
     */
    void gamPoisFactorEM::Estep() {
        m_converged_Estep = false;
        m_globalIter--;
        // Iteration
        int nstab = 0; // number of successive iteration where the normalized gap betwwen two iteration is close to zero (convergence when nstab > rstab)
        int iter = 0;

        while( (iter < m_iterMax_Estep) && (m_converged_Estep==false)) {

            if(m_verbose==true) {
                Rcpp::Rcout << "E-step : iter " << iter << std::endl;
            }

            // Multinomial parameters
            //Rcpp::Rcout << "algorithm: Multinomial parameters" << std::endl;
            this->multinomParam();

            // local parameters
            // U : param phi
            //Rcpp::Rcout << "algorithm: local parameters" << std::endl;
            this->localParam();

            // expectation and log-expectation
            Egam(m_phi1cur, m_phi2cur, m_EU);
            Elgam(m_phi1cur, m_phi2cur, m_ElogU);

            // global parameters
            // V : param theta
            //Rcpp::Rcout << "algorithm: global parameters" << std::endl;
            this->globalParam();

            // expectation and log-expectation
            Egam(m_theta1cur, m_theta2cur, m_EV);
            Elgam(m_theta1cur, m_theta2cur, m_ElogV);

            // Poisson rate
            //Rcpp::Rcout << "algorithm: Poisson rate" << std::endl;
            this->poissonRate();

            // log-likelihood
            //Rcpp::Rcout << "algorithm: loglikelihood" << std::endl;
            this->computeLogLike(m_globalIter);
            // ELBO
            //Rcpp::Rcout << "algorithm: ELBO" << std::endl;
            this->computeELBO(m_globalIter);
            // deviance
            //Rcpp::Rcout << "algorithm: deviance" << std::endl;
            this->computeDeviance(m_globalIter);
            // explained variance
            //Rcpp::Rcout << "algorithm: explained variance" << std::endl;
            this->computeExpVar(m_globalIter);
            // convergence
            //Rcpp::Rcout << "algorithm: convergence ?" << std::endl;
            this->assessConvergenceEstep(iter, nstab);
            // increment values of parameters
            //Rcpp::Rcout << "algorithm: next iteration" << std::endl;
            this->nextIterateEstep();
            // increment iteration
            iter++;
            m_globalIter++;
        }
    }

    /*!
     * \brief compute M-step (prior parameter updates) in EM algo
     */
    void gamPoisFactorEM::Mstep() {
        m_converged_Mstep = false;
        m_globalIter--;
        // Iteration
        int nstab = 0; // number of successive iteration where the normalized gap betwwen two iteration is close to zero (convergence when nstab > rstab)
        int iter = 0;

        while( (iter < m_iterMax_Mstep) && (m_converged_Mstep==false)) {

            if(m_verbose==true) {
                Rcpp::Rcout << "M-step : iter " << iter << std::endl;
            }

            // local parameters
            // U : param phi
            //Rcpp::Rcout << "algorithm: local parameters" << std::endl;
            this->localPriorParam();

            // global parameters
            // V : param theta
            //Rcpp::Rcout << "algorithm: global parameters" << std::endl;
            this->globalPriorParam();

            // log-likelihood
            //Rcpp::Rcout << "algorithm: loglikelihood" << std::endl;
            this->computeLogLike(m_globalIter);
            // ELBO
            //Rcpp::Rcout << "algorithm: ELBO" << std::endl;
            this->computeELBO(m_globalIter);
            // deviance
            //Rcpp::Rcout << "algorithm: deviance" << std::endl;
            this->computeDeviance(m_globalIter);
            // explained variance
            //Rcpp::Rcout << "algorithm: explained variance" << std::endl;
            this->computeExpVar(m_globalIter);
            // convergence
            //Rcpp::Rcout << "algorithm: convergence ?" << std::endl;
            this->assessConvergenceEstep(iter, nstab);
            // increment values of parameters
            //Rcpp::Rcout << "algorithm: next iteration" << std::endl;
            this->nextIterateMstep();
            // increment iteration
            iter++;
            m_globalIter++;
        }
    }

    /*!
     * \brief compute EM algorithm in gamma Poisson factor model
     */
    void gamPoisFactorEM::EMalgorithm() {
        m_converged = false;
        // Iteration
        int nstab = 0; // number of successive iteration where the normalized gap betwwen two iteration is close to zero (convergence when nstab > rstab)
        int iter = 0;

        while( (iter < m_iterMax) && (m_converged==false)) {

            if(m_verbose==true) {
                Rcpp::Rcout << "################" << std::endl;
                Rcpp::Rcout << "iter " << iter << std::endl;
            }

            // E-step
            this->Estep();

            // M-step
            this->Mstep();

            // convergence
            //Rcpp::Rcout << "algorithm: convergence ?" << std::endl;
            this->assessConvergence(iter, nstab);

            // increment iteration
            iter++;
        }
    }

    //-------------------//
    // parameter updates //
    //-------------------//

    // local parameters: alpha (factor U)
    void gamPoisFactorEM::localPriorParam() {
        m_alpha1 = (m_alpha2.mlog().rowwise() + m_ElogU.colwise().mean()).mpsiInv();
        m_alpha2 = m_alpha1.array().rowwise() / m_EU.colwise().mean().array();
    }

    // global parameters: beta (factor V)
    void gamPoisFactorEM::globalPriorParam() {
        m_beta1 = (m_beta2.mlog().rowwise() + m_ElogU.colwise().mean()).mpsiInv();
        m_beta2 = m_beta1.array().rowwise() / m_EU.colwise().mean().array();
    }

    // update parameters between iterations in Estep
    void gamPoisFactorEM::nextIterateEstep() {
        m_phi1old = m_phi1cur;
        m_phi2old = m_phi2cur;

        m_theta1old = m_theta1cur;
        m_theta2old = m_theta2cur;
    }

    // update parameters between iterations in Mstep
    void gamPoisFactorEM::nextIterateMstep() {
        m_alpha1old = m_alpha1;
        m_alpha2old = m_alpha2;

        m_beta1old = m_beta1;
        m_beta2old = m_beta2;
    }

    //-------------------//
    //     algorithm     //
    //-------------------//

    /*!
     * \brief assess convergence in E-step (variational)
     */
    void gamPoisFactorEM::assessConvergenceEstep(int iter, int &nstab) {
        // breaking condition: convergence or not
        double paramNorm = sqrt(parameterNorm2(m_phi1old, m_phi2old) + parameterNorm2(m_theta1old, m_theta2old));
        double diffNorm = sqrt(differenceNorm2(m_phi1old, m_phi2old, m_phi1cur, m_phi2cur) + differenceNorm2(m_theta1old, m_theta2old, m_theta1cur, m_theta2cur));

        m_normGap_Estep(m_iter * m_iterMax_Estep + iter) = diffNorm / paramNorm;

        // derivative order to consider
        double condition = convCondition(m_order, m_normGap_Estep, iter, 0);

        if(std::abs(condition) < m_epsilon) {
            nstab++;
        } else {
            nstab=0;
        }

        if(nstab > m_stabRange) {
            m_converged_Estep=true;
            m_nbIter_Estep(m_iter) = iter;
        } else {
            if(iter == m_iterMax_Estep - 1) {
                m_nbIter_Estep(m_iter) = iter;
            }
        }
    }

    /*!
     * \brief assess convergence in M-step
     */
    void gamPoisFactorEM::assessConvergenceMstep(int iter, int &nstab) {
        // breaking condition: convergence or not
        double paramNorm = sqrt(parameterNorm2(m_alpha1old, m_alpha2old) + parameterNorm2(m_beta1old, m_beta2old));
        double diffNorm = sqrt(differenceNorm2(m_alpha1old, m_alpha2old, m_alpha1, m_alpha2) + differenceNorm2(m_beta1old, m_beta2old, m_beta1, m_beta2));

        m_normGap_Mstep(m_iter * m_iterMax_Estep + iter) = diffNorm / paramNorm;

        // derivative order to consider
        double condition = convCondition(m_order, m_normGap_Mstep, iter, 0);

        if(std::abs(condition) < m_epsilon) {
            nstab++;
        } else {
            nstab=0;
        }

        if(nstab > m_stabRange) {
            m_converged_Mstep=true;
            m_nbIter_Mstep(m_iter) = iter;
        } else {
            if(iter == m_iterMax_Mstep - 1) {
                m_nbIter_Mstep(m_iter) = iter;
            }
        }
    }

    /*!
     * \brief assess convergence of the EM algo
     */
    void gamPoisFactorEM::assessConvergence(int iter, int &nstab) {
        // breaking condition: convergence or not
        double paramNorm = sqrt(parameterNorm2(m_phi1old, m_phi2old) + parameterNorm2(m_theta1old, m_theta2old)
                                + parameterNorm2(m_alpha1old, m_alpha2old) + parameterNorm2(m_beta1old, m_beta2old));
        double diffNorm = sqrt(differenceNorm2(m_phi1old, m_phi2old, m_phi1cur, m_phi2cur) + differenceNorm2(m_theta1old, m_theta2old, m_theta1cur, m_theta2cur)
                                + differenceNorm2(m_alpha1old, m_alpha2old, m_alpha1, m_alpha2) + differenceNorm2(m_beta1old, m_beta2old, m_beta1, m_beta2));

        m_normGap(iter) = diffNorm / paramNorm;

        // derivative order to consider
        double condition = convCondition(m_order, m_normGap, iter, 0);

        if(std::abs(condition) < m_epsilon) {
            nstab++;
        } else {
            nstab=0;
        }

        if(nstab > m_stabRange) {
            m_converged=true;
            m_nbIter=iter;
            m_nbGlobalIter = m_globalIter;
        } else {
            if(iter == m_iterMax - 1) {
                m_nbIter = iter;
                m_nbGlobalIter = m_globalIter;
            }
        }
    }


    //-------------------//
    //       return      //
    //-------------------//

    /*!
    * \brief create list with results to be return
    *
    * @param[out] list containing output
    */
    void gamPoisFactorEM::returnObject(Rcpp::List &results) {

        Rcpp::List logLikelihood = Rcpp::List::create(Rcpp::Named("margLogLike") = m_margLogLike.head(m_nbGlobalIter),
                                                      Rcpp::Named("condLogLike") = m_condLogLike.head(m_nbGlobalIter),
                                                      Rcpp::Named("priorLogLike") = m_priorLogLike.head(m_nbGlobalIter),
                                                      Rcpp::Named("postLogLike") = m_postLogLike.head(m_nbGlobalIter),
                                                      Rcpp::Named("compLogLike") = m_compLogLike.head(m_nbGlobalIter),
                                                      Rcpp::Named("elbo") = m_elbo.head(m_nbGlobalIter));

        Rcpp::List expVariance = Rcpp::List::create(Rcpp::Named("expVar0") = m_expVar0.head(m_nbGlobalIter),
                                                    Rcpp::Named("expVarU") = m_expVarU.head(m_nbGlobalIter),
                                                    Rcpp::Named("expVarV") = m_expVarV.head(m_nbGlobalIter));

        Rcpp::List params = Rcpp::List::create(Rcpp::Named("phi1") = m_phi1cur,
                                               Rcpp::Named("phi2") = m_phi2cur,
                                               Rcpp::Named("theta1") = m_theta1cur,
                                               Rcpp::Named("theta2") = m_theta2cur,
                                               Rcpp::Named("alpha1") = m_alpha1,
                                               Rcpp::Named("alpha2") = m_alpha2,
                                               Rcpp::Named("beta1") = m_beta1,
                                               Rcpp::Named("beta2") = m_beta2);

        Rcpp::List stats = Rcpp::List::create(Rcpp::Named("EU") = m_EU,
                                              Rcpp::Named("EV") = m_EV,
                                              Rcpp::Named("ElogU") = m_ElogU,
                                              Rcpp::Named("ElogV") = m_ElogV);

        Rcpp::List order = Rcpp::List::create(Rcpp::Named("orderDeviance") = m_orderDeviance,
                                              Rcpp::Named("orderExpVar0") = m_orderExpVar0,
                                              Rcpp::Named("orderExpVarU") = m_orderExpVarU,
                                              Rcpp::Named("orderExpVarV") = m_orderExpVarV);

        Rcpp::List criteria_k = Rcpp::List::create(Rcpp::Named("kDeviance") = m_kDeviance,
                                                   Rcpp::Named("kExpVar0") = m_kExpVar0,
                                                   Rcpp::Named("kExpVarU") = m_kExpVarU,
                                                   Rcpp::Named("kExpVarV") = m_kExpVarV);

        Rcpp::List EM = Rcpp::List::create(Rcpp::Named("normGap_Estep") = m_normGap_Estep.head(m_nbIter_Estep.sum()),
                                           Rcpp::Named("normGap_Mstep") = m_normGap_Mstep.head(m_nbIter_Mstep.sum()),
                                           Rcpp::Named("nbIter_Estep") = m_normGap_Estep.head(m_nbIter),
                                           Rcpp::Named("nbIter_Mstep") = m_normGap_Estep.head(m_nbIter));

        Rcpp::List returnObj = Rcpp::List::create(Rcpp::Named("U") = m_EU,
                                                  Rcpp::Named("V") = m_EV,
                                                  Rcpp::Named("logLikelihood") = logLikelihood,
                                                  Rcpp::Named("expVariance") = expVariance,
                                                  Rcpp::Named("params") = params,
                                                  Rcpp::Named("stats") = stats,
                                                  Rcpp::Named("EM") = EM,
                                                  Rcpp::Named("order") = order,
                                                  Rcpp::Named("criteria_k") = criteria_k,
                                                  Rcpp::Named("normGap") = m_normGap.head(m_nbIter),
                                                  Rcpp::Named("deviance") = m_deviance.head(m_nbGlobalIter),
                                                  Rcpp::Named("converged") = m_converged,
                                                  Rcpp::Named("nbIter") = m_nbIter);

        SEXP tmp = Rcpp::Language("c", results, returnObj).eval();

        results = tmp;
    }


}
