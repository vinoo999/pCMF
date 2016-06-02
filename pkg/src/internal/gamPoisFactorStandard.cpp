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

/*!
 * \file gamPoisFactorStandard.cpp
 * \brief class definition for standard Gamma Poisson Factor Model
 * \author Ghislain Durif
 * \version 0.1
 * \date 04/05/2016
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <cstdio>
#include <vector>
#include "gamPoisFactorStandard.h"

#define mexp() unaryExpr(std::ptr_fun<double,double>(std::exp))
#define mlgamma() unaryExpr(std::ptr_fun<double,double>(lgamma))
#define mlog() unaryExpr(std::ptr_fun<double,double>(std::log))

// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision

using std::vector;

// SRC
namespace countMatrixFactor {

    // CONSTRUCTOR
    gamPoisFactorStandard::gamPoisFactorStandard(int n, int p, int K, int iterMax, int order,
                                                 int stabRange, double epsilon, bool verbose,
                                                 const MatrixXi &X,
                                                 const MatrixXd &phi1, const MatrixXd &phi2,
                                                 const MatrixXd &theta1, const MatrixXd &theta2,
                                                 const MatrixXd &alpha1, const MatrixXd &alpha2,
                                                 const MatrixXd &beta1, const MatrixXd &beta2)
                : gamPoisFactor::gamPoisFactor(n, p, K, iterMax, order,
                                              stabRange, epsilon, verbose,
                                              X, phi1, phi2, theta1, theta2,
                                              alpha1, alpha2, beta1, beta2)
    {}

    // DESTRUCTOR
    gamPoisFactorStandard::~gamPoisFactorStandard() {}

    // member functions: documented in src

    /*!
     * \brief Initialization of sufficient statistics
     */
    void gamPoisFactorStandard::Init() {

        // Gamma variational parameter
        Rcpp::Rcout << "Init: Gamma variational parameter" << std::endl;
        for(int k=0; k<m_K; k++) {
            // local parameters
            for(int i=0; i<m_N; i++) {
                double param1 = 0;
                double param2 = 0;
                estimParam(1000, m_alpha1(i,k), m_alpha2(i,k), param1, param2);
                m_phi1cur(i,k) = param1;
                m_phi1old(i,k) = param1;
                m_phi2cur(i,k) = param2;
                m_phi2old(i,k) = param2;
            }
            // global parameters
            for(int j=0; j<m_P; j++) {
                double param1 = 0;
                double param2 = 0;
                estimParam(1000, m_beta1(j,k), m_beta2(j,k), param1, param2);
                m_theta1cur(j,k) = param1;
                m_theta1old(j,k) = param1;
                m_theta2cur(j,k) = param2;
                m_theta2old(j,k) = param2;
            }
        }

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
    //      criteria     //
    //-------------------//

    /*!
     * \brief compute log-likelihood
     *
     * @param[in] iter current iteration
     */
    void gamPoisFactorStandard::computeLogLike(int iter) {
        m_condLogLike(iter) = poisLogLike(m_X, m_lambda);
        m_priorLogLike(iter) = gammaLogLike(m_EU, m_alpha1, m_alpha2) + gammaLogLike(m_EV, m_beta1, m_beta2);
        m_postLogLike(iter) = gammaLogLike(m_EU, m_phi1cur, m_phi2cur) + gammaLogLike(m_EV, m_theta1cur, m_theta2cur);
        m_compLogLike(iter) = m_condLogLike(iter) + m_postLogLike(iter);
        m_margLogLike(iter) = m_condLogLike(iter) + m_priorLogLike(iter) - m_postLogLike(iter);
    }

    /*!
     * \brief compute evidence lower bound
     *
     * @param[in] iter current iteration
     */
    void gamPoisFactorStandard::computeELBO(int iter) {
        m_elbo(iter) = this->ELBO();
    }

    /*!
     * \brief evidence lower bound for the specific gamma Poisson factor model
     *
     * @return value of the ELBO for the current values of estimates
     */
    double gamPoisFactorStandard::ELBO() {

        intermediate::checkExp(m_ElogU);
        intermediate::checkExp(m_ElogV);

        double res1 = (-1) * ( ( (m_X.cast<double>().array() + 1).mlgamma() ).sum() + ( m_EU * m_EV.transpose() ).sum() );
        //Rcpp::Rcout << "ELBO: res1 = " << res1 << std::endl;
        double res2 = ( m_X.cast<double>().array() * (m_ElogU.mexp() * m_ElogV.mexp().transpose()).mlog().array()).sum();
        //Rcpp::Rcout << "ELBO: res2 = " << res2 << std::endl;

        double res3 = ( (m_alpha1.array() - 1) * m_ElogU.array() + m_alpha1.array() * m_alpha2.mlog().array()
                        - m_alpha2.array() * m_EU.array() - m_alpha1.mlgamma().array() ).sum();
        //Rcpp::Rcout << "ELBO: res3 = " << res3 << std::endl;
        double res4 = (-1) * ( (m_phi1cur.array() - 1) * m_ElogU.array() + m_phi1cur.array() * m_phi2cur.mlog().array()
                                   - m_phi2cur.array() * m_EU.array() - m_phi1cur.mlgamma().array() ).sum();
        //Rcpp::Rcout << "ELBO: res4 = " << res4 << std::endl;

        double res5 = ( (m_beta1.array() - 1) * m_ElogV.array() + m_beta1.array() * m_beta2.mlog().array()
                            - m_beta2.array() * m_EV.array() - m_beta1.mlgamma().array() ).sum();
        //Rcpp::Rcout << "ELBO: res5 = " << res5 << std::endl;
        double res6 = (-1) * ( (m_theta1cur.array() - 1) * m_ElogV.array() + m_theta1cur.array() * m_theta2cur.mlog().array()
                                   - m_theta2cur.array() * m_EV.array() - m_theta1cur.mlgamma().array() ).sum();
        //Rcpp::Rcout << "ELBO: res6 = " << res6 << std::endl;

        double res = res1 + res2 + res3 + res4 + res5 + res6;

        return res;
    }

    /*!
     * \brief compute deviance between estimated and saturated model
     *
     * @param[in] iter current iteration
     */
    void gamPoisFactorStandard::computeDeviance(int iter) {
        m_deviance(iter) = this->deviance();
    }

    /*!
     * \brief deviance between estimated and saturated model for Poisson model
     *
     * @return value of the deviance for the current values of estimates
     */
    double gamPoisFactorStandard::deviance() {
        double res = poisDeviance(m_X, m_lambda, m_lambda0);
        return res;
    }

    /*!
     * \brief compute explained variance
     *
     * @param[in] iter current iteration
     */
    void gamPoisFactorStandard::computeExpVar(int iter) {
        m_expVar0(iter) = expVar0(m_X, m_EU, m_EV);
        m_expVarU(iter) = expVarU(m_X, m_EU);
        m_expVarV(iter) = expVarV(m_X, m_EV);
    }

    //-------------------//
    // parameter updates //
    //-------------------//

    /*!
     * \brief update rule for poisson rates in variational inference
     */
    void gamPoisFactorStandard::poissonRate() {
        m_lambda = m_EU * m_EV.transpose();
    }

    /*!
     * \brief update rule for multinomial parameters in variational inference
     */
    void gamPoisFactorStandard::multinomParam() {

        intermediate::checkExp(m_ElogU);
        intermediate::checkExp(m_ElogV);

        m_EZ_j = m_ElogU.mexp().array() * ( (m_X.cast<double>().array() / (m_ElogU.mexp() * m_ElogV.mexp().transpose()).array() ).matrix() * m_ElogV.mexp() ).array();
        m_EZ_i = m_ElogV.mexp().array() * ( (m_X.cast<double>().array() / (m_ElogU.mexp() * m_ElogV.mexp().transpose()).array() ).matrix().transpose() * m_ElogU.mexp() ).array();
    }

    /*!
     * \brief update rule for local parameters phi (factor U) in variational inference
     */
    void gamPoisFactorStandard::localParam() {
        m_phi1cur = m_alpha1.array() + m_EZ_j.array();
        m_phi2cur = m_alpha2.rowwise() + m_EV.colwise().sum();
    }

    /*!
     * \brief update rule for global parameters theta (factor V) in variational inference
     */
    void gamPoisFactorStandard::globalParam() {
        m_theta1cur = m_beta1.array() + m_EZ_i.array();
        m_theta2cur = m_beta2.rowwise() + m_EU.colwise().sum();
    }

    /*!
     * \brief update parameters between iterations
     */
    void gamPoisFactorStandard::nextIterate() {
        m_phi1old = m_phi1cur;
        m_phi2old = m_phi2cur;

        m_theta1old = m_theta1cur;
        m_theta2old = m_theta2cur;
    }


    //-------------------//
    //     algorithm     //
    //-------------------//

    /*!
     * \brief compute algorithm for variational inference in gamma Poisson factor model
     */
    void gamPoisFactorStandard::algorithm() {

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

    /*!
     * \brief assess convergence
     *
     * @param[in] iter current iteration
     * @param[in,out] nstab number of successive iteration respecting the breaking condition
     */
    void gamPoisFactorStandard::assessConvergence(int iter, int &nstab) {
        // breaking condition: convergence or not
        double paramNorm = sqrt(parameterNorm2(m_phi1old, m_phi2old) + parameterNorm2(m_theta1old, m_theta2old));
        double diffNorm = sqrt(differenceNorm2(m_phi1old, m_phi2old, m_phi1cur, m_phi2cur) + differenceNorm2(m_theta1old, m_theta2old, m_theta1cur, m_theta2cur));

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
        } else {
            if(iter == m_iterMax - 1) {
                m_nbIter = iter;
            }
        }
    }

    //-------------------//
    //   order factors   //
    //-------------------//

    /*!
    * \brief order factors according to expVar0
    *
    * @param[out] vector of factor order
    */
    void gamPoisFactorStandard::orderExpVar0(VectorXi &order) {
        vector<int> leftIndex;
        vector<int> chosenIndex;

        MatrixXd Utmp = MatrixXd::Zero(m_EU.rows(), m_EU.cols());
        MatrixXd Vtmp = MatrixXd::Zero(m_EV.rows(), m_EV.cols());

        // init
        for(int k=0; k<m_K; k++) {
            leftIndex.push_back(k);
        }

        // ordering
        for(int k=0; k<m_K; k++) {
            double val_max = -1E6;
            int ind_max = 0;
            int ind_left = 0;
            for(int ind=0; ind<leftIndex.size(); ind++) {
                MatrixXd tmpU(m_EU.rows(), k+1);
                MatrixXd tmpV(m_EV.rows(), k+1);

                tmpU.col(k) = m_EU.col(leftIndex[ind]);
                tmpV.col(k) = m_EV.col(leftIndex[ind]);

                if(k!=0) {
                    tmpU.leftCols(k) = Utmp.leftCols(k);
                    tmpV.leftCols(k) = Vtmp.leftCols(k);
                }

                double res = expVar0(m_X, tmpU, tmpV);
                if(res > val_max) {
                    val_max = res;
                    ind_max = leftIndex[ind];
                    ind_left = ind;
                }
            }
            // max found
            chosenIndex.push_back(ind_max);
            leftIndex.erase(leftIndex.begin() + ind_left);

            Utmp.col(k) = m_EU.col(ind_max);
            Vtmp.col(k) = m_EV.col(ind_max);

            // store it
            m_kExpVar0(k) = val_max;
        }

        // return
        for(int k=0; k<m_K; k++) {
            order(k) = chosenIndex[k] + 1;
        }
    }

    /*!
     * \brief order factors according to expVarU
     *
     * @param[out] vector of factor order
     */
    void gamPoisFactorStandard::orderExpVarU(VectorXi &order) {
        vector<int> leftIndex;
        vector<int> chosenIndex;

        MatrixXd Utmp = MatrixXd::Zero(m_EU.rows(), m_EU.cols());

        // init
        for(int k=0; k<m_K; k++) {
            leftIndex.push_back(k);
        }

        // ordering
        for(int k=0; k<m_K; k++) {
            double val_max = -1E6;
            int ind_max = 0;
            int ind_left = 0;
            for(int ind=0; ind<leftIndex.size(); ind++) {
                MatrixXd tmpU(m_EU.rows(), k+1);
                tmpU.col(k) = m_EU.col(leftIndex[ind]);

                if(k!=0) {
                    tmpU.leftCols(k) = Utmp.leftCols(k);
                }

                double res = expVarU(m_X, tmpU);
                if(res > val_max) {
                    val_max = res;
                    ind_max = leftIndex[ind];
                    ind_left = ind;
                }
            }
            // max found
            chosenIndex.push_back(ind_max);
            leftIndex.erase(leftIndex.begin() + ind_left);

            Utmp.col(k) = m_EU.col(ind_max);

            // store it
            m_kExpVarU(k) = val_max;
        }

        // return
        for(int k=0; k<m_K; k++) {
            order(k) = chosenIndex[k] + 1;
        }
    }

    /*!
     * \brief order factors according to expVarV
     *
     * @param[out] vector of factor order
     */
    void gamPoisFactorStandard::orderExpVarV(VectorXi &order) {
        vector<int> leftIndex;
        vector<int> chosenIndex;

        MatrixXd Vtmp = MatrixXd::Zero(m_EV.rows(), m_EV.cols());

        // init
        for(int k=0; k<m_K; k++) {
            leftIndex.push_back(k);
        }

        // ordering
        for(int k=0; k<m_K; k++) {
            double val_max = -1E6;
            int ind_max = 0;
            int ind_left = 0;
            for(int ind=0; ind<leftIndex.size(); ind++) {
                MatrixXd tmpV(m_EV.rows(), k+1);
                tmpV.col(k) = m_EV.col(leftIndex[ind]);

                if(k!=0) {
                    tmpV.leftCols(k) = Vtmp.leftCols(k);
                }

                double res = expVarV(m_X, tmpV);
                if(res > val_max) {
                    val_max = res;
                    ind_max = leftIndex[ind];
                    ind_left = ind;
                }
            }
            // max found
            chosenIndex.push_back(ind_max);
            leftIndex.erase(leftIndex.begin() + ind_left);

            Vtmp.col(k) = m_EV.col(ind_max);

            // store it
            m_kExpVarV(k) = val_max;
        }

        // return
        for(int k=0; k<m_K; k++) {
            order(k) = chosenIndex[k] + 1;
        }
    }

    /*!
     * \brief order factors according to deviance
     *
     * @param[out] vector of factor order
     */
    void gamPoisFactorStandard::orderDeviance(VectorXi &order) {
        vector<int> leftIndex;
        vector<int> chosenIndex;

        MatrixXd Utmp = MatrixXd::Zero(m_EU.rows(), m_EU.cols());
        MatrixXd Vtmp = MatrixXd::Zero(m_EV.rows(), m_EV.cols());

        // init
        for(int k=0; k<m_K; k++) {
            leftIndex.push_back(k);
        }

        // ordering
        for(int k=0; k<m_K; k++) {
            double val_min = 0;
            int ind_min = 0;
            int ind_left = 0;
            for(int ind=0; ind<leftIndex.size(); ind++) {

                MatrixXd tmpU(m_EU.rows(), k+1);
                MatrixXd tmpV(m_EV.rows(), k+1);

                tmpU.col(k) = m_EU.col(leftIndex[ind]);
                tmpV.col(k) = m_EV.col(leftIndex[ind]);

                if(k!=0) {
                    tmpU.leftCols(k) = Utmp.leftCols(k);
                    tmpV.leftCols(k) = Vtmp.leftCols(k);
                }

                MatrixXd lambda = tmpU * tmpV.transpose();

                Rcpp::Rcout << "min lambda = " << lambda.minCoeff() << std::endl;
                Rcpp::Rcout << "max lambda = " << lambda.maxCoeff() << std::endl;

                Rcpp::Rcout << "min lambda0 = " << m_lambda0.minCoeff() << std::endl;
                Rcpp::Rcout << "max lambda0 = " << m_lambda0.maxCoeff() << std::endl;

                double res = poisDeviance(m_X, lambda, m_lambda0);
                if(ind==0) {
                    val_min = res;
                    ind_min = leftIndex[ind];
                    ind_left = ind;
                }
                if(res < val_min) {
                    val_min = res;
                    ind_min = leftIndex[ind];
                    ind_left = ind;
                }
            }
            // max found
            chosenIndex.push_back(ind_min);
            leftIndex.erase(leftIndex.begin() + ind_left);

            Utmp.col(k) = m_EU.col(ind_min);
            Vtmp.col(k) = m_EV.col(ind_min);

            // store it
            m_kDeviance(k) = val_min;

        }

        // return
        for(int k=0; k<m_K; k++) {
            order(k) = chosenIndex[k] + 1;
        }
    }

    /*!
     * \brief compute factor order
     */
    void gamPoisFactorStandard::computeOrder() {
        this->orderExpVar0(m_orderExpVar0);
        this->orderExpVarU(m_orderExpVarU);
        this->orderExpVarV(m_orderExpVarV);
        this->orderDeviance(m_orderDeviance);
    }


    //-------------------//
    //       return      //
    //-------------------//

    /*!
     * \brief create list with results to be return
     *
     * @param[out] list containing output
     */
    void gamPoisFactorStandard::returnObject(Rcpp::List &results) {
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
                                              Rcpp::Named("orderExpVarV") = m_orderExpVarV);

       Rcpp::List criteria_k = Rcpp::List::create(Rcpp::Named("kDeviance") = m_kDeviance,
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
                                                  Rcpp::Named("criteria_k") = criteria_k,
                                                  Rcpp::Named("normGap") = m_normGap.head(m_nbIter),
                                                  Rcpp::Named("deviance") = m_deviance.head(m_nbIter),
                                                  Rcpp::Named("converged") = m_converged,
                                                  Rcpp::Named("nbIter") = m_nbIter);

        SEXP tmp = Rcpp::Language("c", results, returnObj).eval();

        results = tmp;
    }
}