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
#include "gamDistrib.h"
#include "intermediate.h"

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
    gamPoisFactorPen::gamPoisFactorPen(int n, int p, int K, const MatrixXi &X,
                                       const MatrixXd &phi1, const MatrixXd &phi2,
                                       const MatrixXd &theta1, const MatrixXd &theta2,
                                       const MatrixXd &alpha1, const MatrixXd &alpha2,
                                       const MatrixXd &beta1, const MatrixXd &beta2,
                                       const VectorXd &r_theta2, const VectorXd &r_phi2)
    : gamPoisFactor::gamPoisFactor(n, p, K, X, phi1, phi2, theta1, theta2,
                                    alpha1, alpha2, beta1, beta2)
    {
        // penalty constant
        m_r_theta2 = VectorXd(r_theta2);
        m_r_phi2 = VectorXd(r_phi2);

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
                estimParam(1000, m_alpha1cur(i,k), m_alpha2cur(i,k), param1, param2);
                m_phi1cur(i,k) = param1;
                m_phi1old(i,k) = param1;
                m_phi2inter(i,k) = param2;
            }
            // global parameters
            for(int j=0; j<m_P; j++) {
                double param1 = 0;
                double param2 = 0;
                estimParam(1000, m_beta1cur(j,k), m_beta2cur(j,k), param1, param2);
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
        m_phi1cur = m_alpha1cur.array() + m_EZ_j.array();
        m_phi2inter = m_alpha2cur.rowwise() + m_EV.colwise().sum();

        this->penLocalParam();

        // expectation and log-expectation
        Egam(m_phi1cur, m_phi2cur, m_EU);
        Elgam(m_phi1cur, m_phi2cur, m_ElogU);
    }

    /*!
    * \brief update rule for penalized local parameters phi (factor U) in variational inference
    */
    void gamPoisFactorPen::penLocalParam() {
        m_phi2cur = (m_phi2inter.array() * m_phi1cur.array() + ( (m_phi2inter.array() * m_phi1cur.array()).msquare()
                        + 8 * (m_phi1cur.array().rowwise() * m_r_phi2.transpose().array())).msqrt()).array() / (2 * m_phi1cur.array());

        // test
        //for(int i=0; i<m_N; i++) {
        //    for(int k = 0; k<m_K; k++) {
        //        double test = ( m_phi2inter(i,k) * m_phi1cur(i,k) + sqrt( std::pow(m_phi2inter(i,k) * m_phi1cur(i,k),2)
        //                        + 8 * m_phi1cur(i,k) * m_r_phi2(k)) )/(2*m_phi1cur(i,k));
        //        if(test != m_phi2cur(i,k)) {
        //            Rcpp::Rcout << "error in penalized updates" << std::endl;
        //        }
        //    }
        //}
    }

    /*!
    * \brief update rule for global parameters theta (factor V) in variational inference
    */
    void gamPoisFactorPen::globalParam() {
        m_theta1cur = m_beta1cur.array() + m_EZ_i.array();
        m_theta2inter = m_beta2cur.rowwise() + m_EU.colwise().sum();

        this->penGlobalParam();

        // expectation and log-expectation
        Egam(m_theta1cur, m_theta2cur, m_EV);
        Elgam(m_theta1cur, m_theta2cur, m_ElogV);
    }

    /*!
     * \brief update rule for penalized global parameters theta (factor V) in variational inference
     */
    void gamPoisFactorPen::penGlobalParam() {
        m_theta2cur = (m_theta2inter.array() * m_theta1cur.array() + ( (m_theta2inter.array() * m_theta1cur.array()).msquare()
                    + 8 * (m_theta1cur.array().rowwise() * m_r_theta2.transpose().array())).msqrt()).array() / (2 * m_theta1cur.array());
    }

    //-------------------//
    //       return      //
    //-------------------//

    /*!
    * \brief create list with results to be return
    *
    * @param[out] list containing output
    */
    void gamPoisFactorPen::returnObject(Rcpp::List &results) {

        Rcpp::List pen = Rcpp::List::create(Rcpp::Named("r_theta2") = m_r_theta2,
                                            Rcpp::Named("r_phi2") = m_r_phi2);

        Rcpp::List returnObj;
        gamPoisFactor::returnObject(returnObj);

        SEXP tmp1 = Rcpp::Language("c", returnObj, pen).eval();

        SEXP tmp2 = Rcpp::Language("c", results, tmp1).eval();

        results = tmp2;
    }

}