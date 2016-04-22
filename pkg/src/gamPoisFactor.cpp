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

/*!
 * \file gamPoisFactor.cpp
 * \brief class definition for Gamma Poisson Factor Model
 * \author Ghislain Durif
 * \version 0.1
 * \date 22/04/2016
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/math/special_functions/digamma.hpp>
#include "gamPoisFactor.h"

using namespace Rcpp;

#define coeffExp() unaryExpr(std::ptr_fun<double,double>(exp))
#define coeffSum(X) unaryExpr(std::bind2nd(std::pointer_to_binary_function<double,double,double>(scalsum),X))
#define square() unaryExpr(std::bind2nd(std::pointer_to_binary_function<double,double,double>(pow),2))

// [[Rcpp::depends(BH)]]
using boost::math::digamma;

// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision

// SRC
namespace countMatrixFactor {

    // CONSTRUCTOR
    gamPoisFactor::gamPoisFactor(int n, int p, int K, int iterMax, int order,
                                 int stabRange, double epsilon, bool verbose,
                                 const MatrixXi &X,
                                 const MatrixXd &phi1, const MatrixXd &phi2,
                                 const MatrixXd &theta1, const MatrixXd &theta2,
                                 const MatrixXd &alpha1, const MatrixXd &alpha2,
                                 const MatrixXd &beta1, const MatrixXd &beta2) {

        // dimensions
        m_N = X.rows();
        m_P = X.cols();
        m_K = K;

        // parameters
        m_iterMax = iterMax;
        m_order = order;
        m_stabRange = stabRange;
        m_epsilon = epsilon;
        m_verbose = verbose;

        // data
        m_X = MatrixXi(X);

        // variational parameters
        m_phi1cur = MatrixXd::Zero(n,K);
        m_phi2cur = MatrixXd::Zero(n,K);

        m_phi1old = MatrixXd(phi1);
        m_phi2old = MatrixXd(phi2);

        m_theta1cur = MatrixXd::Zero(p,K);
        m_theta2cur = MatrixXd::Zero(p,K);

        m_theta1old = MatrixXd(theta1);
        m_theta2old = MatrixXd(theta2);

        // sufficient statistics
        m_EU = MatrixXd::Zero(n,K);
        m_ElogU = MatrixXd::Zero(n,K);
        m_EV = MatrixXd::Zero(p,K);
        m_ElogV = MatrixXd::Zero(p,K);

        m_EZ_i = MatrixXd::Zero(p,K);
        m_EZ_j = MatrixXd::Zero(n,K);

        // prior parameter
        m_alpha1 = MatrixXd(alpha1);
        m_alpha2 = MatrixXd(alpha2);
        m_beta1 = MatrixXd(beta1);
        m_beta2 = MatrixXd(beta2);

        // criterion
        m_margLogLike = VectorXd::Zero(iterMax);
        m_condLogLike = VectorXd::Zero(iterMax);
        m_priorLogLike = VectorXd::Zero(iterMax);
        m_postLogLike = VectorXd::Zero(iterMax);
        m_compLogLike = VectorXd::Zero(iterMax);

        m_deviance = VectorXd::Zero(iterMax);
        m_normGap = VectorXd::Zero(iterMax);
        m_expVar0 = VectorXd::Zero(iterMax);
        m_expVarU = VectorXd::Zero(iterMax);
        m_expVarV = VectorXd::Zero(iterMax);

    }

    // DESTRUCTOR
    gamPoisFactor::~gamPoisFactor() {}

    // member functions: documented in src

    // initialization
    /*!
     * \brief Initialization of sufficient statistics
     */
    void gamPoisFactor::Init() {

        this->Egamma(m_alpha1, m_alpha2, m_EU);
        this->Elgamma(m_alpha1, m_alpha2, m_ElogU);
        this->Egamma(m_beta1, m_beta2, m_EV);
        this->Elgamma(m_beta1, m_beta2, m_ElogV);

    }

    //-------------------//
    // sufficient stats  //
    //-------------------//

    // Expectation gamma
    void Egamma(const MatrixXd &alpha, const MatrixXd &beta, MatrixXd &res) {

    }

    // Expectation log gamma
    void Elgamma(const MatrixXd &alpha, const MatrixXd &beta, MatrixXd &res) {

    }

    //-------------------//
    //   convergence     //
    //-------------------//

    // parameter norm
    double parameterNorm(const MatrixXd &phi1, const MatrixXd &phi2,
                         const MatrixXd &theta1, const MatrixXd &theta2) {

    }

    // difference norm (on parameters)
    double differenceNorm(const MatrixXd &phi1old, const MatrixXd &phi2old,
                          const MatrixXd &theta1old, const MatrixXd &theta2old,
                          const MatrixXd &phi1new, const MatrixXd &phi2new,
                          const MatrixXd &theta1new, const MatrixXd &theta2new) {

    }

    // convergence condition
    double convCondition(int order, const VectorXd &normGap, int iter, int drift) {

    }

    //-------------------//
    // parameter updates //
    //-------------------//

    // Poisson intensity
    void poisRate(int n, int p, int K,
                  MatrixXd &EZ_i, MatrixXd &EZ_j,
                  const MatrixXd &ElogU, const MatrixXd &ElogV);

    // local parameters: phi (factor U)
    void localParam(int n, int p, int K,
                    MatrixXd &phi1cur, MatrixXd &phi2cur,
                    const MatrixXd &EZ_j, const MatrixXi &X, const MatrixXd &EV,
                    const MatrixXd &alpha1, const MatrixXd &alpha2);

    // global parameters: theta (factor V)
    void globalParam(int n, int p, int K,
                     MatrixXd &theta1cur, MatrixXd &theta2cur,
                     const MatrixXd &EZ_j, const MatrixXi &X, const MatrixXd &EU,
                     const MatrixXd &beta1, const MatrixXd &beta2);






}