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
 * \file gamPoisFactorStandard.cpp
 * \brief class definition for standard Gamma Poisson Factor Model
 * \author Ghislain Durif
 * \version 0.1
 * \date 04/05/2016
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/math/special_functions/digamma.hpp>
#include "gamPoisFactorStandard.h"

using namespace Rcpp;

#define coeffExp() unaryExpr(std::ptr_fun<double,double>(exp))
#define coeffSum(X) unaryExpr(std::bind2nd(std::pointer_to_binary_function<double,double,double>(scalsum),X))
#define square() unaryExpr(std::bind2nd(std::pointer_to_binary_function<double,double,double>(pow),2))
#define digamma() unaryExpr(std::ptr_fun<double,double>(digamma))
#define log() unaryExpr(std::ptr_fun<double,double>(log))

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
    gamPoisFactorStandard::gamPoisFactorStandard(int n, int p, int K, int iterMax, int order,
                                                 int stabRange, double epsilon, bool verbose,
                                                 const MatrixXi &X,
                                                 const MatrixXd &phi1, const MatrixXd &phi2,
                                                 const MatrixXd &theta1, const MatrixXd &theta2,
                                                 const MatrixXd &alpha1, const MatrixXd &alpha2,
                                                 const MatrixXd &beta1, const MatrixXd &beta2) :
                gamPoisFactor::gamPoisFactor(n, p, K, iterMax, order,
                                              stabRange, epsilon, verbose,
                                              X, phi1, phi2, theta1, theta2,
                                              alpha1, alpha2, beta1, beta2) {}

    // DESTRUCTOR
    gamPoisFactorStandard::~gamPoisFactorStandard() {}

    // member functions: documented in src

    /*!
     * \brief Initialization of sufficient statistics
     */
    void gamPoisFactor::Init() {}

    //-------------------//
    // parameter updates //
    //-------------------//

    // Poisson intensity
    void gamPoisFactor::poisRate() {

    }

    // local parameters: phi (factor U)
    void gamPoisFactor::localParam() {
//         for(int i=0; i<n; i++) {
//             for(int k=0; k<K; k++) {
//                 // sum_j X[i,j] * omega[i,j,k]
//                 double res=0;
//                 for(int j=0; j<p; j++) {
//                     res += (double) X(i,j) * omega[i][j][k];
//                 }
//                 phi1cur(i,k) = alpha1(i,k) + res;
//                 phi2cur(i,k) = alpha2(i,k) + EV.col(k).sum();
//             }
//         }
    }

    // global parameters: theta (factor V)
    void gamPoisFactor::globalParam() {

    }


    //-------------------//
    //     algorithm     //
    //-------------------//
    void algorithm() {

    }



}