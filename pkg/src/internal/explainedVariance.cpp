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
* \file loglikelihood.cpp
* \brief class definition for log-likelihood
* \author Ghislain Durif
* \version 0.1
* \date 22/04/2016
*/

#include <Rcpp.h>
#include <RcppEigen.h>
#include "model/explainedVariance.h"

#define msquare() unaryExpr(std::bind2nd(std::pointer_to_binary_function<double,double,double>(std::pow),2))

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;              // variable size matrix, double precision
using Eigen::MatrixXi;              // variable size matrix, integer
using Eigen::VectorXd;              // variable size matrix, double precision
using Eigen::PartialPivLU;              // for fast matrix inversion

namespace countMatrixFactor {

    //------------------------------------------------------------------------//
    // CONSTRUCTOR
    explainedVariance::explainedVariance(int size) {
        m_expVar0 = VectorXd::Zero(size);
        m_expVarU = VectorXd::Zero(size);
        m_expVarV = VectorXd::Zero(size);
    }

    // DESTRCTOR
    explainedVariance::~explainedVariance() {}

    //------------------------------------------------------------------------//
    // GETTER

    /*!
     * \brief getter for expVar0
     *
     * @param[out] res vector of explained variance as residual sum of squares
     */
    void explainedVariance::getExpVar0(VectorXd &res) {
        res = m_expVar0;
    }

    /*!
     * \brief getter for expVarU
     *
     * @param[out] res vector of explained variance by columns of U
     */
    void explainedVariance::getExpVarU(VectorXd &res) {
        res = m_expVarU;
    }

    /*!
     * \brief getter for expVarV
     *
     * @param[out] res vector of explained variance by columns of V
     */
    void explainedVariance::getExpVarV(VectorXd &res) {
        res = m_expVarV;
    }


    // FUNCTIONS
    /*!
     * \fn computes explained variance as residual sum of squares
     *
     * explained variance as sum_{i,j} (X_{i,j} - (UtV)_{i,j})^2 / sum_{i,j} X_{i,j}^2
     *
     * @param[in] X matrix of count data
     * @param[in] U matrix of observation coordinates in latent subspace
     * @param[in] V matrix of variable coordinates in latent subspace
     *
     * @return effective value
     *
     * @TODO export in R
     */
    double expVar0(const MatrixXi &X, const MatrixXd &U, const MatrixXd &V) {
        double res;
        res = (double) 1 - ( X.cast<double>().array() - (U * V.transpose()).array()).msquare().sum() / ( X.cast<double>().msquare().sum() );
        return res;
    }

    /*!
     * \fn computes explained variance by columns of U
     *
     * explained variance as t(X.proj) %*% X.proj / t(X) %*% X,
     * where X.proj is the projection of X onto the new components
     * Xproj = U %*% ginv(t(U) %*% U) %*% t(U) %*% X
     *
     * @param[in] X matrix of count data
     * @param[in] U matrix of observation coordinates in latent subspace
     *
     * @return effective value
     *
     * @TODO export in R
     */
    double expVarU(const MatrixXi &X, const MatrixXd &U) {
        MatrixXd Xproj(X.rows(), X.cols());
        MatrixXd Xcent(X.rows(), X.cols());
        MatrixXd Xprojcent(X.rows(), X.cols());
        double res;
        MatrixXd prod(U.transpose() * U + 0.00001 * MatrixXd::Identity(U.cols(), U.cols()));
        Xproj = U * PartialPivLU<MatrixXd>(prod).inverse() * U.transpose() * X.cast<double>();
        Xcent = X.cast<double>().rowwise() - X.cast<double>().colwise().mean();
        Xprojcent = Xproj.rowwise() - Xproj.colwise().mean();
        res = (Xprojcent.transpose() * Xprojcent).trace() / (Xcent.transpose() * Xcent).trace();
        return res;
    }

    /*!
     * \fn computes explained variance by columns of V
     *
     * explained variance as t(X.proj) %*% X.proj / t(X) %*% X,
     * where X.proj is the projection of X onto the new components
     * Xproj = X %*% V %*% ginv(t(V) %*% V) %*% t(V)
     *
     * @param[in] X matrix of count data
     * @param[in] V matrix of variable coordinates in latent subspace
     *
     * @return effective value
     *
     * @TODO export in R
     */
    double expVarV(const MatrixXi &X, const MatrixXd &V) {
        MatrixXd Xproj(X.rows(), X.cols());
        MatrixXd Xcent(X.rows(), X.cols());
        MatrixXd Xprojcent(X.rows(), X.cols());
        double res;
        MatrixXd prod(V.transpose() * V + 0.00001 * MatrixXd::Identity(V.cols(), V.cols()));
        Xproj = X.cast<double>() * V * PartialPivLU<MatrixXd>(prod).inverse() * V.transpose();
        Xcent = X.cast<double>().rowwise() - X.cast<double>().colwise().mean();
        Xprojcent = Xproj.rowwise() - Xproj.colwise().mean();
        res = (Xprojcent.transpose() * Xprojcent).trace() / (Xcent.transpose() * Xcent).trace();
        return res;
    }
}










