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
 * \file intermediate.cpp
 * \brief intermediate functions
 * \author Ghislain Durif
 * \version 0.1
 * \date 25/04/2016
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <cstdio>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include "intermediate.h"

#define msquare() unaryExpr(std::bind2nd(std::pointer_to_binary_function<double,double,double>(std::pow),2))

// [[Rcpp::depends(BH)]]
using boost::math::digamma;
using boost::math::trigamma;
using boost::math::binomial_coefficient;
using boost::math::factorial;

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;                  // variable size matrix, double precision

namespace intermediate {

    /*!
     * \fn Dirac function in zero
     *
     * Indicate if the considered value is null
     *
     * @param x value tested
     * @return 0 or 1 indicating if x==0 or not
     *
     * @TODO verify the 0-1 return
     */
    double dirac(double x) {
        if(x==0) {
            return 1;
        } else {
            return 0;
        }
    }

    /*!
     * \fn exponential that throw error when under or overflowing
     *
     * @param[in] x real
     * @return the value of exp(x)
     */
    double myExp(double x) {
        if(std::abs(x) > 700) {
            Rcpp::stop("under or overflow in exponential");
            return(-1);
        } else {
            return exp(x);
        }
    }

    double safeExp(double x) {
        if(x > 700) {
            return(1E12);
        } else if(x < -700) {
            return(1E-12);
        } else {
            return std::exp(x);
        }
    }

    /*!
     * \fn logit function
     *
     * log(x/(1-x))
     *
     * @param[in] x real between 0 and 1
     * @return the value of logit(x)
     */
    double logit(double x) {
        return std::log(x/(1-x));
    }

    /*!
     * \fn logit inverse function
     *
     * exp(x)/(1+exp(x)) = 1/(1+exp(-x))
     *
     * @param[in] x real between 0 and 1
     * @return the value of logit(x)
     */
    double expit(double x) {
        if (x > 30) return 1;
        if (x < -30) return 0;
        return 1 / (1 + std::exp(-x));
    }

    //
    /*!
     * \fn compute E[log(Z!)] where Z follows a binomial distribution
     *
     * Binomial distribution B(N,p)
     *
     * @param[in] N nb of tries in binomial distribution B(N,p)
     * @param[in] p values of probability
     * @return real value
     */
    double lgamBinom(int N, double p) {
        double res=0;
        if(N < 50) {
            for(int l=0; l < N+1; l++) {
                res += std::log(factorial<double>((double) l)) * binomial_coefficient<double>((double) N, (double) l) * std::pow(p,l) * std::pow(1-p, N-l);
            }
        } else {
            res = lgamma(N*p + 1);
        }

        return res;
    }


    /*!
     * \fn erase zeros in matrix
     *
     * @param[in,out] X matrix of data (count)
     */
    void eraseZero(MatrixXd &X) {
        int n = X.rows();
        int p = X.cols();
        for(int i=0; i<n; i++) {
            for(int j=0; j<p; j++) {
                if(X(i,j) == 0) {
                    X(i,j) = X(i,j) + 1;
                }
            }
        }
    }

    /*!
     * \fn threshold at t
     *
     * @param[in] x input value
     * @param[in] t threshold
     * @return max(x,t)
     */
    double threshold(double x, double t) {
        return(std::max(x,t));
    }

    /*!
     * \fn check if exponential can be applied to a matrix
     *
     * throw error when potential under or overflowing
     *
     * @param[in] A matrix to be checked
     */
    void checkExp(const MatrixXd &A) {
        // check for 300 and -300 so a product of such elements does not exceed the limit
        if( A.maxCoeff() > 300) {
            Rcpp::stop("overflow in element-wise exponentiation of matrix");
        }
        if( A.minCoeff() < - 300) {
            Rcpp::stop("underflow in element-wise exponentiation of matrix");
        }
    }

    /*!
     * \fn check matrix (dim, min and max element)
     *
     * @param[in] A matrix to be checked
     */
    void checkMat(const MatrixXd &A) {
        Rcpp::Rcout << "dimension: " << A.rows() << " x " << A.cols() << " / min = " << A.minCoeff() << " / max = " << A.maxCoeff() << std::endl;
    }

    /*!
     * \fn inverse digamma (psi) function
     *
     * find the solution in x with y given of y = psi(x)
     *
     * @param[in] y the value of the digamma function
     * @return the corresponding x
     */
    double psiInv(double y, int nbIter) {
        double x0 = 0;
        double x = 0;

        // init
        if(y >= -2.22) {
            x0 = std::exp(y) + 0.5;
        } else {
            x0 = -1/(y-digamma(1));
        }

        // iter
        for(int i=0; i<nbIter; i++) {
            x = x0 - (digamma(x0) - y)/trigamma(x0);
            x0 = x;
        }

        return x;
    }

    /*!
     * \brief l2 squared norm of all parameters
     *
     * Computation of sum_{ij} param1_{ij}^2 + param2_{ij}^2
     *
     * @param[in] param1 rows x cols, matrix of first parameters
     * @param[in] param2 rows x cols, matrix of second parameters
     * @return res l2 squared norm
     */
    double parameterNorm2(const MatrixXd &param1, const MatrixXd &param2) {
        double res = param1.msquare().sum() + param2.msquare().sum();
        return res;
    }

    /*!
     * \fn difference of squared euclidean norm (on parameters)
     *
     * @param[in] param1a matrix of parameters 1 (state a)
     * @param[in] param2a matrix of parameters 2 (state a)
     * @param[in] param1b matrix of parameters 1 (state b)
     * @param[in] param2b matrix of parameters 2 (state b)
     *
     * @return sum((param1a - param1b)^2) + sum((param2a-param2b)^2)
     */
    double differenceNorm2(const MatrixXd &param1a, const MatrixXd &param2a, const MatrixXd &param1b, const MatrixXd &param2b) {
        double res = parameterNorm2(param1a - param1b, param2a - param2b);
        return(res);
    }

}
