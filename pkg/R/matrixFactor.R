### Copyright 2016-04 Ghislain DURIF
###
### This file is part of the `countMatrixFactor' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


#' @title matrixFactor
#'
#' @description
#' Description
#'
#' @details
#' Wrapper for Cpp function
#'
#' @author
#' Ghislain Durif, \email{ghislain.durif@univ-lyon1.fr}
#'
#'
#' @seealso
#'
#' @import Rcpp
#' @import RcppEigen
#' @useDynLib countMatrixFactor
#'
#' @param
#'
#' @return return
#' \item{}{}
#'

### R wrapper for Gamma-Poisson Factor model

#' @export
matrixFactor = function(X, K, phi01, phi02, theta01, theta02,
                        alpha1, alpha2, beta1, beta2,
                        lambda = NULL, mu = NULL,
                        iterMax=200, iterMax_Estep=50, iterMax_Mstep=50, epsilon=1e-5,
                        order=0, stabRange=5, verbose=TRUE, pen=FALSE, sparse=FALSE, ZI=FALSE,
                        algo="EM") {

    X = apply(X, c(1,2), as.integer)

    results = NULL

    if(algo == "EM") {
        print("EM ok")
        results = gamPoisFactorEM_wrapper(X, K, phi01, phi02, theta01, theta02,
                                          alpha1, alpha2, beta1, beta2,
                                          iterMax, iterMax_Estep, iterMax_Mstep ,epsilon,
                                          order, stabRange, verbose)
    }

    if(algo == "variational") {
        if(ZI) {
            results = gamPoisFactorZI_wrapper(X, K, phi01, phi02, theta01, theta02,
                                              alpha1, alpha2, beta1, beta2,
                                              iterMax, epsilon,
                                              order, stabRange, verbose)
        } else {
            if(pen) {
                if(sparse) {
                    results = gamPoisFactorSparse_wrapper(X, K, phi01, phi02, theta01, theta02,
                                                          alpha1, alpha2, beta1, beta2,
                                                          lambda, mu,
                                                          iterMax, epsilon,
                                                          order, stabRange, verbose)
                } else {
                    results = gamPoisFactorPen_wrapper(X, K, phi01, phi02, theta01, theta02,
                                                       alpha1, alpha2, beta1, beta2,
                                                       lambda, mu,
                                                       iterMax, epsilon,
                                                       order, stabRange, verbose)
                }
            } else {
                results = gamPoisFactor_wrapper(X, K, phi01, phi02, theta01, theta02,
                                                alpha1, alpha2, beta1, beta2,
                                                iterMax, epsilon,
                                                order, stabRange, verbose)
            }
        }
    }

    class(results) = "countMatrixFactor"

    return(results)

}