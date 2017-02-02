### Copyright 2016-04 Ghislain DURIF
###
### This file is part of the `cmf' library for R and related languages.
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
#' R wrapper for Gamma-Poisson Factor model
#'
#' @details
#' Factorization of count matrix
#'
#' @author
#' Ghislain Durif, \email{gd.dev@libertymail.net}
#'
#'
#' @seealso aaa
#'
#' @importFrom Rcpp evalCpp
#' @useDynLib cmf
#'
#' @param X matrix n x p of counts
#' @param K number of factors
#' @param phi01 n x K, initial values of first parameter of Gamma distribution on U
#' @param phi02 n x K, initial values of second parameter of Gamma distribution on U
#' @param theta01 n x K, initial values of first parameter of Gamma distribution on V
#' @param theta02 n x K, initial values of second parameter of Gamma distribution on V
#' @param alpha1 n x K, initial values of first parameter of Gamma prior on U
#' @param alpha2 n x K, initial values of second parameter of Gamma prior on U
#' @param beta1 n x K, initial values of first parameter of Gamma prior on V
#' @param beta2 n x K, initial values of second parameter of Gamma prior on V
#'
#' @return return
#' \item{Y}{Y}
#'
#' @export
matrixFactor <- function(X, K, phi01, phi02, theta01, theta02,
                         alpha1, alpha2, beta1, beta2,
                         lambda = NULL, mu = NULL,
                         iterMax=200, iterMin=1, epsilon=1e-5,
                         order=0, stabRange=5, verbose=TRUE, sparse=FALSE, ZI=FALSE,
                         algo="EM") {

    X = apply(X, c(1,2), as.integer)

    results = NULL

    if(algo == "EM") {
        # print("EM ok")
        if(ZI) {
            if(sparse) {
                results = gamPoisFactorEM_Sparse_ZI_wrapper(X, K, ZI, phi01, phi02, theta01, theta02,
                                                            alpha1, alpha2, beta1, beta2,
                                                            iterMax, iterMin, epsilon,
                                                            order, stabRange, verbose)
            } else {
                results = gamPoisFactorEM_ZI_wrapper(X, K, ZI, phi01, phi02, theta01, theta02,
                                                     alpha1, alpha2, beta1, beta2,
                                                     iterMax, iterMin, epsilon,
                                                     order, stabRange, verbose)
            }
        } else {
            if(sparse) {
                results = gamPoisFactorEM_Sparse_wrapper(X, K, ZI, phi01, phi02, theta01, theta02,
                                                         alpha1, alpha2, beta1, beta2,
                                                         iterMax, iterMin, epsilon,
                                                         order, stabRange, verbose)
            } else {
                results = gamPoisFactorEM_wrapper(X, K, ZI, phi01, phi02, theta01, theta02,
                                                  alpha1, alpha2, beta1, beta2,
                                                  iterMax, iterMin, epsilon,
                                                  order, stabRange, verbose)
            }
        }
    }

    if(algo == "variational") {
        if(ZI) {
            results = gamPoisFactorZI_wrapper(X, K, ZI, phi01, phi02, theta01, theta02,
                                              alpha1, alpha2, beta1, beta2,
                                              iterMax, iterMin, epsilon,
                                              order, stabRange, verbose)
        } else {
            if(sparse) {
                results = gamPoisFactorSparse_wrapper(X, K, ZI, phi01, phi02, theta01, theta02,
                                                      alpha1, alpha2, beta1, beta2,
                                                      iterMax, iterMin, epsilon,
                                                      order, stabRange, verbose)
            } else {
                results = gamPoisFactor_wrapper(X, K, ZI, phi01, phi02, theta01, theta02,
                                                alpha1, alpha2, beta1, beta2,
                                                iterMax, iterMin, epsilon,
                                                order, stabRange, verbose)
            }
        }
    }

    class(results) = "cmf"

    return(results)

}
