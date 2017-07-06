### Copyright 2017-02 Ghislain DURIF
###
### This file is part of the `pCMF' library for R and related languages.
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


#' @title pCMF
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
#' @useDynLib pCMF
#'
#' @param X matrix n x p of counts
#' @param K number of factors
#' @param alpha1 n x K, initial values of first parameter of Gamma prior on U
#' @param alpha2 n x K, initial values of second parameter of Gamma prior on U
#' @param beta1 n x K, initial values of first parameter of Gamma prior on V
#' @param beta2 n x K, initial values of second parameter of Gamma prior on V
#'
#' @return return
#' \item{Y}{Y}
#'
#' @export
pCMF <- function(X, K, alpha1=NULL, alpha2=NULL, beta1=NULL, beta2=NULL,
                iterMax=200, iterMin=100, epsilon=1e-5,
                verbose=TRUE, sparse=FALSE, ZI=FALSE, ncores=1,
                nbInit=1, iterMaxInit=50, noise=0.5, seed=NULL) {

    ncomp <- K
    n <- nrow(X)
    p <- ncol(X)

    if(is.null(alpha1)) {
        alpha01 <- matrix(1, nrow=n, ncol=ncomp)
        alpha02 <- matrix(1, nrow=n, ncol=ncomp)

        beta01 <- matrix(1, nrow=p, ncol=ncomp)
        beta02 <- matrix(1, nrow=p, ncol=ncomp)
    }

    phi01 <- matrix(1, nrow=n, ncol=ncomp)
    phi02 <- matrix(1, nrow=n, ncol=ncomp)

    theta01 <- matrix(1, nrow=p, ncol=ncomp)
    theta02 <- matrix(1, nrow=p, ncol=ncomp)

    results <- matrixFactor(X=X, K=ncomp,
                            phi01=phi01, phi02=phi02, theta01=theta01, theta02=theta02,
                            alpha1=alpha01, alpha2=alpha02, beta1=beta01, beta2=beta02,
                            iterMax=iterMax, iterMin=iterMin, epsilon=1e-5,
                            order=0, stabRange=5, verbose=verbose, sparse=sparse, ZI=ZI,
                            algo="EM", ncores=ncores)

    class(results) = "pCMF"

    return(results)

}
