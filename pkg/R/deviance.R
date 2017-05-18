### Copyright 2017-05 Ghislain DURIF
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


#' @title expDev
#'
#' @description
#' Get the percentage of deviance explained by a Gamma-Poisson Factor model
#'
#' @details
#' see pCMF function output
#'
#' @author
#' Ghislain Durif, \email{gd.dev@libertymail.net}
#'
#'
#' @seealso \code{\link{pCMF}}
#'
#' @useDynLib pCMF
#'
#' @param model a Gamma-Poisson factor model output by the
#' function \code{\link{pCMF}}
#'
#' @return the matrix U of individuals coordinates in the lower
#' dimensional sub-space
#'
#' @export
expDev <- function(model, X) {

    ## computation of the deviance percentage

    n <- nrow(X)
    p <- ncol(X)

    lambda <- model$lambda

    lambda0 <- matrix(rep(apply(X,2,mean), times=n),nrow=n, ncol=p, byrow=TRUE)

    res1 = poisLogLike(X, lambda)
    res2 = poisLogLike(X, X)
    res0 = poisLogLike(X, lambda0)

    return((res1-res0)/(res2-res0))
}


#' @title dev
#'
#' @description
#' Get the associated to a Gamma-Poisson Factor model
#'
#' @details
#'
#'
#' @author
#' Ghislain Durif, \email{gd.dev@libertymail.net}
#'
#'
#' @seealso \code{\link{expDev}}
#'
#' @useDynLib pCMF
#'
#' @param model a Gamma-Poisson factor model output by the
#' function \code{\link{pCMF}}
#'
#' @return the matrix U of individuals coordinates in the lower
#' dimensional sub-space
#'
#' @export
dev <- function(model, X) {

    ## computation of the deviance

    lambda <- model$lambda

    res1 = poisLogLike(X, lambda)
    res2 = poisLogLike(X, X)

    return(-2*(res1 - res2))
}


#' @title poisLogLike
#' @keywords internal
#'
#' @description
#' Poisson log-likelihood
#'
#' @details
#'
#'
#' @author
#' Ghislain Durif, \email{gd.dev@libertymail.net}
#'
#'
#' @useDynLib pCMF
#'
#' @param X matrix n x p of count data
#' @param lambda matrix n x p of Poisson rate
#'
#' @return value of the Poisson log-likelihood
#'
#' @export
poisLogLike = function(X, lambda) {
    X.ij = as.vector(X)
    lambda.ij = as.vector(lambda)
    lambda.ij.nnull = ifelse(lambda.ij<1E-12, 1, lambda.ij)

    return(sum(X.ij * log(lambda.ij.nnull) - lambda.ij - lfactorial(X.ij)))
}
