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


#' @title getU
#'
#' @description
#' Getter for the matrix U in the Gamma-Poisson Factor model
#'
#' @details
#' Rrepresentation of individuals in in the lower
#' dimensional sub-space in the Euclidean geometry, corresponding to the
#' matrix U, or in the geometry related to the Gamma distribution in the
#' exponential family (log), corresponding to matrix logU
#'
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
#' @param log_representation boolean, indicating if the representation is in
#' the natural geometry associated to the Gamma distribution (log) or in the
#' Euclidean space, default is TRUE
#'
#' @return the matrix U of individuals coordinates in the lower
#' dimensional sub-space
#'
#' @export
getU <- function(model, log_representation=TRUE) {

    if(class(model) != "pCMF")
        stop("wrong model in input")

    if(log_representation) {
        U <- as.matrix(model$stats$ElogU[,model$order$orderDeviance])
    } else {
        U <- as.matrix(model$U[,model$order$orderDeviance])
    }

    return(U)

}


#' @title getV
#'
#' @description
#' Getter for the matrix V in the Gamma-Poisson Factor model
#'
#' @details
#' Rrepresentation of individuals in in the lower
#' dimensional sub-space in the Euclidean geometry, corresponding to the
#' matrix U, or in the geometry related to the Gamma distribution in the
#' exponential family (log), corresponding to matrix logU
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
#' @param log_representation boolean, indicating if the representation is in
#' the natural geometry associated to the Gamma distribution (log) or in the
#' Euclidean space, default is TRUE
#'
#'
#' @return the matrix V of variables \(features\) contributions to the lower
#' dimensional sub-space
#'
#' @export
getV <- function(model, log_representation=TRUE) {

    if(class(model) != "pCMF")
        stop("wrong model in input")

    if(log_representation) {
        V <- as.matrix(model$stats$ElogV[,model$order$orderDeviance])
    } else {
        V <- as.matrix(model$V[,model$order$orderDeviance])
    }

    return(V)

}
