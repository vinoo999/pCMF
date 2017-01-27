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


icl <- function(model, X) {

    n <- nrow(X)
    p <- ncol(X)

    U <- model$U
    V <- model$V
    lambda <- model$lambda

    K <- ncol(U)

    alpha1 <- model$params$alpha1
    alpha2 <- model$params$alpha2

    beta1 <- model$params$beta1
    beta2 <- model$params$beta2

    ### compute criterion
    res <- poisLoglike(X, lambda) + gammaLoglike(U, alpha1, alpha2) + gammaLoglike(V, beta1, beta2) - K * log(n*p)

    return(res)
}



poisLoglike = function(X, lambda) {
    X.ij = as.vector(X)
    lambda.ij = as.vector(lambda)

    return(sum(X.ij * log(lambda.ij) - lambda.ij - lfactorial(X.ij)))
}


gammaLoglike = function(X, shape, rate) {
    X.ij = as.vector(X)
    shape.ij = as.vector(shape)
    rate.ij = as.vector(rate)

    return(sum( (shape.ij - 1)*log(X.ij)  + shape.ij*log(rate.ij) - rate.ij * X.ij - lgamma(shape.ij) ))
}