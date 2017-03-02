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

    K0 <- ncol(U)

    # alpha1 <- model$params$alpha1
    # alpha2 <- model$params$alpha2
    #
    # beta1 <- model$params$beta1
    # beta2 <- model$params$beta2

    alpha1 <- model$params$phi1
    alpha2 <- model$params$phi2

    beta1 <- model$params$theta1
    beta2 <- model$params$theta2

    nu_K <- K0 * (n+p)

    cst <- log(K0)
    # cst <- log(K0)

    pen1 <- (nu_K/2) * (cst)
    pen2 <- K0 * log(n) + jeffreys(alpha1[1,],alpha2[1,]) + jeffreys(beta1[1,],beta2[1,]) - (nu_K/2) * (cst)


    ### compute criterion
    res1 <- poisLoglike(X, lambda) + gammaLoglike(U, alpha1, alpha2) + gammaLoglike(V, beta1, beta2)
    res2 <- res1 - pen1
    res3 <- res1 + pen2

    return(c(res1, res2, res3))
}


bic <- function(model, X) {

    n <- nrow(X)
    p <- ncol(X)

    U <- model$U
    V <- model$V
    lambda <- model$lambda

    K0 <- ncol(U)

    nu_K <- K0 * (n+p)

    cst <- log(n)

    ### compute criterion
    res1 <- poisLoglike(X, lambda)
    res2 <- res1 - (nu_K/2) * (cst) # + 2*log(2) - log(2*pi))

    return(c(res1, res2))
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

jeffreys <- function(shape,rate) {
    return(sum(0.5*log(sqrt(trigamma(shape)*shape -1)) - log(rate)))
}



icl2 <- function(model, X) {

    n <- nrow(X)
    p <- ncol(X)

    myOrder <- model$order$orderDeviance

    U <- model$U[,myOrder]
    V <- model$V[,myOrder]

    K <- ncol(U)

    # alpha1 <- model$params$alpha1
    # alpha2 <- model$params$alpha2
    #
    # beta1 <- model$params$beta1
    # beta2 <- model$params$beta2

    alpha1 <- model$params$phi1
    alpha2 <- model$params$phi2

    beta1 <- model$params$theta1
    beta2 <- model$params$theta2

    criterion <- sapply(1:K, function(k) {

        nu_K <- k * (n+p)

        cst <- log(n*p)
        # cst <- log(n)

        Utmp <- as.matrix(U[,1:k])
        Vtmp <- as.matrix(V[,1:k])
        alpha1tmp <- as.matrix(alpha1[,1:k])
        alpha2tmp <- as.matrix(alpha2[,1:k])
        beta1tmp <- as.matrix(beta1[,1:k])
        beta2tmp <- as.matrix(beta2[,1:k])

        lambda <- Utmp %*% t(Vtmp)

        ### compute criterion
        res1 <- poisLoglike(X, lambda) + gammaLoglike(Utmp, alpha1tmp, alpha2tmp) + gammaLoglike(Vtmp, beta1tmp, beta2tmp) - (nu_K/2) * (cst + 2*log(2) + log(2*pi))
        res2 <- poisLoglike(X, lambda) + gammaLoglike(Utmp, alpha1tmp, alpha2tmp) + gammaLoglike(Vtmp, beta1tmp, beta2tmp)

        return(c(res1,res2))
    })

    return(t(criterion))
}


bic2 <- function(model, X) {

    n <- nrow(X)
    p <- ncol(X)

    myOrder <- model$order$orderDeviance

    U <- model$U[,myOrder]
    V <- model$V[,myOrder]

    K <- ncol(U)

    # alpha1 <- model$params$alpha1
    # alpha2 <- model$params$alpha2
    #
    # beta1 <- model$params$beta1
    # beta2 <- model$params$beta2

    alpha1 <- model$params$phi1
    alpha2 <- model$params$phi2

    beta1 <- model$params$theta1
    beta2 <- model$params$theta2

    criterion <- sapply(1:K, function(k) {

        nu_K <- k * (n+p)

        cst <- log(n*p)
        # cst <- log(n)

        Utmp <- as.matrix(U[,1:k])
        Vtmp <- as.matrix(V[,1:k])
        alpha1tmp <- as.matrix(alpha1[,1:k])
        alpha2tmp <- as.matrix(alpha2[,1:k])
        beta1tmp <- as.matrix(beta1[,1:k])
        beta2tmp <- as.matrix(beta2[,1:k])

        lambda <- Utmp %*% t(Vtmp)

        ### compute criterion
        res1 <- poisLoglike(X, lambda) - (nu_K/2) * (cst + 2*log(2) + log(2*pi))

        return(res1)
    })

    return(criterion)
}






