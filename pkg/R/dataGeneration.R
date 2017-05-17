### Copyright 2016-05 Ghislain DURIF
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

#' @title dataGeneration
#'
#' @description
#' Simulates a data matrix X according to the (possibly
#' zero-inflated) Gamma Poisson Factor Model
#'
#' @details
#' ## generative process (without zero-inflation)
#'   X | U,V ~ Poisson(U t(V))
#'   U ~ Gamma(alpha)
#'   V ~ Gamma(beta)
#'
#' ## Generative process (with zero-inflation)
#'  X_{ij} = sum_k Z_{ijk}
#'  X_{ij} | U,V ~ Y_{ij} * Poisson(U_i t(V_j))
#'  Z_{ijk} | U,V ~ Y_{ij} * Poisson(U_{ik} V_{jk})
#'  Y_{ij} ~ Bernoulli(p_j)
#'  U ~ Gamma(alpha)
#'  V ~ Gamma(beta)
#'
#' @author
#' Ghislain Durif, \email{gd.dev@libertymail.net}
#'
#' @param n number of observations (nb of rows)
#' @param p number of variables (nb of columns)
#' @param K number of latent factors (dimension of latent subspace)
#' @param alpha1 matrix n x K, first parameters for the prior distribution on U
#' @param alpĥa2 matrix n x K, second parameters for the prior distribution on U
#' @param beta1 matrix n x K, first parameters for the prior distribution on V
#' @param beta2 matrix n x K, second parameters for the prior distribution on U
#' @param ZI boolean, indicating if the data are zero-inflated (default is FALSE)
#' @param prob1 vector of Bernoulli probability for zero-inflation (default is NULL),
#' length p, one probability per variable, probability to not get a drop-out event
#' @param rate0 vector of Bernoulli probability for zero-inflation (default is NULL),
#' enforce the ZI probabibility to depend on the count mean, the higher the smaller probability
#'
#' @return list containing the following
#' \item{X}{data matrix n x p of counts}
#' \item{U}{matrix n x K of factor coordinates (in observation space)}
#' \item{V}{matrix p x K of factor loadings (in variable space)}
#' \item{n}{number of observations (or rows in X)}
#' \item{p}{number of variables (or columns in X)}
#' \item{K}{number of latent factors (dimension of latent subspace)}
#' \item{alpha1}{matrix n x K, first parameters for the prior distribution on U}
#' \item{alpĥa2}{matrix n x K, second parameters for the prior distribution on U}
#' \item{beta1}{matrix n x K, first parameters for the prior distribution on V}
#' \item{beta2}{matrix n x K, second parameters for the prior distribution on U}
#' \item{ZI}{boolean, indicating if the data are zero-inflated (default is FALSE)}
#' \item{prob1}{vector of Bernoulli probability for zero-inflation, probability
#' to not get a drop-out event (default is NULL)}
#' \item{rate0}{vector of Bernoulli probability for zero-inflation, enforce the ZI probabibility
#' to depend on the count mean, the higher the smaller probability (default is NULL)}
#'
#' @export
dataGeneration <- function(n, p, K, alpha1, alpha2, beta1, beta2, ZI=FALSE, prob1=NULL, rate0=NULL, reorder=FALSE) {

    if(ZI & (is.null(prob1) & is.null(rate0))) {
        stop("message from dataGeneration: zero-inflated model is asked but prob1 are not set in input")
    }

    ## generating the components
    U <- sapply(1:K, function(k) rgamma(n, shape=alpha1[,k]/sqrt(K), rate=alpha2[,k])) # matrix n x K
    V <- sapply(1:K, function(k) rgamma(p, shape=beta1[,k]/sqrt(K), rate=beta2[,k])) # matrix p x K

    ## reordering individuals and genes
    orderInd <- NULL
    orderVar <- NULL
    if(reorder) {
        orderInd <- sample.int(n,n)
        orderVar <- sample.int(p,p)
    } else {
        orderInd <- 1:n
        orderVar <- 1:p
    }
    U <- U[orderInd,]
    V <- V[orderVar,]

    ## generating the count
    Xnzi <- matrix(rpois(n=n*p, lambda=as.vector(U %*% t(V))), nrow=n, ncol=p)

    ## generating the Bernoulli variables if necessary
    if(ZI) {
        if(!is.null(rate0)) {
            meanX <- apply(Xnzi, 2, mean)
            prob1 <- 1 - exp(-rate0*meanX)
            Y <- sapply(prob1, function(pi) return(rbinom(n=n,size=1,prob=pi)))
        } else {
            Y <- sapply(prob1, function(pi) return(rbinom(n=n,size=1,prob=pi)))
        }

    } else {
        Y <- matrix(1, nrow=n, ncol=p)
    }

    ## counts model (with zero-inflation if so)
    X <- Xnzi * Y

    if(!ZI) {
        Xnzi <- NULL
        Y <- NULL
    }

    ## return
    return(list(X=X, U=U, V=V, n=n, p=p, K=K,
                alpha1=alpha1, alpha2=alpha2, beta1=beta1, beta2=beta2,
                ZI=ZI, prob1=prob1, rate0=rate0, Xnzi=Xnzi, ZIind=Y,
                orderInd=orderInd, orderVar=orderVar))
}
