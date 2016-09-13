### Copyright 2016-05 Ghislain DURIF
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
#' Ghislain Durif, \email{ghislain.durif@univ-lyon1.fr}
#'
#' @param n number of observations (nb of rows)
#' @param p number of variables (nb of columns)
#' @param K number of latent factors (dimension of latent subspace)
#' @param alpha1 matrix n x K, first parameters for the prior distribution on U
#' @param alpĥa2 matrix n x K, second parameters for the prior distribution on U
#' @param beta1 matrix n x K, first parameters for the prior distribution on V
#' @param beta2 matrix n x K, second parameters for the prior distribution on U
#' @param ZI boolean, indicating if the data are zero-inflated (default is FALSE)
#' @param prob0 vector of Bernoulli probability for zero-inflation (default is NULL), length p, one probability per variable
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
#' \item{prob0}{vector of Bernoulli probability for zero-inflation}
#'
#' @export
dataGeneration = function(n, p, K, alpha1, alpha2, beta1, beta2, ZI=FALSE, prob0=NULL) {

    if(ZI && is.null(prob0)) {
        stop("message from dataGeneration: zero-inflated model is asked but prob0 are not set in input")
    }

    ## generating the components
    U = sapply(1:K, function(k) rgamma(n, shape=alpha1[,k], rate=alpha2[,k])) # matrix n x K
    V = sapply(1:K, function(k) rgamma(p, shape=beta1[,k], rate=beta2[,k])) # matrix p x K

    ## generating the count
    Xnzi = matrix(rpois(n=n*p, lambda=as.vector(U %*% t(V))), nrow=n, ncol=p)

    ## generating the Bernoulli variables if necessary
    if(ZI) {
        Y = sapply(prob0, function(p) return(rbinom(n=n,size=1,prob=p)))
    } else {
        Y = matrix(1, nrow=n, ncol=p)
    }

    ## counts model (with zero-inflation if so)
    X = Xnzi * Y

    ## return
    return(list(X=X, U=U, V=V, n=n, p=p, K=K,
                alpha1=alpha1, alpha2=alpha2, beta1=beta1, beta2=beta2,
                ZI=ZI, prob0=prob0))
}
