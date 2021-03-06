

######################
### Cpp version
######################

rm(list=ls())

RDIR <- system("git rev-parse --show-toplevel", intern=TRUE)
source(paste0(RDIR, "/set_working_dir.R"))
library(pCMF, lib.loc=paste0(WORKINGDIR, "/", myLib))

## generating the data
n <- 100
p <- 100
K <- 10



## need of a priori values for gamma distribution
signalBlock <- matrix(c(1,3,4,2), nrow=2, ncol=2)
blockAlpha1 <- blockMatrix(nrow=n, ncol=K, nRowBlock=2, nColBlock=2, signalBlock=signalBlock)
alpha1 <- blockAlpha1$mat
alpha2 <- matrix(1, nrow=n, ncol=K)

signalBlock <- matrix(rev(c(1,3,4,2)), nrow=2, ncol=2)
blockBeta1 <- blockMatrix(nrow=p, ncol=K, nRowBlock=2, nColBlock=2, signalBlock=signalBlock)
beta1 <- blockBeta1$mat
beta2 <- matrix(1, nrow=p, ncol=K)

## generating the data
data1 <- dataGeneration(n=n, p=p, K=K, alpha1=alpha1, alpha2=alpha2, beta1=beta1, beta2=beta2, ZI=FALSE)
str(data1)

X <- data1$X

## heatmap
matrixHeatmap(alpha1, xlab="k = 1...K", ylab="i = 1...n")
matrixHeatmap(beta1, xlab="k = 1...K", ylab="j = 1...p")
matrixHeatmap(data1$X, xlab="j = 1...p", ylab="i = 1...n")
matrixHeatmap(data1$U, xlab="k = 1...K", ylab="i = 1...n")

####### TESTING ALGO

ncomp <- 2

res1 <- pCMF(X, ncomp, iterMax=500, iterMin=100, epsilon=1e-3, verbose=TRUE, ncores=16,
             nbInit=10, iterMaxInit=30, noise=0.5, seed=NULL)

str(res1)

# elbo
plot(res1$logLikelihood$elbo, xlab="iteration", ylab="elbo", col="blue", type="l")

# convergence criterion
plot(res1$normGap[-1], xlab="iteration", ylab="norm. gap", col="blue", type="l")

## deviance
plot(res1$criteria_k$kDeviance, type="l", ylab="deviance")

## percentage of explained deviance
expDev(res1, X)

# individuals
U <- getU(res1, log_representation=TRUE)
str(U)
graphU(res1, axes=1:2, labels=factor(blockAlpha1$idRows),
       log_representation=TRUE, edit_theme=TRUE, graph=TRUE)

matrixHeatmap(U)


## test multi init
layout(matrix(1:10, nrow=2))
for(i in 1:10) {
    res1 <- pCMF(X, ncomp, iterMax=500, iterMin=100, epsilon=1e-3, verbose=FALSE, ncores=16,
                 nbInit=10, iterMaxInit=30, noise=0.5, seed=NULL)
    U <- getU(res1, log_representation=FALSE)
    plot(U, col=blockAlpha1$idRows)

}
