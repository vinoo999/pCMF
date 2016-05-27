### PCA on count

# testing model 1

rm(list=ls())

source("/home/durif/source_code/countMatrixFactor/set_working_dir.R")
source("sources/projectSources.R")

library(fields)
library(ggplot2)


######################
### Cpp version
######################

## generating the data
n = 100
p = 50
K = 10



## need of a priori values for gamma distribution
signalBlock = matrix(c(1,3,4,2), nrow=2, ncol=2)
blockAlpha1 = blockMatrix(nrow=n, ncol=K, nRowBlock=2, nColBlock=2, signalBlock=signalBlock)
alpha1 = blockAlpha1$mat
image.plot(alpha1, xaxt="n", yaxt="n", xlab="i", ylab="k")
alpha2 = matrix(1, nrow=n, ncol=K)

blockBeta1 = blockMatrixGamma(nrow=p, ncol=K, nRowBlock=2, nColBlock=2, signalBlock=signalBlock)
beta1 = blockBeta1$mat
beta2 = matrix(1, nrow=p, ncol=K)


## generating the data
data1 = dataGeneration(n=n, p=p, K=K, alpha1=alpha1, alpha2=alpha2, beta1=beta1, beta2=beta2)
str(data1)



####### TESTING ALGO

ncomp=K

alpha01 = matrix(1, nrow=n, ncol=ncomp)
alpha02 = matrix(1, nrow=n, ncol=ncomp)

beta01 = matrix(1, nrow=p, ncol=ncomp)
beta02 = matrix(1, nrow=p, ncol=ncomp)

phi01 = matrix(1, nrow=n, ncol=ncomp)
phi02 = matrix(1, nrow=n, ncol=ncomp)

theta01 = matrix(1, nrow=p, ncol=ncomp)
theta02 = matrix(1, nrow=p, ncol=ncomp)

lambda = rep(0.01, ncomp)
mu = rep(0.01, ncomp)

res1 = matrixFactor(data1$X, ncomp, phi01, phi02, theta01, theta02,
                    alpha01, alpha02, beta01, beta02,
                    lambda, mu,
                    iterMax=300, epsilon=1e-4, pen=TRUE)

str(res1)
myOrder1 = res1$order$orderExpVarU

U = res1$U[, myOrder1]
plot(U, col=blockAlpha1$idRows)

V = res1$V[, myOrder1]
plot(V, col=blockBeta1$idRows)

library(fields)
image.plot(U, xaxt="n", yaxt="n", xlab="i", ylab="k")

image.plot(data1$U, xaxt="n", yaxt="n", xlab="i", ylab="k")

image.plot(V, xaxt="n", yaxt="n", xlab="i", ylab="k")

image.plot(data1$V, xaxt="n", yaxt="n", xlab="i", ylab="k")

res3 = prcomp(data1$X)
str(res3)
U = data1$X %*% res3$rotation[,1:10]
plot(U, col=blockAlpha1$idRows)
image.plot(U, xaxt="n", yaxt="n", xlab="i", ylab="k")


### ELBO
plot(res1$logLikelihood$elbo)
