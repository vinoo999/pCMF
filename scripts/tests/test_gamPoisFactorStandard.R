### PCA on count

# testing model 1

rm(list=ls())

source("/home/durif/source_code/countMatrixFactor/set_working_dir.R")
source("sources/projectSources.R")


######################
### Cpp version
######################

## generating the data
n = 50
p = 20
K = 10



## need of a priori values for gamma distribution
blockAlpha1 = blockMatrixGamma(nrow=n, ncol=K, nblock=5)
alpha1 = blockAlpha1$mat
alpha2 = matrix(1, nrow=n, ncol=K)
blockBeta1 = blockMatrixGamma(nrow=n, ncol=K, nblock=5)
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

res1 = matrixFactor(data1$X, ncomp, phi01, phi02, theta01, theta02, alpha01, alpha02, beta01, beta02, iterMax=100, epsilon=1e-4)

str(res1)

myOrder = res1$order$orderExpVarU

###### comparison of the results

setwd(FIGUREDIR)

## log-likelihood
plot(res1$logLikelihood$condLogLike, xlab="iteration", ylab="conditional log likelihood", col="blue", type="l")

plot(res1$logLikelihood$margLogLike, xlab="iteration", ylab="complete log likelihood", col="blue", type="l")

plot(res1$logLikelihood$elbo, xlab="iteration", ylab="elbo", col="blue", type="l")

## norm gap
plot(res1$normGap[-1], xlab="iteration", ylab="normalized gap", col="blue", type="b")

## exp var
plot(res1$expVariance$expVar0, xlab="iteration", ylab="expVar0", col="blue", type="b")
plot(res1$expVariance$expVarU, xlab="iteration", ylab="expVarU", col="blue", type="b")
plot(res1$expVariance$expVarV, xlab="iteration", ylab="expVarV", col="blue", type="b")

## depending on K
plot(res1$order$kDeviance, xlab="k", ylab="deviance", col="blue", type="b")
plot(res1$order$kExpVar0, xlab="k", ylab="expVar0", col="blue", type="b")
plot(res1$order$kExpVarU, xlab="k", ylab="expVarU", col="blue", type="b")
plot(res1$order$kExpVarV, xlab="k", ylab="expVarV", col="blue", type="b")

## order
res1$order$orderDeviance
res1$order$orderExpVar0
res1$order$orderExpVarU
res1$order$orderExpVarV


