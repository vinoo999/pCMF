### PCA on count

# testing model 1

rm(list=ls())

source("/home/durif/source_code/countMatrixFactor/set_working_dir.R")
source("sources/projectSources.R")


######################
### Cpp version
######################

## generating the data
n = 100
p = 50
K = 10



## need of a priori values for gamma distribution
blockAlpha1 = blockMatrixGamma(nrow=n, ncol=K, nblock=5)
alpha1 = blockAlpha1$mat
alpha2 = matrix(1, nrow=n, ncol=K)
blockBeta1 = blockMatrixGamma(nrow=p, ncol=K, nblock=5)
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

lambda = rep(1, ncomp)
mu = rep(1, ncomp)

res1 = matrixFactor(data1$X, ncomp, phi01, phi02, theta01, theta02,
                    alpha01, alpha02, beta01, beta02,
                    lambda, mu,
                    iterMax=100, epsilon=1e-4, pen=TRUE, sparse=TRUE)

str(res1)
myOrder1 = res1$order$orderDeviance



res2 = matrixFactor(data1$X, ncomp, phi01, phi02, theta01, theta02,
                    alpha01, alpha02, beta01, beta02,
                    lambda, mu,
                    iterMax=100, epsilon=1e-4, pen=TRUE)

str(res2)
myOrder2 = res2$order$orderDeviance


###### comparison of the results

setwd(FIGUREDIR)

## log-likelihood
plot(res1$logLikelihood$condLogLike, xlab="iteration", ylab="conditional log likelihood", col="blue", type="l")
points(res2$logLikelihood$condLogLike, col="red")

plot(res1$logLikelihood$margLogLike, xlab="iteration", ylab="complete log likelihood", col="blue", type="l")
points(res2$logLikelihood$margLogLike, col="red")

plot(res1$logLikelihood$elbo, xlab="iteration", ylab="elbo", col="blue", type="l")
points(res2$logLikelihood$elbo, col="red")

## norm gap
plot(res1$normGap[-1], xlab="iteration", ylab="normalized gap", col="blue", type="b")
points(res2$normGap[-1], col="red")

## exp var
plot(res1$expVariance$expVar0, xlab="iteration", ylab="expVar0", col="blue", type="b")
points(res2$expVariance$expVar0, col="red")
plot(res1$expVariance$expVarU, xlab="iteration", ylab="expVarU", col="blue", type="b")
points(res2$expVariance$expVarU, col="red")
plot(res1$expVariance$expVarV, xlab="iteration", ylab="expVarV", col="blue", type="b")
points(res2$expVariance$expVarV, col="red")

## depending on K
plot(res1$criteria_k$kDeviance, xlab="k", ylab="deviance", col="blue", type="b")
points(res2$criteria_k$kDeviance, col="red")
plot(res1$criteria_k$kExpVar0, xlab="k", ylab="expVar0", col="blue", type="b")
points(res2$criteria_k$kExpVar0, col="red")
plot(res1$criteria_k$kExpVarU, xlab="k", ylab="expVarU", col="blue", type="b")
points(res2$criteria_k$kExpVarU, col="red")
plot(res1$criteria_k$kExpVarV, xlab="k", ylab="expVarV", col="blue", type="b")
points(res2$criteria_k$kExpVarV, col="red")

## order
res1$order$orderDeviance
res1$order$orderExpVar0
res1$order$orderExpVarU
res1$order$orderExpVarV


