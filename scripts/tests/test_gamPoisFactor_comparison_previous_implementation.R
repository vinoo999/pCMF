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






source("~/source_code/countMatrixFactorArchives/sources/light.sources.R")

init0 = initialize1(n, p, alpha1=alpha01, alpha2=alpha02, beta1=beta01, beta2=beta02, ncomp=ncomp, noisePos=0, noiseLevel=0.2)

res2 = algo1.varInf(X=data1$X, ncomp=ncomp, phi01=init0$phi01, phi02=init0$phi02,
                                        theta01=init0$theta01, theta02=init0$theta02,
                                        alpha1=init0$alpha01, alpha2=init0$alpha02, beta1=init0$beta01,
                                        beta2=init0$beta02,
                                        iterMax=100, epsilon=1e-4, order=0, rstab=5, verbose=T,
                                        criterion="expVar1")

res1 = matrixFactor(data1$X, ncomp, init0$phi01, init0$phi02, init0$theta01, init0$theta02, init0$alpha01, init0$alpha02, init0$beta01, init0$beta02, iterMax=100, epsilon=1e-4)

str(res1)

myOrder = res1$order$orderExpVarU

###### comparison of the results

setwd(FIGUREDIR)

## variational parameters
plot(as.vector(res1$params$phi1[,myOrder]), as.vector(res2$phi1))
abline(a=0,b=1)

plot(as.vector(res1$params$phi2[,myOrder]), as.vector(res2$phi2))
abline(a=0,b=1)

plot(as.vector(res1$params$theta1[,myOrder]), as.vector(res2$theta1))
abline(a=0,b=1)

plot(as.vector(res1$params$theta2[,myOrder]), as.vector(res2$theta2))
abline(a=0,b=1)

## log-likelihood
plot(res1$logLikelihood$condLogLike, xlab="iteration", ylab="log likelihood", col="blue", type="l")
points(res2$loglikeCond[1:(res2$nbIter-1)], col="red")

plot(res1$logLikelihood$margLogLike, xlab="iteration", ylab="log likelihood", col="blue", type="l")
points(res2$loglike[1:(res2$nbIter-1)], col="red")
