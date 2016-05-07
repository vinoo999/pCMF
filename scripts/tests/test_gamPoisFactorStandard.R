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

res1 = matrixFactor(data1$X, ncomp, phi01, phi02, theta01, theta02, alpha01, alpha02, beta01, beta02, iterMax=100)


str(res1)


# ## INITvalues
# ## need of a priori values for gamma distribution
# alpha1=1
# beta1=1
#
# init0 = initialization1random(data1$X, n, p, alpha1=alpha1, beta1=beta1, noiseAlpha2=1, noiseBeta2=1, ncomp=ncomp)
#
# ## RUNNNING
# res1 = algo1ToC.varInf(X=data1$X, ncomp=ncomp, phi01=init0$phi01, phi02=init0$phi02,
#                        theta01=init0$theta01, theta02=init0$theta02,
#                        alpha1=init0$alpha01, alpha2=init0$alpha02, beta1=init0$beta01,
#                        beta2=init0$beta02,
#                        iterMax=100, epsilon=1e-4, order=0, rstab=5, verbose=T,
#                        criterion="expVar1")
#
# res2 = algo1.varInf(X=data1$X, ncomp=ncomp, phi01=init0$phi01, phi02=init0$phi02,
#                     theta01=init0$theta01, theta02=init0$theta02,
#                     alpha1=init0$alpha01, alpha2=init0$alpha02, beta1=init0$beta01,
#                     beta2=init0$beta02,
#                     iterMax=100, epsilon=1e-4, order=0, rstab=5, verbose=T,
#                     criterion="expVar1")
#
# str(res1)
# str(res2)
#
# ###### comparison of the results
#
# setwd(FIGUREDIR)
#
# ## variational parameters
# plot(as.vector(res1$phi1), as.vector(res2$phi1))
# abline(a=0,b=1)
#
# plot(as.vector(res1$phi2), as.vector(res2$phi2))
# abline(a=0,b=1)
#
# plot(as.vector(res1$theta1), as.vector(res2$theta1))
# abline(a=0,b=1)
#
# plot(as.vector(res1$theta2), as.vector(res2$theta2))
# abline(a=0,b=1)
#
# ## vs true values
# plot(data1$U[,c(1,3)])
# plot(res1$U[,c(1,3)])
#
# plot(data1$alpha1[,c(1,2)])
# plot(res1$phi1[,c(2,5)])
#
#
# ## log-likelihood
# plot(res1$loglike, xlab="iteration", ylab="log likelihood", col="blue", type="l")
# points(res2$loglike[1:(res2$nbIter-1)], col="red")
#
# plot(res1$loglike, res2$loglike[1:(res2$nbIter)])
# abline(a=0,b=1)
#
# plot(res1$loglikeCond, xlab="iteration", ylab="log likelihood", col="blue", type="l")
# points(res2$loglikeCond[1:(res2$nbIter-1)], col="red")
#
# plot(res1$loglikePrior, xlab="iteration", ylab="log likelihood", col="blue", type="l")
# points(res2$loglikePrior[1:(res2$nbIter-1)], col="red")
#
# plot(res1$loglikePosterior, xlab="iteration", ylab="log likelihood", col="blue", type="l")
# points(res2$loglikePosterior[1:(res2$nbIter-1)], col="red")
#
#
#
# ## noramlized gap
# plot(res1$normGap[-1], xlab="iteration", ylab="norm. gap")
# points(res2$normGap[2:(res2$nbIter)], col="red")
#
# lim=0.001
# plot(res1$normGap, res2$normGap, log="xy")
# abline(a=0,b=1)
#
# ## explained variance
# plot(res1$explainedVar1, xlab="iteration", ylab="exp var", col="blue", type="l")
# points(res1$explainedVar1, col="red", type="l")
#
#
# ## output
# layout(matrix(1:2, ncol=2))
# plot(res1$loglike, xlab="iteration", ylab = "", main="log likelihood", col="blue")
# plot(res1$normGap[-1], xlab="iteration", ylab="", main="norm. gap")
#
# ##### benchmark
# # library(microbenchmark)
# #
# # microbenchmark(
# #     algo1ToC.varInf(X=data1$X, K=K, phi0=phi0, theta0=theta0, alpha=alpha, beta=beta, iterMax=200),
# #     #algo1.varInf(X=data1$X, K=K, phi0=phi0, theta0=theta0, alpha=alpha, beta=beta, iterMax=200),
# #     times=10L
# #     )
