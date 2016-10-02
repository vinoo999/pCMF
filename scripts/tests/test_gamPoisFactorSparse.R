### PCA on count

# testing model 1

rm(list=ls())

source("/home/durif/source_code/countMatrixFactor/set_working_dir.R")
source("sources/projectSources.R")


######################
### Cpp version
######################

## generating the data
n = 20
p = 10
K = 4
nblock=2

## need of a priori values for gamma distribution
signalBlock <- matrix(4, nblock, nblock) + diag(8, nblock)
colOrder <- sample.int(n=nblock, size=nblock, replace=FALSE)
signalBlockA <- signalBlock[,colOrder]
blockAlpha1 <- blockMatrix(nrow=n, ncol=K, nRowBlock=nblock, nColBlock=nblock, signalBlock=signalBlockA)
alpha1 <- blockAlpha1$mat
alpha2 <- matrix(1, nrow=n, ncol=K)

signalBlock <- matrix(4, nblock, nblock) + diag(8, nblock)
colOrder <- sample.int(n=nblock, size=nblock, replace=FALSE)
signalBlockB <- signalBlock[,colOrder]
blockBeta1 <- blockMatrix(nrow=p, ncol=K, nRowBlock=nblock, nColBlock=nblock, signalBlock=signalBlockB)
beta1 <- blockBeta1$mat
beta2 <- matrix(1, nrow=p, ncol=K)

## generating the data
data1 = dataGeneration(n=n, p=p, K=K, alpha1=alpha1, alpha2=alpha2, beta1=beta1, beta2=beta2)
str(data1)

X = cbind(data1$X, matrix(rpois(n*10, lambda=2), nrow=n))
n=n
p=p+10

# X=data1$X

## heatmap
# matrixHeatmap(alpha1, xlab="k = 1...K", ylab="i = 1...n")
# matrixHeatmap(beta1, xlab="k = 1...K", ylab="j = 1...p")
# matrixHeatmap(data1$X, xlab="j = 1...p", ylab="i = 1...n")
# matrixHeatmap(data1$U, xlab="k = 1...K", ylab="i = 1...n")
# matrixHeatmap(data1$V, xlab="k = 1...K", ylab="i = 1...n")

####### TESTING ALGO
ncomp=4

alpha01 = matrix(1, nrow=n, ncol=ncomp)
alpha02 = matrix(1, nrow=n, ncol=ncomp)

beta01 = matrix(1, nrow=p, ncol=ncomp)
beta02 = matrix(1, nrow=p, ncol=ncomp)

phi01 = matrix(1, nrow=n, ncol=ncomp)
phi02 = matrix(1, nrow=n, ncol=ncomp)

theta01 = matrix(1, nrow=p, ncol=ncomp)
theta02 = matrix(1, nrow=p, ncol=ncomp)

res1 = matrixFactor(X, ncomp, phi01, phi02, theta01, theta02, alpha01, alpha02, beta01, beta02,
                    iterMax=200, iterMin=100, epsilon=1e-4, algo="EM", verbose=TRUE, sparse=TRUE)

res2 = matrixFactor(X, ncomp, phi01, phi02, theta01, theta02, alpha01, alpha02, beta01, beta02,
                    iterMax=200, iterMin=100, epsilon=1e-4, algo="EM", verbose=TRUE, sparse=FALSE)

str(res1)

print(res1$criteria_k$kDeviance)
myOrder = res1$order$orderDeviance
print(myOrder)

res1$sparseParams$probSparse
res1$sparseParams$probSparsePrior
res1$sparseParams$sparseIndic

layout(matrix(1:2,ncol=2))
plot(res1$logLikelihood$elbo[-(1:2)], xlab="iteration", ylab="elbo", col="blue", type="l")
plot(res1$normGap[-1], xlab="iteration", ylab="normalized gap", col="blue", type="b", log="y")

matrixHeatmap(res1$sparseParams$probSparse)
matrixHeatmap(res1$sparseParams$sparseIndic)

matrixHeatmap(res1$V)
matrixHeatmap(res2$V)

matrixHeatmap(res1$U)
matrixHeatmap(res2$U)

matrixHeatmap(res1$params$theta1)
matrixHeatmap(res2$params$theta1)

matrixHeatmap(res1$params$theta2)
matrixHeatmap(res2$params$theta2)

matrixHeatmap(res1$params$beta1)
matrixHeatmap(res2$params$beta1)

matrixHeatmap(res1$params$beta2)
matrixHeatmap(res2$params$beta2)

matrixHeatmap(res1$params$phi1)
matrixHeatmap(res2$params$phi1)

matrixHeatmap(res1$params$phi2)
matrixHeatmap(res2$params$phi2)

matrixHeatmap(res1$params$alpha1)
matrixHeatmap(res2$params$alpha1)

matrixHeatmap(res1$params$alpha2)
matrixHeatmap(res2$params$alpha2)



# res1$V[,3]
#
# ###### comparison of the results
#
# setwd(FIGUREDIR)
#
# ## log-likelihood
# plot(res1$logLikelihood$condLogLike, xlab="iteration", ylab="conditional log likelihood", col="blue", type="l")
#
# plot(res1$logLikelihood$margLogLike, xlab="iteration", ylab="complete log likelihood", col="blue", type="l")
#
# plot(res1$logLikelihood$elbo[-(1:2)], xlab="iteration", ylab="elbo", col="blue", type="l")
#
# ## norm gap
# plot(res1$normGap[-1], xlab="iteration", ylab="normalized gap", col="blue", type="b", log="y")
#
# ## exp var
# plot(res1$expVariance$expVar0, xlab="iteration", ylab="expVar0", col="blue", type="b")
# plot(res1$expVariance$expVarU, xlab="iteration", ylab="expVarU", col="blue", type="b")
# plot(res1$expVariance$expVarV, xlab="iteration", ylab="expVarV", col="blue", type="b")
#
# ## depending on K
# plot(res1$criteria_k$kDeviance, xlab="k", ylab="deviance", col="blue", type="b")
# plot(res1$criteria_k$kExpVar0, xlab="k", ylab="expVar0", col="blue", type="b")
# plot(res1$criteria_k$kExpVarU, xlab="k", ylab="expVarU", col="blue", type="b")
# plot(res1$criteria_k$kExpVarV, xlab="k", ylab="expVarV", col="blue", type="b")
#
# ## order
# res1$order$orderDeviance
# res1$order$orderExpVar0
# res1$order$orderExpVarU
# res1$order$orderExpVarV
#
#
