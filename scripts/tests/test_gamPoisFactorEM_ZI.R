
rm(list=ls())

source("/home/durif/source_code/countMatrixFactor/set_working_dir.R")
source("sources/projectSources.R")

# library(fields)

######################
### Cpp version
######################

## generating the data
n = 100
p = 50
K = 10



## need of a priori values for gamma distribution
signalBlock = matrix(c(1,2,4,2), nrow=2, ncol=2)
blockAlpha1 = blockMatrix(nrow=n, ncol=K, nRowBlock=2, nColBlock=2, signalBlock=signalBlock)
alpha1 = blockAlpha1$mat
alpha2 = matrix(1, nrow=n, ncol=K)

signalBlock = matrix(rev(c(1,3,4,2)), nrow=2, ncol=2)
blockBeta1 = blockMatrix(nrow=p, ncol=K, nRowBlock=2, nColBlock=2, signalBlock=signalBlock)
beta1 = blockBeta1$mat
beta2 = matrix(1, nrow=p, ncol=K)
prob1 = round(runif(p, 0.3, 0.7), digits=2)

## generating the data
data1 = dataGeneration(n=n, p=p, K=K, alpha1=alpha1, alpha2=alpha2, beta1=beta1, beta2=beta2, ZI=TRUE, prob1=prob1)
str(data1)

cbind(apply(data1$X,2, function(x) sum(x!=0)), prob1)

## heatmap
# matrixHeatmap(alpha1, xlab="k = 1...K", ylab="i = 1...n")
# matrixHeatmap(beta1, xlab="k = 1...K", ylab="j = 1...p")
# matrixHeatmap(data1$X, xlab="j = 1...p", ylab="i = 1...n")
# matrixHeatmap(data1$U, xlab="k = 1...K", ylab="i = 1...n")

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

res1 = matrixFactor(data1$X, ncomp,
                    phi01, phi02, theta01, theta02,
                    alpha01, alpha02, beta01, beta02,
                    iterMax=500, epsilon=1e-5,
                    ZI=TRUE, algo = "EM")

str(res1)

# myOrder = res1$order$orderExpVarU

cbind(apply(data1$X,2,function(x) sum(x!=0)), prob1, res1$ZIparams$probPrior)


# elbo
plot(res1$logLikelihood$elbo[-(1:10)], xlab="iteration", ylab="elbo", col="blue", type="l")
plot(res1$normGap[-(1:10)], xlab="iteration", ylab="norm. gap", col="blue", type="l")

# graph
U = res1$U[, res1$order$orderDeviance]
plot(U, col=blockAlpha1$idRows)

# ZI indic
matrixHeatmap(data1$ZIind)
matrixHeatmap(res1$ZIparams$prob)
matrixHeatmap(data1$Xnzi>0)

## higher k
ncomp=6

alpha01 = matrix(1, nrow=n, ncol=ncomp)
alpha02 = matrix(1, nrow=n, ncol=ncomp)

beta01 = matrix(1, nrow=p, ncol=ncomp)
beta02 = matrix(1, nrow=p, ncol=ncomp)

phi01 = matrix(1, nrow=n, ncol=ncomp)
phi02 = matrix(1, nrow=n, ncol=ncomp)

theta01 = matrix(1, nrow=p, ncol=ncomp)
theta02 = matrix(1, nrow=p, ncol=ncomp)

res1 = matrixFactor(data1$Xnzi, ncomp,
                    phi01, phi02, theta01, theta02,
                    alpha01, alpha02, beta01, beta02,
                    iterMax=500, epsilon=1e-5,
                    ZI=TRUE, algo = "EM")

res2 = matrixFactor(data1$Xnzi, ncomp,
                    phi01, phi02, theta01, theta02,
                    alpha01, alpha02, beta01, beta02,
                    iterMax=500, epsilon=1e-5,
                    ZI=FALSE, algo = "EM")

# U
matrixHeatmap(res1$U)
matrixHeatmap(res2$U)
# V
matrixHeatmap(res1$V)
matrixHeatmap(res2$V)

# graph
U = res2$U[, res2$order$orderDeviance]
plot(U, col=blockAlpha1$idRows)

### ordering
orderInd <- sample.int(n,n)
orderVar <- sample.int(p,p)
X = data1$Xnzi[orderInd,]
X = X[,orderVar]
res2 = matrixFactor(X, ncomp,
                    phi01, phi02, theta01, theta02,
                    alpha01, alpha02, beta01, beta02,
                    iterMax=500, epsilon=1e-5,
                    ZI=FALSE, algo = "EM")

# graph
U = res2$U[, res2$order$orderDeviance]
plot(U, col=blockAlpha1$idRows[orderInd])
