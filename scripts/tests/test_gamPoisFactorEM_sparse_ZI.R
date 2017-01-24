
rm(list=ls())

source("/home/durif/source_code/countMatrixFactor/set_working_dir.R")
source("sources/projectSources.R")

# library(fields)

######################
### Cpp version
######################


shift <- function(x) {return(c(tail(x, length(x)-1), head(x,1)))}
compos <- function(x, n, fun) {
    res <- matrix(NA, nrow=length(x), ncol=n)
    res[,1] = x
    if(n>=2) {
        for(i in 2:n) {
            res[,i] = fun(res[,i-1])
        }
    }
    return(res)
}

## generating the data
n = 100
p = 100
K = 10

nblockU = 2
nblockV = 2
epsilonU = 4
espilonV = 4
p0 = 50



## organization of U
signalBlockU <- matrix(rep(round(runif(K, 2, 6), digits=2), each=nblockU), nrow=nblockU, ncol=K)
epsilonBlockU <- compos(seq(0,1, length.out=nblockU), K, shift)
signalBlockU <- signalBlockU + epsilonU * epsilonBlockU
blockAlpha1 <- blockMatrix(nrow=n, ncol=K, nRowBlock=nblockU, nColBlock=K, signalBlock=signalBlockU)
alpha1 <- blockAlpha1$mat
alpha2 <- matrix(1, nrow=n, ncol=K)

## organization of V
signalBlockV <- matrix(rep(round(runif(K, 2, 6), digits=2), each=nblockV), nrow=nblockV, ncol=K)
epsilonBlockV <- compos(seq(0,1, length.out=nblockV), K, shift)
signalBlockV <- signalBlockV + epsilonU * epsilonBlockV
blockBeta1 <- blockMatrix(nrow=p0, ncol=K, nRowBlock=nblockV, nColBlock=K, signalBlock=signalBlockV)
beta1 <- NULL
if(p-p0 > 0) {
    beta1 <- rbind(blockBeta1$mat, matrix(0.7, nrow=p-p0, ncol=K))
} else {
    beta1 <- blockBeta1$mat
}
beta2 <- matrix(1, nrow=p, ncol=K)

## ZI
prob1 <- round(runif(p, 0.5, 0.8), digits=2)

## generating the data
data1 <- dataGeneration(n=n, p=p, K=K, alpha1=alpha1, alpha2=alpha2, beta1=beta1, beta2=beta2,
                        ZI=TRUE, prob1=prob1, reorder=TRUE)
str(data1)

## heatmap
matrixHeatmap(alpha1, xlab="k = 1...K", ylab="i = 1...n")
matrixHeatmap(beta1, xlab="k = 1...K", ylab="j = 1...p")
matrixHeatmap(data1$X, xlab="j = 1...p", ylab="i = 1...n")
matrixHeatmap(data1$U, xlab="k = 1...K", ylab="i = 1...n")
matrixHeatmap(data1$V, xlab="k = 1...K", ylab="j = 1...p")

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


res1 = matrixFactor(data1$X, ncomp, phi01, phi02, theta01, theta02, alpha01, alpha02, beta01, beta02,
                    iterMax=800, iterMin=100, epsilon=1e-4, algo="EM", verbose=TRUE, ZI=TRUE, sparse=TRUE)

res2 = matrixFactor(data1$X, ncomp, phi01, phi02, theta01, theta02, alpha01, alpha02, beta01, beta02,
                    iterMax=500, iterMin=100, epsilon=1e-4, algo="EM", verbose=TRUE, ZI=TRUE, sparse=FALSE)

str(res1)

print(res1$criteria_k$kDeviance)
myOrder = res1$order$orderDeviance
print(myOrder)

# elbo
plot(res1$logLikelihood$elbo[-(1:10)], xlab="iteration", ylab="elbo", col="blue", type="l")
plot(res1$normGap[-(1:10)], xlab="iteration", ylab="norm. gap", col="blue", type="l")

# graph
U1 = res1$U[, res1$order$orderDeviance]
U2 = res2$U[, res2$order$orderDeviance]
plot(U1, col=blockAlpha1$idRows[data1$orderInd])
plot(U2, col=blockAlpha1$idRows[data1$orderInd])

# selection
print(res2$params$theta1)

matrixHeatmap(res1$sparseParams$sparseIndic[,myOrder])
matrixHeatmap(data1$beta1[data1$orderVar,])
matrixHeatmap((res1$V * res1$sparseParams$sparseIndic)[,myOrder])
matrixHeatmap(res2$V[,myOrder])
