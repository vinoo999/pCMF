
library(pCMF)

##########################
### Example of pCMF use
##########################

## data specifications
n = 100
p = 500
K = 10

## need of a priori values for gamma distribution
signalBlock = matrix(c(1,3,4,2), nrow=2, ncol=2)
blockAlpha1 = blockMatrix(nrow=n, ncol=K, nRowBlock=2, nColBlock=2, signalBlock=signalBlock)
alpha1 = blockAlpha1$mat
alpha2 = matrix(1, nrow=n, ncol=K)

signalBlock = matrix(rev(c(1,3,4,2)), nrow=2, ncol=2)
blockBeta1 = blockMatrix(nrow=p, ncol=K, nRowBlock=2, nColBlock=2, signalBlock=signalBlock)
beta1 = blockBeta1$mat
beta2 = matrix(1, nrow=p, ncol=K)

## generating the data
data1 = dataGeneration(n=n, p=p, K=K, alpha1=alpha1, alpha2=alpha2, beta1=beta1, beta2=beta2, ZI=FALSE)
str(data1)

X = data1$X

## heatmap
matrixHeatmap(alpha1, xlab="k = 1...K", ylab="i = 1...n")
matrixHeatmap(beta1, xlab="k = 1...K", ylab="j = 1...p")
matrixHeatmap(data1$X, xlab="j = 1...p", ylab="i = 1...n")
matrixHeatmap(data1$U, xlab="k = 1...K", ylab="i = 1...n")

####### TESTING ALGO

ncomp=2

res1 <- pCMF(X, ncomp, iterMax=500, iterMin=100, epsilon=1e-3, verbose=TRUE, sparse=FALSE, ZI=FALSE, ncores=16)

str(res1)

# elbo
plot(res1$logLikelihood$elbo, xlab="iteration", ylab="elbo", col="blue", type="l")

# convergence criterion
plot(res1$normGap[-1], xlab="iteration", ylab="norm. gap", col="blue", type="l")

plot(res1$criteria_k$kDeviance, type="l")

### INDIVIDUALS

# representation of individuals in log space (matrix logU)
U <- getU(res1)
plot(U, col=blockAlpha1$idRows)

# heatmap
matrixHeatmap(U)

# representation of individuals (matrix U)
U <- getU(res1, log_representation=FALSE)
plot(U, col=blockAlpha1$idRows)

# heatmap
matrixHeatmap(U)

### VARIABLES

# representation of variables in log space (matrix logV)
V <- getV(res1)
plot(V, col=blockBeta1$idRows)

# heatmap
matrixHeatmap(V)

# representation of variables (matrix V)
V <- getV(res1, log_representation=FALSE)
plot(V, col=blockBeta1$idRows)

# heatmap
matrixHeatmap(V)



# percentage of explained deviance
expDev(res1, X)
