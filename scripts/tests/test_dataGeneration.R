### PCA on count

# testing model 1

rm(list=ls())

source("/home/durif/source_code/countMatrixFactor/set_working_dir.R")
source("sources/projectSources.R")

library(fields)

######################
### Cpp version
######################

## generating the data
n = 100
p = 50
K = 10


## need of a priori values for gamma distribution
signalBlock = matrix(1:(10*5), nrow=10, ncol=5)
blockAlpha1 = blockMatrix(nrow=n, ncol=K, nRowBlock=10, nColBlock=5, signalBlock=signalBlock)
alpha1 = blockAlpha1$mat

image.plot(alpha1, xaxt="n", yaxt="n", xlab="i", ylab="k")

alpha2 = matrix(1, nrow=n, ncol=K)
blockBeta1 = blockMatrixGamma(nrow=p, ncol=K, nRowBlock=5, nColBlock=5)
beta1 = blockBeta1$mat
beta2 = matrix(1, nrow=p, ncol=K)


## generating the data
data1 = dataGeneration(n=n, p=p, K=K, alpha1=alpha1, alpha2=alpha2, beta1=beta1, beta2=beta2)
str(data1)