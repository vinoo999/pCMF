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
p = 500
K = 5

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
data1 = dataGeneration(n=n, p=p, K=K, alpha1=alpha1, alpha2=alpha2, beta1=beta1, beta2=beta2)
str(data1)

## heatmap
# matrixHeatmap(alpha1, xlab="k = 1...K", ylab="i = 1...n")
# matrixHeatmap(beta1, xlab="k = 1...K", ylab="j = 1...p")
# matrixHeatmap(data1$X, xlab="j = 1...p", ylab="i = 1...n")
# matrixHeatmap(data1$U, xlab="k = 1...K", ylab="i = 1...n")

####### TESTING ALGO
ncomp=10

alpha01 = matrix(1, nrow=n, ncol=ncomp)
alpha02 = matrix(1, nrow=n, ncol=ncomp)

beta01 = matrix(1, nrow=p, ncol=ncomp)
beta02 = matrix(1, nrow=p, ncol=ncomp)

phi01 = matrix(1, nrow=n, ncol=ncomp)
phi02 = matrix(1, nrow=n, ncol=ncomp)

theta01 = matrix(1, nrow=p, ncol=ncomp)
theta02 = matrix(1, nrow=p, ncol=ncomp)

res1 = matrixFactor(data1$X, ncomp, phi01, phi02, theta01, theta02, alpha01, alpha02, beta01, beta02, iterMax=1000, epsilon=1e-4, algo="EM", verbose=TRUE)

str(res1)

print(res1$criteria_k$kDeviance)
myOrder = res1$order$orderDeviance
print(myOrder)


#### plot norms
resU = data.frame(param=rep("U", times=ncomp*res1$nbIter),
                  iter=rep(1:res1$nbIter, each=ncomp),
                  comp=rep(1:ncomp, times=res1$nbIter),
                  norm=as.vector(t(res1$norms$normU[,myOrder])))

resV = data.frame(param=rep("V", times=ncomp*res1$nbIter),
                  iter=rep(1:res1$nbIter, each=ncomp),
                  comp=rep(1:ncomp, times=res1$nbIter),
                  norm=as.vector(t(res1$norms$normV[,myOrder])))

resPhi1 = data.frame(param=rep("Phi1", times=ncomp*res1$nbIter),
                     iter=rep(1:res1$nbIter, each=ncomp),
                     comp=rep(1:ncomp, times=res1$nbIter),
                     norm=as.vector(t(res1$norms$normPhi1[,myOrder])))

resPhi2 = data.frame(param=rep("Phi2", times=ncomp*res1$nbIter),
                     iter=rep(1:res1$nbIter, each=ncomp),
                     comp=rep(1:ncomp, times=res1$nbIter),
                     norm=as.vector(t(res1$norms$normPhi2[,myOrder])))

resTheta1 = data.frame(param=rep("Theta1", times=ncomp*res1$nbIter),
                       iter=rep(1:res1$nbIter, each=ncomp),
                       comp=rep(1:ncomp, times=res1$nbIter),
                       norm=as.vector(t(res1$norms$normTheta1[,myOrder])))

resTheta2 = data.frame(param=rep("Theta2", times=ncomp*res1$nbIter),
                       iter=rep(1:res1$nbIter, each=ncomp),
                       comp=rep(1:ncomp, times=res1$nbIter),
                       norm=as.vector(t(res1$norms$normTheta2[,myOrder])))

resToPlot = rbind(resU, resV, resPhi1, resPhi2, resTheta1, resTheta2)

library(ggplot2)
g1 = ggplot(resToPlot, aes(x=iter, y=norm)) + facet_grid(. ~ param, scales="free")
g1 = g1 + geom_line(aes(col=factor(comp)))
g1 = g1 + xlab(paste0("iter")) + ylab(paste0("dist")) + scale_y_log10()

g1 = g1 + theme(legend.text=element_text(size=12), legend.title=element_text(size=12), axis.text.x=element_text(size=12), axis.title.y = element_text(size=12, vjust=1.5), axis.text.y= element_text(size=12), axis.title.x=element_text(size=12), plot.title=element_text(size=12), strip.text.x=element_text(size=12), strip.text.y=element_text(size=12)) + theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid.major=element_line(color = "grey90"), panel.grid.minor= element_line(color="grey90")) #, legend.position = "none")

g1


### Alpha and Beta
resAlpha1 = data.frame(param=rep("Alpha1", times=ncomp*res1$nbIter),
                       iter=rep(1:res1$nbIter, each=ncomp),
                       comp=rep(1:ncomp, times=res1$nbIter),
                       norm=as.vector(t(res1$norms$normAlpha1[,myOrder])))

resAlpha2 = data.frame(param=rep("Alpha2", times=ncomp*res1$nbIter),
                       iter=rep(1:res1$nbIter, each=ncomp),
                       comp=rep(1:ncomp, times=res1$nbIter),
                       norm=as.vector(t(res1$norms$normAlpha2[,myOrder])))

resBeta1 = data.frame(param=rep("Beta1", times=ncomp*res1$nbIter),
                      iter=rep(1:res1$nbIter, each=ncomp),
                      comp=rep(1:ncomp, times=res1$nbIter),
                      norm=as.vector(t(res1$norms$normBeta1[,myOrder])))

resBeta2 = data.frame(param=rep("Beta2", times=ncomp*res1$nbIter),
                      iter=rep(1:res1$nbIter, each=ncomp),
                      comp=rep(1:ncomp, times=res1$nbIter),
                      norm=as.vector(t(res1$norms$normBeta2[,myOrder])))

resToPlot = rbind(resAlpha1, resAlpha2, resBeta1, resBeta2)

library(ggplot2)
g1 = ggplot(resToPlot, aes(x=iter, y=norm)) + facet_grid(. ~ param, scales="free")
g1 = g1 + geom_line(aes(col=factor(comp)))
g1 = g1 + xlab(paste0("iter")) + ylab(paste0("dist")) + scale_y_log10()

g1 = g1 + theme(legend.text=element_text(size=12), legend.title=element_text(size=12), axis.text.x=element_text(size=12), axis.title.y = element_text(size=12, vjust=1.5), axis.text.y= element_text(size=12), axis.title.x=element_text(size=12), plot.title=element_text(size=12), strip.text.x=element_text(size=12), strip.text.y=element_text(size=12)) + theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid.major=element_line(color = "grey90"), panel.grid.minor= element_line(color="grey90")) #, legend.position = "none")

g1

# tmp = sapply(1:ncomp, function(k) {
#     lambda = res1$U[,k] %*% t(res1$V[,k])
#     summary(log(lambda))
#     return(sum(data1$X * log(lambda) + data1$X * log(data1$X) - lambda - data1$X))
# })
#
# library(fields)
#
# lambda = res1$U[,1] %*% t(res1$V[,1])
# lambda0 = data1$X
# lambda0[lambda0==0] = 1
# dim(lambda0)
# min(lambda0)
# summary(log(lambda))
# image.plot(log(lambda))
#
# sum(is.na(log(lambda)))
#
# hist(res1$U[,1])
# hist(res1$V[,1])
# hist(lambda)
# sum(lambda==0)
#
# image.plot(data1$X * log(data1$X))
# image.plot(data1$X * log(lambda))
# hist(colSums(data1$X * log(data1$X)))
#
#
# image.plot(data1$X * lambda)
#
# min(lambda0 * lambda)
#
# sum(is.na(data1$X * lambda))
#
# image.plot(log(lambda0 * lambda))
#
# sum(data1$X * log(lambda) + data1$X * log(lambda0) - lambda - data1$X)
#
# min(data1$X * (log(lambda0 * lambda) -1) - lambda)
# max(data1$X * (log(lambda0 * lambda) -1) - lambda)
# hist(data1$X * (log(lambda0 * lambda) -1) - lambda)
#
# hist(log(lambda0 / lambda) -1)
# hist(data1$X * (log(lambda0 / lambda) -1))
# sum((data1$X * (log(lambda0 / lambda) -1))<0)
# sum((data1$X * (log(lambda0 / lambda) -1))>0)
# sum(data1$X * (log(lambda0 / lambda) -1))
#
# sum(data1$X * (log(lambda0 / lambda) -1) - lambda)


# tmp/1E6
#
# tmp = data1$X - lambda
# hist(tmp)
# hist(data1$X)
# hist(lambda)

# plot(res1$criteria_k$kExpVarU)

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
# plot(res1$logLikelihood$elbo, xlab="iteration", ylab="elbo", col="blue", type="l")
#
# ## norm gap
# plot(res1$normGap[-1], xlab="iteration", ylab="normalized gap", col="blue", type="b")
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
