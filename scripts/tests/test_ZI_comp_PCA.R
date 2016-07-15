### pca on ZI vs non ZI data


source("/home/durif/source_code/countMatrixFactor/set_working_dir.R")
source("sources/lightSources.R")


######################
### Cpp version
######################

## generating the data
n = 100
p = 10
K = 4

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

X = data1$X

## heatmap
matrixHeatmap(alpha1, xlab="k = 1...K", ylab="i = 1...n")
matrixHeatmap(beta1, xlab="k = 1...K", ylab="j = 1...p")
matrixHeatmap(data1$X, xlab="j = 1...p", ylab="i = 1...n")
matrixHeatmap(data1$U, xlab="k = 1...K", ylab="i = 1...n")

# png(filename=paste0(WD, "count_matrix.png"), width=180, height=120, res=500, units="mm")
# matrixHeatmap(X)
# dev.off()

## ZI model

ncomp=2

alpha01 = matrix(1, nrow=n, ncol=ncomp)
alpha02 = matrix(1, nrow=n, ncol=ncomp)

beta01 = matrix(1, nrow=p, ncol=ncomp)
beta02 = matrix(1, nrow=p, ncol=ncomp)

phi01 = matrix(1, nrow=n, ncol=ncomp)
phi02 = matrix(1, nrow=n, ncol=ncomp)

theta01 = matrix(1, nrow=p, ncol=ncomp)
theta02 = matrix(1, nrow=p, ncol=ncomp)

res1a = matrixFactor(X, ncomp,
                     phi01, phi02, theta01, theta02,
                     alpha01, alpha02, beta01, beta02,
                     iterMax=500, epsilon=1e-4,
                     ZI=FALSE, algo = "variational")

str(res1a)

U1a = res1a$U
V1a = res1a$V


## PCA

res1b = prcomp(X)
str(res1b)
V1b = res1b$rotation[,1:ncomp]
U1b = X %*% V1b
# plot(U1b, col=blockAlpha1$idRows)
# plot(V1b, col=blockBeta1$idRows)
# matrixHeatmap(U1, xlab="k = 1...K", ylab="i = 1...n")

### ZI
prob0 = round(runif(p, 0.3, 0.6), digits=2)
Y = sapply(prob0, function(p) return(rbinom(n=n,size=1,prob=p)))

X_ZI = X * Y

res2a = matrixFactor(X_ZI, ncomp,
                     phi01, phi02, theta01, theta02,
                     alpha01, alpha02, beta01, beta02,
                     iterMax=500, epsilon=1e-4,
                     ZI=TRUE, algo = "variational")

str(res2a)

U2a = res2a$U
V2a = res2a$V

res2b = prcomp(X_ZI)
str(res2b)
V2b = res2b$rotation[,1:ncomp]
U2b = X_ZI %*% V2b


### graph
U1as = scale(U1a)
U1bs = scale(U1b)
U2as = scale(U2a)
U2bs = scale(U2b)

dataToPlot1a = data.frame(model = rep("non ZI", n), method = rep("varinf", n), comp1 = U1as[,1], comp2 = U1as[,2], group=blockAlpha1$idRows)
dataToPlot1b = data.frame(model = rep("non ZI", n), method = rep("pca", n), comp1 = U1bs[,1], comp2 = U1bs[,2], group=blockAlpha1$idRows)
dataToPlot2a = data.frame(model = rep("ZI", n), method = rep("varinf", n), comp1 = U2as[,1], comp2 = U2as[,2], group=blockAlpha1$idRows)
dataToPlot2b = data.frame(model = rep("ZI", n), method = rep("pca", n), comp1 = U2bs[,1], comp2 = U2bs[,2], group=blockAlpha1$idRows)

dataToPlot = rbind(dataToPlot1a, dataToPlot1b, dataToPlot2a, dataToPlot2b)
dataToPlot$group = as.factor(dataToPlot$group)

g1 = ggplot(dataToPlot, aes(x=comp1, y=comp2, col=group)) + geom_point() + facet_grid(model ~ method)

g1 = g1 + theme(legend.text=element_text(size=10), legend.title=element_text(size=10), axis.text.x=element_text(size=10), axis.title.y = element_text(size=10, vjust=1.5), axis.text.y= element_text(size=10), axis.title.x=element_text(size=10), plot.title=element_text(size=10), strip.text.x=element_text(size=10), strip.text.y=element_text(size=10)) + theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid.major=element_line(color = "grey90"), panel.grid.minor= element_line(color="grey90"), legend.position = "none")

g1

cumsum(res2$sdev)/sum(res2$sdev)

ggsave(filename = paste0(path, "/error_pred_box_", data_set, ".eps"), plot = g1, width=234, height=80, dpi=1200, scale=1, units="mm")

#### graphs
dataToPlot1$group = as.factor(dataToPlot1$group)
e1 = round(100*res1$sdev/sum(res1$sdev), 2)
g1 = ggplot(dataToPlot1, aes(x=comp1, y=comp2, col=group)) + geom_point()
g1 = g1 + xlab(paste0("Comp 1 (", e1[1], "%)")) + ylab(paste0("Comp 2 (", e1[2], "%)"))

g1 = g1 + theme(legend.text=element_text(size=10), legend.title=element_text(size=10), axis.text.x=element_text(size=10), axis.title.y = element_text(size=10, vjust=1.5), axis.text.y= element_text(size=10), axis.title.x=element_text(size=10), plot.title=element_text(size=10), strip.text.x=element_text(size=10), strip.text.y=element_text(size=10)) + theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid.major=element_line(color = "grey90"), panel.grid.minor= element_line(color="grey90"), legend.position = "none")

g1

ggsave(filename = paste0(WD, "pca_non_zi.eps"), plot = g1, width=100, height=100, dpi=1200, scale=1, units="mm")



dataToPlot2$group = as.factor(dataToPlot2$group)
e2 = round(100*res2$sdev/sum(res2$sdev), 2)
g2 = ggplot(dataToPlot2, aes(x=comp1, y=comp2, col=group)) + geom_point()
g2 = g2 + xlab(paste0("Comp 1 (", e2[1], "%)")) + ylab(paste0("Comp 2 (", e2[2], "%)"))

g2 = g2 + theme(legend.text=element_text(size=10), legend.title=element_text(size=10), axis.text.x=element_text(size=10), axis.title.y = element_text(size=10, vjust=1.5), axis.text.y= element_text(size=10), axis.title.x=element_text(size=10), plot.title=element_text(size=10), strip.text.x=element_text(size=10), strip.text.y=element_text(size=10)) + theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid.major=element_line(color = "grey90"), panel.grid.minor= element_line(color="grey90"), legend.position = "none")

g2

ggsave(filename = paste0(WD, "pca_zi.eps"), plot = g2, width=100, height=100, dpi=1200, scale=1, units="mm")


#### MODEL
source("/home/durif/source_code/countMatrixFactor/set_working_dir.R")
source("sources/projectSources.R")

## TESTING ALGO
ncomp=2

alpha01 = matrix(1, nrow=n, ncol=ncomp)
alpha02 = matrix(1, nrow=n, ncol=ncomp)

beta01 = matrix(1, nrow=p, ncol=ncomp)
beta02 = matrix(1, nrow=p, ncol=ncomp)

phi01 = matrix(1, nrow=n, ncol=ncomp)
phi02 = matrix(1, nrow=n, ncol=ncomp)

theta01 = matrix(1, nrow=p, ncol=ncomp)
theta02 = matrix(1, nrow=p, ncol=ncomp)

res3 = matrixFactor(X, ncomp,
                    phi01, phi02, theta01, theta02,
                    alpha01, alpha02, beta01, beta02,
                    iterMax=200, epsilon=1e-4,
                    ZI=FALSE, algo = "variational")

str(res3)
myOrder1 = res3$order$orderExpVarV

U3 = res3$U[, myOrder1]
plot(U3, col=blockAlpha1$idRows)

res4 = matrixFactor(X_ZI, ncomp,
                    phi01, phi02, theta01, theta02,
                    alpha01, alpha02, beta01, beta02,
                    iterMax=200, epsilon=1e-4,
                    ZI=FALSE, algo = "variational")

str(res4)
myOrder1 = res4$order$orderExpVarV

U4 = res4$U[, myOrder1]
plot(U4, col=blockAlpha1$idRows)

### graph
U3s = scale(U3)
U4s = scale(U4)
dataToPlot3 = data.frame(method = rep("non ZI", n), comp1 = U3s[,1], comp2 = U3s[,2], group=blockAlpha1$idRows)
dataToPlot4 = data.frame(method = rep("ZI", n), comp1 = U4s[,1], comp2 = U4s[,2], group=blockAlpha1$idRows)


dataToPlot3$group = as.factor(dataToPlot3$group)
g1 = ggplot(dataToPlot3, aes(x=comp1, y=comp2, col=group)) + geom_point()

g1 = g1 + theme(legend.text=element_text(size=10), legend.title=element_text(size=10), axis.text.x=element_text(size=10), axis.title.y = element_text(size=10, vjust=1.5), axis.text.y= element_text(size=10), axis.title.x=element_text(size=10), plot.title=element_text(size=10), strip.text.x=element_text(size=10), strip.text.y=element_text(size=10)) + theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid.major=element_line(color = "grey90"), panel.grid.minor= element_line(color="grey90"), legend.position = "none")

g1

ggsave(filename = paste0(WD, "varinf_non_zi.eps"), plot = g1, width=100, height=100, dpi=1200, scale=1, units="mm")



dataToPlot4$group = as.factor(dataToPlot4$group)
g2 = ggplot(dataToPlot4, aes(x=comp1, y=comp2, col=group)) + geom_point()

g2 = g2 + theme(legend.text=element_text(size=10), legend.title=element_text(size=10), axis.text.x=element_text(size=10), axis.title.y = element_text(size=10, vjust=1.5), axis.text.y= element_text(size=10), axis.title.x=element_text(size=10), plot.title=element_text(size=10), strip.text.x=element_text(size=10), strip.text.y=element_text(size=10)) + theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid.major=element_line(color = "grey90"), panel.grid.minor= element_line(color="grey90"), legend.position = "none")

g2

ggsave(filename = paste0(WD, "varinf_zi.eps"), plot = g2, width=100, height=100, dpi=1200, scale=1, units="mm")
