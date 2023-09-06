library(vegan)
library(bipartite)
library(RInSp)
library(ape)
library(phangorn)

df =as.matrix(dataset[,-1])
row.names(df) = t(dataset[,1])
data = as.phyDat(df, type="USER", levels = c("0", "1"), , ambiguity="-")

# cumulative cultural evolution
## Matrix temperature
out <- nestedtemp(df)
out
plot(out,names = TRUE)
plot(out,names = TRUE, kind="incid",ylab = 'NNs', xlab = 'features',cex = 0.5, main = 'Plot of Temperature')
plot(nestednodf(dataset),names = TRUE)


plot(out,names = TRUE,kind="temperature",cex = 0.5, col=rev(heat.colors(100)),ylab = 'NNs', xlab = 'features',main = 'Plot of Temperature')
plot(nestednodf(df), col = "red", names = TRUE,main = 'Plot of Temperature')

par(mfrow=c(1,1))
a = as.matrix(dataset)
#The entire matrix order is rearranged
random_T <- matrix(sample(a), nrow = nrow(a))
#Rearrange without changing the number of 1s in each row
random_b <- t(apply(a, 1, sample))
#Rearrange without changing the number of 1s in each colume
random_c <- (apply(a, 2, sample))


out_T = nestedtemp(random_T)
out_T 
plot(out_T)
plot(out_T, kind="incid",ylab = 'NNs',xlab = 'features',main = 'Plot of Temperature in random order-T')

out_b = nestedtemp(random_b)
out_b
plot(out_b)
plot(out_b, kind="incid",ylab = 'NNs',xlab = 'features',main = 'Plot of Temperature in random order-R')

out_c = nestedtemp(random_c)
out_c
plot(out_c)
plot(out_c, kind="incid",ylab = 'NNs',xlab = 'features',main = 'Plot of Temperature in random order-C')


plot(nestednodf(df), names= TRUE)



dataset <- import.RInSp(df)
result <- NODF(dataset)

random_NODF <- oecosimu(df, nestednodf, "r00", nsimul = 1000,alternative ="greater")
random_NODF
random_T <- oecosimu(df, nestedtemp, "r00", nsimul = 1000,alternative ="less")
random_T 





#####Phylogenetic tree

# Distance based methods tree-Jaccard

dm1  <- dist(df, method = "binary")
treeNJ1  <- NJ(dm1)

par(mar=c(1, 1, 1, 1))
plot(treeNJ1, "unrooted", cex = 0.4,adj = 0.3,lwd=0.5)

# bootstrap
fun1 <- function(x) NJ(dm1)
bs_NJ1 <- bootstrap.phyDat(data,  fun1)
plotBS(treeNJ1, bs_NJ1,'unrooted', cex = 0.5, main="NJ - Jaccard")


setEPS()

postscript("NJ tree_Jaccard.eps", width = 50, height = 20)
plot(treeNJ1, "unrooted", cex = 0.8,adj = 0.3)
dev.off()


#neighbor joining tree
nnet1 <- neighborNet(dm1)
par("mar" = rep(1, 4))
plot(nnet1, show.edge.label=FALSE,cex=0.5,adj = 0.5, main ='NeighborNet - Jaccard' )


# bootstrap
fun1 <- function(x) NJ(dm1)
bs_NJ1 <- bootstrap.phyDat(data,  fun1)
plotBS(bs_NJ1,'unrooted', cex = 0.5, main="NeighborNet - Jaccard")



# Distance based methods tree-euclidean

dm2  <- dist(df, method = "euclidean")
treeNJ2  <- NJ(dm2)

par(mar=c(1, 1, 1, 1))
plot(treeNJ2, "unrooted", cex = 0.5, main="NJ - Euclidean")

# bootstrap
fun2 <- function(x) NJ(dm2)
bs_NJ2 <- bootstrap.phyDat(data,  fun2)
plotBS(treeNJ2, bs_NJ2,'unrooted', cex = 0.5, main="NJ - Euclidean")

nnet2 <- neighborNet(dm2)
par("mar" = rep(1, 4))
plot(nnet2)
nnet3 <- addConfidences(nnet3, tree)
plot(nnet2, show.edge.label=FALSE,cex=0.5,main ='NeighborNet' )


# Distance based methods tree-Hamming
dm3  <- dist.hamming(data)
treeNJ3  <- NJ(dm3)
plot(treeNJ3 ,cex = 0.5, pch = '-', "unrooted", main="NJ - Hamming")

# bootstrap
fun3 <- function(x) NJ(dm3)
bs_NJ3 <- bootstrap.phyDat(data,  fun3)
plotBS(treeNJ3, bs_NJ3,'unrooted', cex = 0.5, main="NJ - Hamming")


nnet3 <- neighborNet(dm3)
par("mar" = rep(1, 4))
plot(nnet3)
nnet3 <- addConfidences(nnet3, tree)
par("mar" = rep(1, 4))
plot(nnet3, show.edge.label=TRUE,adj = 0.5, cex=0.4,main ='NeighborNet' )




setEPS()
postscript("output.eps", width = 20, height = 20)
plot(NJ(dist.hamming(data,ratio = TRUE)),'unrooted',lwd=0.5)
dev.off()


# bootstrap
fun <- function(x) upgma(dist.hamming(data))
bs_upgma <- bootstrap.phyDat(data,  fun)
plotBS(treeUPGMA, bs_upgma,main="UPGMA")

fun <- function(x) NJ(dist(df, method = "binary"))
bs_NJ <- bootstrap.phyDat(data,  fun)
plotBS(treeNJ, bs_NJ,'unrooted', cex = 0.5, main="NJ")



# neighborNet

cnet <- nnls.networx(cnet, dm2)
par("mar" = rep(1, 4))
plot(cnet, show.edge.label=TRUE)

############################
# use hierarchical clustering to generate phylogenic tree
library(cluster)
m <- c("average","single","complete")
names(m) <- c("average","single","complete")
# function to compute coefficient
ac <- function(x){
  agnes(dm, method = x)$ac
}
library(tidyverse)
map_dbl(m,ac)


dd <- dist.hamming(data)
hc <- hclust(dm2, method = "complete")
hc
dendros <- as.dendrogram(hc)
plot(dendros, main = "NNs - Complete linkage")

library(factoextra)

fviz_nbclust(df, hcut, method = "silhouette")+
  labs(title = "Hierarchical")

fviz_dend(hc, k = 4, k_colors = "jco",main = "NNs - Complete linkage",type = "phylogenic")



# JS linkage
dataset2 = t(dataset)
df2 = as.matrix(dataset2[-1,])
col.names(df2) = t(dataset[1,])

dm0 <- dist(df2, method = "binary")
dm_JS = 1-dm0
dm_JS = as.matrix(dm_JS)

write.csv(dm_JS,"C:/Users/PC/Desktop/dm_JS.csv")

library(reshape2)
library(ggplot2)


longData<-melt(dm_JS)
longData<-longData[longData$value!=0,]

p = ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red")

p = p+theme(axis.text.x = element_text(
  angle = 45,hjust = 1,vjust = 1))

p


df_dependency = as.matrix(dependency[,-1])

row.names(df_dependency) = t(dependency[,1])
longData2<-melt(df_dependency)
longData2<-longData2[longData2$value!=0.5,]

p2 = ggplot(longData2, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red")

p2 = p2+theme(axis.text.x = element_text(
  angle = 45,hjust = 1,vjust = 1))

p2

