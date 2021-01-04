# Tensor decomposition on fish and benthos communities

# analysis is carried out on the median of the resamplings
# the name of species have been altered to A:XX


## Load dataset -----------------------------------------
# Load needed packages
library(PTAk)
library(RColorBrewer)

# Load data 
load("TensorNorthSea.Rdata")

#Check the diensions of the tensor
dim(TensorNS) #114, 15, 6
nsp<-dim(TensorNS)[1]
nye<-dim(TensorNS)[2]
nar<-dim(TensorNS)[3]


## Scaling -----------------------------------------------
scaleT<-array(0,dim=dim(TensorNS))
for (i in 1:dim(TensorNS)[1]){
  ma<-mean(TensorNS[i,,])
  sa<-sd(TensorNS[i,,])
  scaleT[i,,]<-(TensorNS[i,,]-ma)/sa
}
dimnames(scaleT)<-dimnames(TensorNS)


## Principal tensor ------------------------------
#1 Simple PCA 3D with only key tensor
PTAcompo<-PTA3(scaleT, nbPT = 4, nbPT2 = 3, minpct = 0.1)
summary(PTAcompo)

summary.PTAk(PTAcompo,testvar = 0)
#First principal tensor explain 28% of variability (v111)
out <- !substr(PTAcompo[[3]]$vsnam, 1, 1) == "*"
gct<-(PTAcompo[[3]]$pct*PTAcompo[[3]]$ssX/PTAcompo[[3]]$ssX[1])[out]

barplot(sort(gct, decreasing = TRUE),
        border=NA, xlab="", ylab="Percentage of variance")

nkeep <- 5 #X1. 7, X2. 6
keepCompo <- (1:length(out))[out][order(gct, decreasing = TRUE)[1:nkeep]]
eigCompo <- gct[keepCompo]


# Description of selected PT --------------------
# similar to figure 2
labperc <- round((100 * (PTAcompo[[3]]$d[keepCompo])^2)/PTAcompo[[3]]$ssX[1],1)
labkeep <- paste0(paste0("PT", 1:nkeep), " - ", labperc, "%")
labbox<-gsub("BOX ", "", dimnames(scaleT)[[3]])
topyr<-seq(2000,2015, by=5)
labyr<-ifelse(dimnames(scaleT)[[2]]%in%topyr, dimnames(scaleT)[[2]], "")

op <- par(no.readonly = TRUE)
par(mfrow=c(2,3), mar = c(3,3,3,1), oma=c(1,1,0,0))
for (i in seq_along(keepCompo)){
  lab <- labkeep[i]
  temp <- PTAcompo[[3]]$v[keepCompo[i],] %o% PTAcompo[[2]]$v[keepCompo[i],]
  dimnames(temp) <- list(labbox, labyr)
  temp<-temp[rev(1:nrow(temp)),]
  # if (i ==4) temp <- -1 * temp
  myHeatmap(temp, pal="BrBG", title=lab, colscale = FALSE, 
            mary = -3,cex.x = 1, cex.y = 1, rm.empty = TRUE,
            tck.x = -0.04, tck.y = -0.04, padj.x = -0.5, hadj.y = 0.5) 
  if (i%%3==1)
    mtext("Box", side = 2, line = 2, xpd=NA, cex=0.8)
  if(i>2)
    mtext("Year", side = 1, line = 2, xpd=NA, cex=0.7)
  axis(1, at = 1:nye,tck = -0.02, labels = FALSE)
}
par(mar=c(4,7,4,7))
pal <- brewer.pal(8, "BrBG")
plot(0, xlim=c(0,1), ylim=c(0.1,length(pal)-0.1), type="n", 
     xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
for (i in 1:length(pal)){
  rect(0,i-1,1,i,col = pal[i])
}
axis(2, at=c(0, length(pal)/2, length(pal)), labels = c("Low", "0", "High"))
#mtext("Loadings", 3, line = 0.5, adj = 1)
mtext("Loadings", 2, line=2)
par(op) #reset graphical parameters

## Classification -------------------------------------------
# classify species based on their scores
# get Figure SX and 3

# Get the score of the species on the selected PTs
cooCompo <- t(PTAcompo[[1]]$v[c(keepCompo),])

# Compute the Euclidean distance between species' score 
distCoo <- dist(cooCompo, method = "euclidean")

# Create hierarchical clustering based on Ward's criteria
den <- hclust(distCoo,method = "ward.D2")

# Plot the dendogram
plot(den, hang=-1, ax = T, ann=F, xlab="", sub="",labels = FALSE)

# Select 6 clusters
nclust<-6

# Visualize the clusters on the dendogram
rect.hclust(den, k=nclust)

# Create the clusters
clustCompo <- as.factor(cutree(den, k=nclust))

#Alternative way to visualize the number of clusters
par(mar=c(4,4,1,1))
barp<-barplot(sort(den$height, decreasing = TRUE)[1:30],
              ylab="height", xlab="# clusters",
              col=c(rep("black", nclust-1), rep("grey", 31-nclust)))
axis(1, barp, labels = 2:(length(barp)+1))
par(op) #reset graphical parameters


colo <- brewer.pal(nclust,"Set1")
levels(clustCompo) <- c("clustA", "clustM", "clustM+",
                      "clustC", "clustD", "clustA+")
clustCompo <- factor(clustCompo, levels = levels(clustCompo)[c(6,1,4,5,2,3)])

## Interpret the clusters -----------------------
# reconstruct the dynamics from the selected PTs
PTDyn <- REBUILD(PTAcompo, nTens = keepCompo, testvar = 1)

# compute the average dynamics per cluster
meanDyn<-array(0, dim=c(nclust, nye, nar))
for (c in 1:nclust){
  meanDyn[c,,]<-apply(PTDyn[as.numeric(clustCompo)==c,,],c(2,3),mean)
}
dimnames(meanDyn)<-list(1:nclust, dimnames(scaleT)[[2]], dimnames(scaleT)[[3]]) 

# Set the color scale
ab <- max(abs(meanDyn))
colcut <- seq(-ab, ab, length.out = 9)
colpal <- brewer.pal(8,"BrBG")

par(mfrow=c(3,nclust+1),oma=c(1,3,0,0))
for (c in 1:nclust){
  temp <- t(meanDyn[c,,])
  if (c==1){
    dimnames(temp) <- list(labbox, labyr)
  } else {
    dimnames(temp) <- list(rep("", nar), labyr)
  }
  temp<-temp[rev(1:nrow(temp)),]
  myHeatmap(temp, colscale = FALSE, pal = colpal, 
            title=levels(clustCompo)[c], mary = -5, breaks = colcut, 
            rm.empty = TRUE, tck.x = -0.07, tck.y = -0.05, 
            cex.y = 1, cex.x = 1, padj.x = -0.8, hadj.y=0.5)
  axis(1, at = 1:nye,tck = -0.03, labels = FALSE)
  axis(2,tck = -0.05, labels = FALSE)
  if (c==1){
    mtext("Box", side = 2, line = 2, xpd=NA, cex=0.8)
  }
}
par(mar=c(1,3,1,1))
plot(0, xlim=c(0,1), ylim=c(0.1,7.9), type="n", 
     xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
for (i in 1:8){
  rect(0,i-1,1,i,col = colpal[i])
}
axis(2,at=c(0.5,4, 7.5), labels = c("low","0","high"), las=1, cex=0.4)
mtext("Anomaly", side = 2, line=2, xpd=NA, cex=0.8)
for (c in 1:nclust){
  par(mar=c(2,0,0,1))
  maploadings(coordinatesBox, apply(meanDyn[c,,],2,mean), 
              xlim=c(-1.5, 8.5), ylim=c(53.5, 61.5),
              xlab="", yaxt="n")
  if (c==1){
    axis(2, tck = -0.05, las=1, cex=0.9, hadj = 0.7)
    mtext("Latitude", side = 2, line = 2, xpd=NA, cex=0.8)
  } else {
    axis(2, tck = -0.05, cex=0.9, hadj = 0.5, label=FALSE)
  }
  axis(1)
}
plot.new()
for (c in 1:nclust){
  par(mar=c(2,0,0,1))
  plot(dimnames(scaleT)[[2]],apply(meanDyn[c,,],1,mean), type = "l", 
       ylim=c(-0.6, 0.5), xaxt="n", yaxt="n")
  axis(1, at = dimnames(scaleT)[[2]], tck = -0.03, labels=FALSE)
  axis(1, at = seq(1990, 2015,5), tck = -0.05, cex=0.9, padj = -0.7)
  if (c==1){
    axis(2, tck = -0.05, las=1, cex=0.9, hadj = 0.7)
  } else {
    axis(2, tck = -0.05, cex=0.9, hadj = 0.5, label=FALSE)
  }
  abline(h=0, lty=3, col="grey")
  mtext("Year", side = 1, line = 1.5, xpd = NA, cex = 0.8)
  if (c==1){
    mtext("Anomaly", side = 2, line = 2, xpd=NA, cex=0.8)
  }
}
par(op) #reset graphical parameters

