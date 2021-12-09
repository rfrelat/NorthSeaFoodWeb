# Co-inertia of tensor decomposition
# The analysis is carried out on the median of the resamplings
# The name of species have been changed to S1-S114
# The food web metrics is loaded from the file PTA_results.Rdata
# which is the output from the script 2_TensorStructure.R

# Author: R. Frelat, last update: 08.01.2021

## 1. Load dataset ----------------------------------------
# Load needed packages
require(ade4)

# Load data 
load("TensorNorthSea.Rdata") 
load("PTA_results.Rdata")

nar<-dim(TensorNS)[3]
nye<-dim(TensorNS)[2]

## 2. Transform data into dudi object ---------------------
# for composition - scaleT
info <- expand.grid(dimnames(TensorNS)[2:3])
tab2D <- as.data.frame(apply(scaleT, 1, c))
names(tab2D) <- dimnames(scaleT)[[1]]
dudiPTAABO <- list(
  "tab"= tab2D,
  "rank"=89,
  "cw"=rep(1, nrow(scaleT)),
  "c1"=t(PTAcompo[[1]]$v[c(keepCompo),]),
  "lw"=rep(1/(nar*nye), nar*nye),
  "nf"=length(keepCompo),
  "eig"=eigCompo
)
class(dudiPTAABO) <- "dudi"
names(dudiPTAABO$tab) <- dimnames(scaleT)[[1]]

# for structure - scaleFW
tab2D <- as.data.frame(apply(scaleFW, 1, c))
names(tab2D) <- dimnames(scaleFW)[[1]]
dudiPTAFW <- list(
  "tab"= tab2D,
  "rank"=20,
  "cw"=rep(1, nrow(scaleFW)),
  "c1"=t(PTAfw[[1]]$v[c(keepFW),]),
  "lw"=rep(1/(nar*nye), nar*nye),
  "nf"=length(keepFW),
  "eig"=eigFW
)
class(dudiPTAFW) <- "dudi"

## 3. Co-intertia analysis --------------------------------
# Run co-inertia analysis
coi <- coinertia(dudiPTAABO, dudiPTAFW, scannf = FALSE, nf = 3)

# default plot
plot(coi)

# RV test, checking the significance of co-variance
RV.rtest(dudiPTAABO$tab, dudiPTAFW$tab, nrepet = 1000)

# Number of PC to be view
nkeep <- 2

# customize the labels
labperc <- round(coi$eig/sum(coi$eig)*100, 1)[1:nkeep]
labkeep <- paste0(paste0("PC", 1:nkeep), " - ", labperc, "%")
labbox<-gsub("BOX ", "", dimnames(TensorNS)[[3]])
topyr<-seq(2000,2015, by=5)
labyr<-ifelse(dimnames(TensorNS)[[2]]%in%topyr, dimnames(TensorNS)[[2]], "")

#create plot for each PC
op <- par(no.readonly = TRUE)
par(mfrow=c(3,nkeep), mar = c(3,3,3,1), oma=c(1,1,0,0))
for (i in 1:nkeep){
  lab <- labkeep[i]
  PCi <- t(matrix(c(coi$lX[,i]), ncol=nlevels(info$Var2)))
  #dimnames(PCi) <- list(levels(info$Var2), levels(info$Var1))
  dimnames(PCi) <- list(labbox, labyr)
  PCi<-PCi[rev(1:nrow(PCi)),]
  # if (i ==4) temp <- -1 * temp
  myHeatmap(PCi, pal="BrBG", title=lab, colscale = FALSE, 
            mary = -3,cex.x = 1, cex.y = 1, rm.empty = TRUE,
            tck.x = -0.04, tck.y = -0.04, padj.x = -0.5, hadj.y = 0.5) 
  if (i%%4==1)
    mtext("Box", side = 2, line = 2, xpd=NA, cex=0.8)
  if(i%/%4==1)
    mtext("Year", side = 1, line = 2, xpd=NA, cex=0.7)
  axis(1, at = 1:nye,tck = -0.02, labels = FALSE)
}
par(mar = c(3,4,1,1))
for (i in 1:nkeep){
  boxplot(coi$co[,i]~clustCompo, col=colo, las=1, 
          horizontal = TRUE, ylab="")
  abline(h=0, lty=2, col="grey")
  if (i==1){
    mtext("Species clusters", side = 2, line = 4, xpd=NA, cex=0.8)
  }
}
par(mar = c(2,2,0,1))
for (i in 1:nkeep){
  pal <- brewer.pal(9,"BrBG")
  minmax <- max(abs(coi$li[,i]))
  bk <- seq(-minmax, minmax, length=10)
  coli <- cut(coi$li[,i], breaks = bk, labels = pal,
              include.lowest = TRUE)
  coli <- as.character(coli)
  ordi <- order(coi$li[,i])
  barp <- barplot(coi$li[ordi,i], horiz = TRUE, 
                  col=coli[ordi])
  #las=1, names.arg = row.names(coi$li)[ordi])
  text((-1)*sign(coi$li[ordi,i]), barp, row.names(coi$li)[ordi], 
       cex=0.7, font=2)
  if (i==1){
    mtext("Food web", side = 2, line = 2, xpd=NA, cex=0.8)
  }
}
par(op)

## 4. Co-inertia per box ----------------------------------
yr <- levels(info$Var1)
coiBox <- list()
rvBox <- c()
# loop for each box (3rd dimension)
for (i in 1:dim(TensorNS)[3]){
  abuI <- t(TensorNS[,,i])
  fwI <- t(fwmet[,,i])
  
  # remove species absent from box i
  abuI <- abuI[,apply(abuI,2,sum)>0]
  # run PCA on abundance of box i
  pcaAI <- dudi.pca(abuI, scannf = FALSE, nf = 3)
  # run PCA on food web metrics of box i
  pcaFI <- dudi.pca(fwI, scannf = FALSE, nf = 3)
  # run co-intertia analysis
  coI <- coinertia(pcaAI, pcaFI, scannf = FALSE, nf=2)
  
  # change sign to help comparison of results
  # temporal scores have preferable a positive slope
  if (lm(coI$lX[,1]~as.numeric(yr))$coefficients[2]<0){
    coI$lX[,1] <- (-1)*coI$lX[,1]
    coI$lY[,1] <- (-1)*coI$lY[,1]
    coI$li[,1] <- (-1)*coI$li[,1]
    coI$co[,1] <- (-1)*coI$co[,1]
  }
  # save result of co-intertia analysis in a list
  coiBox[[i]] <- coI
  
  # RV test
  rvI <- RV.rtest(pcaAI$li, pcaFI$li, nrepet = 1000)
  rvBox <- c(rvBox, rvI$pvalue)
}

# similar to figure 6A
par(mar=c(3,4,1,1))
plot(yr, coiBox[[i]]$lX[,1], ylim=c(-9, 8),
     ylab="PC1", type="n")
for (i in 1:length(coiBox)){
  lines(yr, coiBox[[i]]$lX[,1], type="l",
        col=rainbow(6)[i])
}
legend("topleft", legend = dimnames(TensorNS)[[3]],
       lty=1, col=rainbow(6), bty="n")
par(op)

# similar to figure 6B
met <- c()
for (i in 1:length(coiBox)){
  met <- cbind(met, coiBox[[i]]$li[,1])
}
met <- as.data.frame(met)
row.names(met) <- row.names(coiBox[[i]]$li)
colnames(met) <- dimnames(abuM)[[3]]

myHeatmap(met, pal="BrBG")
par(op)
