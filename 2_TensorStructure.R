# Tensor decomposition on food web metrics
# The analysis is carried out on the median of the resamplings
# The name of species have been changed to S1-S114

# Author: R. Frelat, last update: 08.01.2021

## 1. Load dataset ----------------------------------------
# Load needed packages
library(PTAk)
library(corrplot)
library(RColorBrewer)
library(igraph)
library(ade4)
library(abind)

# Load data 
load("TensorNorthSea.Rdata")

## Add abundance for non sampled taxa
nonsampled <- c("Detritus",  "Microalgae", "Macroalgae",
                "Phytoplankton", "Zooplankton")  

#add constant abundance for non-sampled species
constabu<-array(mean(TensorNS), dim = c(length(nonsampled), dim(TensorNS)[-1]),
                dimnames = list(nonsampled, dimnames(TensorNS)[[2]],
                                dimnames(TensorNS)[[3]]))

TensorNS <- abind(TensorNS, constabu, along=1)
dim(TensorNS) #119 species, 15 years, 6 boxes

# match the species order
TensorNS <- TensorNS[match(V(netNS)$name, dimnames(TensorNS)[[1]]),,]
# dimnames(TensorNS)[[1]]==V(netNS)$name

## Plot food web
# similar to figure 1b
# Average abundance per species
abuM <- apply(TensorNS, 1, mean)
# scale abundance
nsizea<-10
nsizeb<-5
scaleabu<-abuM/max(abuM)*nsizea + nsizeb
sizeN<-scaleabu[match(V(netNS)$name, names(scaleabu))]

op <- par(no.readonly = TRUE)
par(mar=c(1,4,1,1))
coonet <- plotfwtl(netNS, size = sizeN, nylevel = 6, retcoo = TRUE)
par(op) #reset graphical parameters

#Compute metrics on the metaweb
metamet <- fwind(netNS, ab = abuM)
print(metamet)

## 2. Compute food web metrics ----------------------------
fwmet<-array(NA,dim = c(length(metamet), dim(TensorNS)[2:3]))
#The loop might take some minutes to compute
for (i in 1:dim(TensorNS)[2]){
  for (j in 1:dim(TensorNS)[3]){
      abI <- TensorNS[,i,j]
      spI <- dimnames(TensorNS)[[1]][abI>0]
      abI <- abI[abI>0]
      netI <- delete_vertices(netNS, V(netNS)$name[!V(netNS)$name %in% spI])
      #Compute indicators
      tmp <- fwind(netI, ab=abI)
      fwmet[,i,j]<-as.numeric(tmp)
  }
} 
dimnames(fwmet)<-list(names(metamet),
                      dimnames(TensorNS)[[2]],
                      dimnames(TensorNS)[[3]])


# 3. Tensor decomposition ---------------------------------

# 3a Pre-processing 
dim(fwmet) #16 metrics, 15 years, 6 boxes
nar<-dim(fwmet)[3]
nye<-dim(fwmet)[2]
nme<-dim(fwmet)[1]

# correlation plot
# similar to Figure S8
matfwM <- apply(fwmet, 1, c)
corrplot(cor(matfwM), type = "upper", method = "ellipse",
         diag = FALSE, tl.col = "black")

# Scaling
scaleFW<-array(0,dim=dim(fwmet))
for (i in 1:dim(fwmet)[1]){
  ma<-mean(fwmet[i,,])
  sa<-sd(fwmet[i,,])
  scaleFW[i,,]<-(fwmet[i,,]-ma)/sa
}
dimnames(scaleFW)<-dimnames(fwmet)

# 3b Run PTA and select relevant PTs
PTAfw<-PTA3(scaleFW, nbPT = 4, nbPT2 = 3, minpct = 0.1)
# See the summary
summary(PTAfw)

# Select the relevant PTs (identified with a '*')
out <- !substr(PTAfw[[3]]$vsnam, 1, 1) == "*"

# Get the percentage of variance explained by each PT
gct<-(PTAfw[[3]]$pct*PTAfw[[3]]$ssX/PTAfw[[3]]$ssX[1])[out]

# Barplot of the successive explained variance
barplot(sort(gct, decreasing = TRUE),
        border=NA, xlab="", ylab="Percentage of variance")

# Select the number of PTs to be kept
nkeep <- 4 
# Identify the selected PTs and their contribution
keepFW <- (1:length(out))[out][order(gct, decreasing = TRUE)[1:nkeep]]
eigFW <- gct[keepFW]


# 4. Description of selected PT --------------------
# similar to figure 4

# customize the labels
cooFW<-t(PTAfw[[1]]$v[c(keepFW),])
row.names(cooFW)<-dimnames(fwmet)[[1]]
labperc <- round((100 * (PTAfw[[3]]$d[keepFW])^2)/PTAfw[[3]]$ssX[1],1)
labkeep <- paste0(paste0("PT", 1:nkeep), " - ", labperc, "%")
labbox<-gsub("BOX ", "", dimnames(fwmet)[[3]])
topyr<-seq(2000,2015, by=5)
labyr<-ifelse(dimnames(fwmet)[[2]]%in%topyr, dimnames(fwmet)[[2]], "")

#inverse PT sign to simplify interpretation
PTAfw[[1]]$v[keepFW[4],]<- (-1)*PTAfw[[1]]$v[keepFW[4],]
PTAfw[[2]]$v[keepFW[4],]<- (-1)*PTAfw[[2]]$v[keepFW[4],]
PTAfw[[3]]$v[keepFW[4],]<- (-1)*PTAfw[[3]]$v[keepFW[4],]

#create plot for each PT
par(mfrow=c(2,nkeep), mar = c(3,3,3,1), oma=c(1,1,0,0))
for (i in seq_along(keepFW)){
  lab <- labkeep[i]
  temp <- PTAfw[[3]]$v[keepFW[i],] %o% PTAfw[[2]]$v[keepFW[i],]
  dimnames(temp) <- list(labbox, labyr)
  temp<-temp[rev(1:nrow(temp)),]
  # if (i ==4) temp <- -1 * temp
  myHeatmap(temp, pal="BrBG", title=lab, colscale = FALSE, 
            mary = -3,cex.x = 1, cex.y = 1, rm.empty = TRUE,
            tck.x = -0.04, tck.y = -0.04, padj.x = -0.5, hadj.y = 0.5) 
  if (i%%4==1)
    mtext("Box", side = 2, line = 2, xpd=NA, cex=0.8)
  if(i%/%4==1)
    mtext("Year", side = 1, line = 2, xpd=NA, cex=0.7)
  axis(1, at = 1:nye,tck = -0.02, labels = FALSE)
}
for (i in 1:nkeep){
  pal <- brewer.pal(9,"BrBG")
  minmax <- max(abs(cooFW[,i]))
  bk <- seq(-minmax, minmax, length=10)
  coli <- cut(cooFW[,i], breaks = bk, labels = pal,
              include.lowest = TRUE)
  coli <- as.character(coli)
  ordi <- order(cooFW[,i])
  barp <- barplot(cooFW[ordi,i], horiz = TRUE, 
                  col=coli[ordi], 
                  names.arg = rep("", length(ordi)))
  #las=1, names.arg = row.names(coi$li)[ordi])
  text((-0.15)*sign(cooFW[ordi,i]), barp, row.names(cooFW)[ordi], 
       cex=0.7, font=2)
  mtext("Loadings", side = 1, line = 2, xpd=NA, cex=0.7)
}
par(op) #reset graphical parameters

