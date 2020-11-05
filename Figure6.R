### Figure 6a

library(openxlsx)
bands <- read.xlsx("tables/Supplementary Table 16.xlsx", startRow = 4, sheet = 5)
bands <- bands[c(1:3, 4, 5:7, 24:26), ] #select the 3 most significant bands per group (if available)
bands[,6:13] <- apply(bands[,6:13], 2, as.numeric)

M <- as.matrix(bands[, c("p-value.MTC", "p-value.GPM", "p-value.NEU", "p-value.PPR")])
Mstat <- as.matrix(bands[, c("odds.ratio.MTC", "odds.ratio.GPM", "odds.ratio.NEU", "odds.ratio.PPR")])
rownames(M) <- rownames(Mstat) <- bands$Chromosomal.Band

M2 <- -log10(M)
M2[Mstat < 1] <- -M2[Mstat < 1]

colnames(M2) <- c("MTC", "GPM", "NEU", "PPR")

library(corrplot)
library(RColorBrewer)
par(mar=c(8, 4.1, 8, 2.1))
corrplot(M2, method="circle", is.corr=F, p.mat = M, sig.level = c(0.02, 0.05, 0.1), insig = "label_sig", pch.cex = 2, pch.col = "white",
         col=colorRampPalette(c("#8D0000",  "white", "blue4"))(7), tl.srt = 45, tl.col = "black", tl.cex = 1, cl.lim = c(-2, 2))


### Figure 6b
seg <- read.table("RData/gbm.seg", stringsAsFactors = F, header = T)
ss <- unique(seg$Sample)
ss <- ss[substr(ss, 14, 16) == "01A"]
ss <- ss[!ss %in% ss[duplicated(substr(ss, 1, 12))]]
seg <- seg[seg$Sample %in% ss, ]
seg$Sample <- substr(seg$Sample, 1, 12)

load("RData/TCGA_geData.RData", verbose = T)

ss <- intersect(unique(seg$Sample), colnames(geData))
seg <- seg[seg$Sample %in% ss, ]
seg <- as.data.frame(seg)
seg <- cbind(seg, Segment_Length = (seg$End - seg$Start + 1))

karyotypeBand <- read.table("RData/karyotypeBandExons.txt", stringsAsFactors = F, header = T, sep = "\t")
karyotypeBand <- karyotypeBand[karyotypeBand$Gene.type == "protein_coding", ]
karyotypeBand <- karyotypeBand[-grep("GL|H|MT", karyotypeBand$Chromosome.scaffold.name), ]
karyotypeBand <- cbind(karyotypeBand, Gene.Length = (karyotypeBand$Gene.end..bp. - karyotypeBand$Gene.start..bp. + 1))


ssamples <- sort(unique(seg$Sample))
ggenes <- sort(unique(karyotypeBand$Gene.name))

focAmp <- matrix(0, nrow = length(ggenes), ncol = length(ssamples))
rownames(focAmp) <- ggenes
colnames(focAmp) <- ssamples
recAmp <- focDel <- recDel <- focAmp

library(doMC)
C <- 64
registerDoMC(C)
ans <- foreach(j = 1:length(ssamples)) %dopar% {
  message(paste(j, "of", length(ssamples)))
  tmpSeg <- seg[seg$Sample == ssamples[j], ]
  tmpSeg <- tmpSeg[abs(tmpSeg$Segment_Mean) > 0.3, ]
  
  tPos <- abs(sum(tmpSeg$Segment_Mean[tmpSeg$Segment_Mean > 0]))
  tNeg <- abs(sum(tmpSeg$Segment_Mean[tmpSeg$Segment_Mean < 0]))
  
  fAmp <- fDel <- rAmp <- rDel <- NULL
  for(i in 1:nrow(tmpSeg)){
    tmpKB <- karyotypeBand[karyotypeBand$Chromosome.scaffold.name == tmpSeg$Chromosome[i], ]
    tmpKB <- tmpKB[tmpKB$Exon.region.start..bp. <= tmpSeg$End[i] & tmpKB$Exon.region.end..bp. >= tmpSeg$Start[i], ]
    tmpKB <- tmpKB[!duplicated(tmpKB$Gene.name), ]
    if(nrow(tmpKB) == 0) next
    
    ansF <- abs(tmpSeg$Segment_Mean[i])/nrow(tmpKB)
    ansF <- rep(ansF, nrow(tmpKB))
    names(ansF) <- tmpKB$Gene.name
    if(tmpSeg$Segment_Mean[i] > 0) fAmp <- c(fAmp, ansF) else
      fDel <- c(fDel, ansF)
    
    if(tmpSeg$Segment_Mean[i] > 0){
      ansR <- abs(tmpSeg$Segment_Mean[i])/tPos
      ansR <- rep(ansR, nrow(tmpKB))
      names(ansR) <- tmpKB$Gene.name
      rAmp <- c(rAmp, ansR)
    } else {
      ansR <- abs(tmpSeg$Segment_Mean[i])/tNeg
      ansR <- rep(ansR, nrow(tmpKB))
      names(ansR) <- tmpKB$Gene.name
      rDel <- c(rDel, ansR)
    }
  }
  
  if(length(fAmp) > 0){
    fAmp <- tapply(fAmp, names(fAmp), mean)
    rAmp <- tapply(rAmp, names(rAmp), mean)
  }
  if(length(fDel) > 0){
    fDel <- tapply(fDel, names(fDel), mean)
    rDel <- tapply(rDel, names(rDel), mean)
  }
  
  ans <- list(fAmp = fAmp, rAmp = rAmp, fDel = fDel, rDel = rDel)
  return(ans)
  
}
names(ans) <- ssamples


for(i in 1:length(ans)){
  focAmp[names(ans[[i]]$fAmp), i] <- as.numeric(ans[[i]]$fAmp)
  recAmp[names(ans[[i]]$rAmp), i] <- as.numeric(ans[[i]]$rAmp)
  focDel[names(ans[[i]]$fDel), i] <- as.numeric(ans[[i]]$fDel)
  recDel[names(ans[[i]]$rDel), i] <- as.numeric(ans[[i]]$rDel)
  print(i)
}

focAmpScore <- rowSums(focAmp)
recAmpScore <- rowSums(recAmp)
focDelScore <- rowSums(focDel)
recDelScore <- rowSums(recDel)
rawAmpScore <- focAmpScore * recAmpScore
rawDelScore <- focDelScore * recDelScore
ampScore <- rawAmpScore/sum(rawAmpScore)
delScore <- rawDelScore/sum(rawDelScore)
focAmpRank <- recAmpRank <- ampRank <- focDelRank <- recDelRank <- delRank <- rep(NA, length(ampScore))
focAmpRank[focAmpScore != 0] <- length(focAmpScore[focAmpScore != 0]) - rank(focAmpScore[focAmpScore != 0]) + 1
recAmpRank[recAmpScore != 0] <- length(recAmpScore[recAmpScore != 0]) - rank(recAmpScore[recAmpScore != 0]) + 1
ampRank[ampScore != 0] <- length(ampScore[ampScore != 0]) - rank(ampScore[ampScore != 0]) + 1
focDelRank[focDelScore != 0] <- length(focDelScore[focDelScore != 0]) - rank(focDelScore[focDelScore != 0]) + 1
recDelRank[recDelScore != 0] <- length(recDelScore[recDelScore != 0]) - rank(recDelScore[recDelScore != 0]) + 1
delRank[delScore != 0] <- length(delScore[delScore != 0]) - rank(delScore[delScore != 0]) + 1

cfScores <- data.frame(focAmpScore, focAmpRank, recAmpScore, recAmpRank, rawAmpScore, ampScore, ampRank, focDelScore, focDelRank, recDelScore, recDelRank, rawDelScore, delScore, delRank)

karyotypeBand <- read.table("RData/karyotypeBand.txt", stringsAsFactors = F, header = T, sep = "\t")
karyotypeBand <- karyotypeBand[karyotypeBand$Gene.type == "protein_coding", ]
karyotypeBand <- karyotypeBand[-grep("GL|H|MT", karyotypeBand$Chromosome.scaffold.name), ]
karyotypeBand <- karyotypeBand[!duplicated(karyotypeBand$Gene.name), ]
rownames(karyotypeBand) <- karyotypeBand$Gene.name

band <- paste0(karyotypeBand[rownames(cfScores), "Chromosome.scaffold.name"], karyotypeBand[rownames(cfScores), "Karyotype.band"])
cfScores <- cbind(band = band, cfScores, stringsAsFactors = F)

x <- cfScores[cfScores$band == "1p36.23", ]
x <- x[order(x$delRank), ]
x

gg <- c("VAMP3", "PER3", "UTS2", "TNFRSF9", "PARK7", "ERRFI1", "SLC45A1", "RERE", "ENO1", "CA6", "SLC2A7", "SLC2A5", "GPR157")
gg <- rev(gg)

#left panel
cfScores <- cfScores[gg, ]
plot(cfScores$delScore, 1:nrow(cfScores), type = "l", yaxt = "n", xlab = "", ylab = "")
points(cfScores$delScore, 1:nrow(cfScores), pch = 16)

#lower panel
load("RData/TCGA_CNVfunc.RData", verbose = T)
load("RData/TCGA_CNVthres.RData", verbose = T)
CNVthres <- CNVthres[rownames(CNVfunc), colnames(CNVfunc)]

karyotypeBand <- read.table("RData/karyotypeBand.txt", stringsAsFactors = F, header = T, sep = "\t")
karyotypeBand <- karyotypeBand[karyotypeBand$Gene.type == "protein_coding", ]
karyotypeBand <- karyotypeBand[-grep("GL|H|MT", karyotypeBand$Chromosome.scaffold.name), ]
karyotypeBand <- karyotypeBand[!duplicated(karyotypeBand$Gene.name), ]
rownames(karyotypeBand) <- karyotypeBand$Gene.name
band <- paste0(karyotypeBand[, "Chromosome.scaffold.name"], karyotypeBand[, "Karyotype.band"])
karyotypeBand <- cbind(karyotypeBand, band = band, stringsAsFactors = F)

karyotypeBand <- karyotypeBand[karyotypeBand$band == "1p36.23", ]
karyotypeBand <- karyotypeBand[order(karyotypeBand$Gene.start..bp., decreasing = T), ]

ggenesToConsider <- rownames(karyotypeBand)
CNVthres <- CNVthres[ggenesToConsider, ]
CNVfunc <- CNVfunc[ggenesToConsider, ]

CNVthres[CNVfunc == -1] <- -3
CNVthres[CNVfunc == -2] <- -4
CNVthres[CNVthres > 0] <- 0

CNVthres <- CNVthres[, colSums(CNVthres != 0) > 0]
CNVthres <- CNVthres[, order(CNVthres["SLC45A1", ])]


library(heatmap3)
heatmap3(CNVthres, showRowDendro = F, showColDendro = F,
         Rowv = NA,
         Colv = NA,
         col = colorRampPalette(brewer.pal(n=9, name="Blues")[c(9, 1)])(5),
         labCol = NA,
         scale = 'none', margins = c(8, 5))


### Figure 6e
# upload the dataframe relative to the differential expression analysis between MTC compared to the GPM
# the dataframe should include a column called "log2FC" in which log2FC between two conditions is reported

load("RData/sc_dataset1_MTCvsGPM_DEGs.RData")
dataset1_FC <- DEGs_MTCvsGPM_dataset1[, "FC"]
load("RData/sc_dataset2_MTCvsGPM_DEGs.RData")
dataset2_FC <- DEGs_MTCvsGPM_dataset2[, "FC"]
load("RData/sc_dataset3_MTCvsGPM_DEGs.RData")
dataset3_FC <- DEGs_MTCvsGPM_dataset3[, "FC"]
load("RData/TCGA_MTCvsGPM_DEGs.RData")
dataset_TCGA <- DEGs_MTCvsGPM_TCGA[, "FC"]

# upload the dataframe relative significantly differential expressed genes (see methods for thresholds used) between the subtype of interest compared to the others
# the rownames of dataframe should be the gene names of significantly differential expressed genes (see supplementary table 8c)
library(openxlsx)
SuppTab3 <- read.xlsx("tables/Supplementary Table 8.xlsx", sheet = 3)

sigDEGs_datasetTCGA <- SuppTab3[SuppTab3$X6 == "TCGA",]
colnames(sigDEGs_datasetTCGA) <- c("Gene", "Gene family", "fold change", 	"p-value", "FDR", "Dataset")
sigDEGs_datasetTCGA <- sigDEGs_datasetTCGA[-1, ]
rownames(sigDEGs_datasetTCGA) <- sigDEGs_datasetTCGA$Gene

sigDEGs_dataset1 <- SuppTab3[SuppTab3$X6 == "Dataset 1",]
colnames(sigDEGs_dataset1) <- c("Gene", "Gene family", "fold change", 	"p-value", "FDR", "Dataset")
sigDEGs_dataset1 <- sigDEGs_dataset1[-1, ]
rownames(sigDEGs_dataset1) <- sigDEGs_dataset1$Gene

sigDEGs_dataset2 <-SuppTab3[SuppTab3$X6 == "Dataset 2",]
colnames(sigDEGs_dataset2) <- c("Gene", "Gene family", "fold change", 	"p-value", "FDR", "Dataset")
sigDEGs_dataset2 <- sigDEGs_dataset2[-1, ]
rownames(sigDEGs_dataset2) <- sigDEGs_dataset2$Gene

sigDEGs_dataset3 <-SuppTab3[SuppTab3$X6 == "Dataset 3",]
colnames(sigDEGs_dataset3) <- c("Gene", "Gene family", "fold change", 	"p-value", "FDR", "Dataset")
sigDEGs_dataset3 <- sigDEGs_dataset3[-1, ]
rownames(sigDEGs_dataset3) <- sigDEGs_dataset3$Gene

# intersect genes which expression has been detected in all three single cell datasets
gg <- intersect(intersect(intersect(names(dataset1_FC), names(dataset2_FC)), names(dataset3_FC)), names(dataset_TCGA))
dataset_TCGA <- dataset_TCGA[gg]
dataset1_FC <- dataset1_FC[gg]
dataset2_FC <- dataset2_FC[gg]
dataset3_FC <- dataset3_FC[gg]

# order the vector of the FC derived for all genes in decreasing order
dataset_TCGA <- dataset_TCGA[order(dataset_TCGA, decreasing = F)]
dataset1_FCord <- dataset1_FC[order(dataset1_FC, decreasing = F)]
dataset2_FCord <- dataset2_FC[order(dataset2_FC, decreasing = F)]
dataset3_FCord <- dataset3_FC[order(dataset3_FC, decreasing = F)]

# generate plots Fig. 6e
plot(x = 1:length(dataset1_FCord), y = dataset1_FCord, col="gray", pch=16, type = 'n', axes = F, ylim = c(-4, 4))
points(which(!names(dataset1_FCord)%in%rownames(sigDEGs_dataset1)), dataset1_FCord[which(!names(dataset1_FCord)%in%rownames(sigDEGs_dataset1))], pch = 16, col="gray78", cex =0.2)
points(which(!names(dataset2_FCord)%in%rownames(sigDEGs_dataset2)), dataset2_FCord[which(!names(dataset2_FCord)%in%rownames(sigDEGs_dataset2))], pch = 16, col="gray38", cex = 0.2)
points(which(!names(dataset3_FCord)%in%rownames(sigDEGs_dataset3)), dataset3_FCord[which(!names(dataset3_FCord)%in%rownames(sigDEGs_dataset3))], pch = 16, col="gray58", cex = 0.2)
points(which(!names(dataset_TCGA)%in%rownames(sigDEGs_datasetTCGA)), dataset_TCGA[which(!names(dataset_TCGA)%in%rownames(sigDEGs_datasetTCGA))], pch = 16, col="black", cex = 0.2)
abline(h = -0.3,  col = "darkgray", lty = 3)
points(which(names(dataset1_FCord)%in%rownames(sigDEGs_dataset1)), dataset1_FCord[which(names(dataset1_FCord)%in%rownames(sigDEGs_dataset1))], pch = 16, col="blue", cex = 0.45)
points(which(names(dataset2_FCord)%in%rownames(sigDEGs_dataset2)), dataset2_FCord[which(names(dataset2_FCord)%in%rownames(sigDEGs_dataset2))], pch = 16, col="blue", cex = 0.45)
points(which(names(dataset3_FCord)%in%rownames(sigDEGs_dataset3)), dataset3_FCord[which(names(dataset3_FCord)%in%rownames(sigDEGs_dataset3))], pch = 16, col="blue", cex = 0.45)
points(which(names(dataset_TCGA)%in%rownames(sigDEGs_datasetTCGA)), dataset_TCGA[which(names(dataset_TCGA)%in%rownames(sigDEGs_datasetTCGA))], pch = 16, col="blue", cex = 0.45)
axis(side = 2, at= -4:4, labels = -4:4)
