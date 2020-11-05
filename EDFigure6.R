### ED Figure 6a
# upload the gene signatures in TCGA subtypes
library(openxlsx)

geneSig <- read.xlsx("tables/Supplementary Table 6.xlsx", sheet = 10)
geneSig <- geneSig[2:nrow(geneSig), ]
colnames(geneSig) <- geneSig[1, ]
geneSig <- geneSig[-1, ]
geneSigGPM <- geneSig$GPM
geneSigMTC <- geneSig$MTC
geneSigPPR <- geneSig$PPR
geneSigNEU <- geneSig$NEU

# upload the TCGA-RNAseq classification
library(openxlsx)
TCGAclassification <- read.xlsx("tables/Supplementary Table 10.xlsx", sheet = 1)
TCGAclassification <- TCGAclassification[2:nrow(TCGAclassification), ]
colnames(TCGAclassification) <- TCGAclassification[1, ]
TCGAclassification <- TCGAclassification[-1, ]
rownames(TCGAclassification) <- TCGAclassification$`TCGA ID`
TCGAclassification2 <- TCGAclassification$`GBM subtype`
names(TCGAclassification2) <- rownames(TCGAclassification)
table(TCGAclassification2)
TCGAclassification2 <- c(TCGAclassification2[TCGAclassification2 == "NEU"], TCGAclassification2[TCGAclassification2 == "PPR"], TCGAclassification2[TCGAclassification2 == "MTC"], TCGAclassification2[TCGAclassification2 == "GPM"])

# upload the TCGA-RNAseq classification color annotation
colClass <- c("red", "green3", "cyan", "blue")
names(colClass) <- c("GPM", "MTC", "PPR", "NEU")
colClass <- colClass[TCGAclassification2]
names(colClass) <- names(TCGAclassification2)
ColSideColors = colClass
colnames(ColSideColors) <- NULL

load("RData/exprsMatrixRNAseq_TCGA.RData", verbose = T)

# create a list with the gene signatures in each subtype and color annotation
listOfDiffs <- vector("list", 4)
names(listOfDiffs) <- c("GPM", "MTC", "PPR", "NEU")
listOfDiffs$GPM <- geneSigGPM
listOfDiffs$MTC <-  geneSigMTC
listOfDiffs$PPR <-  geneSigPPR
listOfDiffs$NEU <-  geneSigNEU

gg <- intersect(as.character(unlist(listOfDiffs)), rownames(geData_TCGARNAseq))
listOfDiffs <- lapply(listOfDiffs, function(x) x[x %in% gg])

listOfDiffs <- lapply(listOfDiffs, function(x) x)
tmp <- lapply(listOfDiffs, function(x) length(x))
aa <- rep("red", tmp$GPM)
bb <- rep("green3", tmp$MTC)
cc <- rep("cyan", tmp$PPR)
dd <- rep("blue", tmp$NEU)
names(aa) <- listOfDiffs$GPM
names(bb) <- listOfDiffs$MTC
names(cc) <- listOfDiffs$PPR
names(dd) <- listOfDiffs$NEU

colGenes <- c(aa, bb, cc, dd)
table(colGenes)
annot3 <- cbind(colGenes)
colnames(annot3) <- c("colPATH")
rownames(annot3) <- as.character(unlist(listOfDiffs))
RowSideColors = annot3

geData_TCGARNAseq <- geData_TCGARNAseq[gg, names(TCGAclassification2)]
toPlot <- geData_TCGARNAseq 

#plot the ED. Fig. 6a
toPlotStand <- toPlot
meanEachGene <- apply(toPlot, 1, function(x) mean(x))
stDevEachGene <- apply(toPlot, 1, function(x) sd(x))
for (i in 1:nrow(toPlot)){
  gene <- toPlot[i,]
  toPlotStand[i,] <- (gene-meanEachGene[i])/stDevEachGene[i]
}
toPlot <- toPlotStand
summary(as.vector(toPlot))
quantile(toPlot, 0.15); quantile(toPlot, 0.85)
toPlot[toPlot <= quantile(toPlot, 0.15)] <- quantile(toPlot, 0.15)
toPlot[toPlot >= quantile(toPlot, 0.85)] <- quantile(toPlot, 0.85)

geneExpr <- t(toPlot)

library(heatmap3)
heatmap3(t(geneExpr), showRowDendro = F, showColDendro = F,
         Rowv = NA,
         Colv = NA,
         ColSideColors = ColSideColors, ColSideLabs = NA,
         RowSideColors = RowSideColors, RowSideLabs = NA,
         col = colorRampPalette(c("#053061", "#4393C3", "cornsilk", "#D6604D", "#67001F"))(75),
         labCol = NA,
         cexCol = 0.5,
         cexRow = 0.3,
         scale = 'none', useRaster = F)



### ED Figure 6b
# upload the gene signatures in TCGA subtypes
library(openxlsx)

geneSig <- read.xlsx("tables/Supplementary Table 6.xlsx", sheet = 10)
geneSig <- geneSig[2:nrow(geneSig), ]
colnames(geneSig) <- geneSig[1, ]
geneSig <- geneSig[-1, ]
geneSigGPM <- geneSig$GPM
geneSigMTC <- geneSig$MTC
geneSigPPR <- geneSig$PPR
geneSigNEU <- geneSig$NEU

# upload the CGGA datasaset classification
library(openxlsx)
CGGAclassification <- read.xlsx("tables/Supplementary Table 10.xlsx", sheet = 2)
CGGAclassification <- CGGAclassification[2:nrow(CGGAclassification), ]
colnames(CGGAclassification) <- CGGAclassification[1, ]
CGGAclassification <- CGGAclassification[-1, ]
rownames(CGGAclassification) <- CGGAclassification$`CGGA ID`
CGGAclassification2 <- CGGAclassification$`GBM subtype`
names(CGGAclassification2) <- rownames(CGGAclassification)
table(CGGAclassification2)
CGGAclassification2 <- c(CGGAclassification2[CGGAclassification2 == "NEU"], CGGAclassification2[CGGAclassification2 == "PPR"], CGGAclassification2[CGGAclassification2 == "MTC"], CGGAclassification2[CGGAclassification2 == "GPM"])

# upload the CGGA datasaset classification color annotation
colClass <- c("red", "green3", "cyan", "blue")
names(colClass) <- c("GPM", "MTC", "PPR", "NEU")
colClass <- colClass[CGGAclassification2]
names(colClass) <- names(CGGAclassification2)
ColSideColors = colClass
colnames(ColSideColors) <- NULL


load("RData/exprsMatrixRNAseq_CGGA.RData", verbose = T)

# create a list with the gene signatures in each subtype and color annotation
listOfDiffs <- vector("list", 4)
names(listOfDiffs) <- c("GPM", "MTC", "PPR", "NEU")
listOfDiffs$GPM <- geneSigGPM
listOfDiffs$MTC <-  geneSigMTC
listOfDiffs$PPR <-  geneSigPPR
listOfDiffs$NEU <-  geneSigNEU

gg <- intersect(as.character(unlist(listOfDiffs)), rownames(geData_CGGA))
listOfDiffs <- lapply(listOfDiffs, function(x) x[x %in% gg])

listOfDiffs <- lapply(listOfDiffs, function(x) x)
tmp <- lapply(listOfDiffs, function(x) length(x))
aa <- rep("red", tmp$GPM)
bb <- rep("green3", tmp$MTC)
cc <- rep("cyan", tmp$PPR)
dd <- rep("blue", tmp$NEU)
names(aa) <- listOfDiffs$GPM
names(bb) <- listOfDiffs$MTC
names(cc) <- listOfDiffs$PPR
names(dd) <- listOfDiffs$NEU

colGenes <- c(aa, bb, cc, dd)
table(colGenes)
annot3 <- cbind(colGenes)
colnames(annot3) <- c("colPATH")
rownames(annot3) <- as.character(unlist(listOfDiffs))
RowSideColors = annot3

geData_CGGA <- geData_CGGA[gg, names(CGGAclassification2)]
toPlot <- geData_CGGA 

#plot the ED. Fig. 6b
toPlotStand <- toPlot
meanEachGene <- apply(toPlot, 1, function(x) mean(x))
stDevEachGene <- apply(toPlot, 1, function(x) sd(x))
for (i in 1:nrow(toPlot)){
  gene <- toPlot[i,]
  toPlotStand[i,] <- (gene-meanEachGene[i])/stDevEachGene[i]
}
toPlot <- toPlotStand
summary(as.vector(toPlot))
quantile(toPlot, 0.15); quantile(toPlot, 0.85)
toPlot[toPlot <= quantile(toPlot, 0.15)] <- quantile(toPlot, 0.15)
toPlot[toPlot >= quantile(toPlot, 0.85)] <- quantile(toPlot, 0.85)

geneExpr <- t(toPlot)

library(heatmap3)
heatmap3(t(geneExpr), showRowDendro = F, showColDendro = F,
         Rowv = NA,
         Colv = NA,
         ColSideColors = ColSideColors, ColSideLabs = NA,
         RowSideColors = RowSideColors, RowSideLabs = NA,
         col = colorRampPalette(c("#053061", "#4393C3", "cornsilk", "#D6604D", "#67001F"))(75),
         labCol = NA,
         cexCol = 0.5,#3,
         cexRow = 0.3,
         scale = 'none', useRaster = F)


### ED Figure 6c
# upload the gene signatures in TCGA subtypes
library(openxlsx)

geneSig <- read.xlsx("tables/Supplementary Table 6.xlsx", sheet = 10)
geneSig <- geneSig[2:nrow(geneSig), ]
colnames(geneSig) <- geneSig[1, ]
geneSig <- geneSig[-1, ]
geneSigGPM <- geneSig$GPM
geneSigMTC <- geneSig$MTC
geneSigPPR <- geneSig$PPR
geneSigNEU <- geneSig$NEU

# upload the Lee et al. datasaset classification
library(openxlsx)
Leeclassification <- read.xlsx("tables/Supplementary Table 10.xlsx", sheet = 3)
Leeclassification <- Leeclassification[2:nrow(Leeclassification), ]
colnames(Leeclassification) <- Leeclassification[1, ]
Leeclassification <- Leeclassification[-1, ]
rownames(Leeclassification) <- Leeclassification$`GEO accession number`
Leeclassification2 <- Leeclassification$`GBM subtype`
names(Leeclassification2) <- rownames(Leeclassification)
table(Leeclassification2)
Leeclassification2 <- c(Leeclassification2[Leeclassification2 == "NEU"], Leeclassification2[Leeclassification2 == "PPR"], Leeclassification2[Leeclassification2 == "MTC"], Leeclassification2[Leeclassification2 == "GPM"])

# upload the Lee et al. datasaset classification color annotation
colClass <- c("red", "green3", "cyan", "blue")
names(colClass) <- c("GPM", "MTC", "PPR", "NEU")
colClass <- colClass[Leeclassification2]
names(colClass) <- names(Leeclassification2)
ColSideColors = colClass
colnames(ColSideColors) <- NULL

load("RData/exprsMatrix_LeeEtAl.RData", verbose = T)

# create a list with the gene signatures in each subtype and color annotation
listOfDiffs <- vector("list", 4)
names(listOfDiffs) <- c("GPM", "MTC", "PPR", "NEU")
listOfDiffs$GPM <- geneSigGPM
listOfDiffs$MTC <-  geneSigMTC
listOfDiffs$PPR <-  geneSigPPR
listOfDiffs$NEU <-  geneSigNEU

gg <- intersect(as.character(unlist(listOfDiffs)), rownames(geDataLeeEtAl))
listOfDiffs <- lapply(listOfDiffs, function(x) x[x %in% gg])

listOfDiffs <- lapply(listOfDiffs, function(x) x)
tmp <- lapply(listOfDiffs, function(x) length(x))
aa <- rep("red", tmp$GPM)
bb <- rep("green3", tmp$MTC)
cc <- rep("cyan", tmp$PPR)
dd <- rep("blue", tmp$NEU)
names(aa) <- listOfDiffs$GPM
names(bb) <- listOfDiffs$MTC
names(cc) <- listOfDiffs$PPR
names(dd) <- listOfDiffs$NEU

colGenes <- c(aa, bb, cc, dd)
table(colGenes)
annot3 <- cbind(colGenes)
colnames(annot3) <- c("colPATH")
rownames(annot3) <- as.character(unlist(listOfDiffs))
RowSideColors = annot3

geDataLeeEtAl <- geDataLeeEtAl[gg, names(Leeclassification2)]
toPlot <- geDataLeeEtAl 

#plot the ED. Fig. 6c
toPlotStand <- toPlot
meanEachGene <- apply(toPlot, 1, function(x) mean(x))
stDevEachGene <- apply(toPlot, 1, function(x) sd(x))
for (i in 1:nrow(toPlot)){
  gene <- toPlot[i,]
  toPlotStand[i,] <- (gene-meanEachGene[i])/stDevEachGene[i]
}
toPlot <- toPlotStand
summary(as.vector(toPlot))
quantile(toPlot, 0.15); quantile(toPlot, 0.85)
toPlot[toPlot <= quantile(toPlot, 0.15)] <- quantile(toPlot, 0.15)
toPlot[toPlot >= quantile(toPlot, 0.85)] <- quantile(toPlot, 0.85)

geneExpr <- t(toPlot)

library(heatmap3)
heatmap3(t(geneExpr), showRowDendro = F, showColDendro = F,
         Rowv = NA,
         Colv = NA,
         ColSideColors = ColSideColors, ColSideLabs = NA,
         RowSideColors = RowSideColors, RowSideLabs = NA,
         col = colorRampPalette(c("#053061", "#4393C3", "cornsilk", "#D6604D", "#67001F"))(75),
         labCol = NA,
         cexCol = 0.5,#3,
         cexRow = 0.3,
         scale = 'none', useRaster = F)



### ED Figure 6d
# upload the TCGA-RNAseq classification
library(openxlsx)
TCGAclassification <- read.xlsx("tables/Supplementary Table 10.xlsx", sheet = 1)
TCGAclassification <- TCGAclassification[2:nrow(TCGAclassification), ]
colnames(TCGAclassification) <- TCGAclassification[1, ]
TCGAclassification <- TCGAclassification[-1, ]
rownames(TCGAclassification) <- TCGAclassification$`TCGA ID`
TCGAclassification2 <- TCGAclassification$`GBM subtype`
names(TCGAclassification2) <- rownames(TCGAclassification)
table(TCGAclassification2)
TCGAclassification2 <- c(TCGAclassification2[TCGAclassification2 == "NEU"], TCGAclassification2[TCGAclassification2 == "PPR"], TCGAclassification2[TCGAclassification2 == "MTC"], TCGAclassification2[TCGAclassification2 == "GPM"])
groups <- TCGAclassification2

#upload survival data
load("RData/TCGAGBM_clinicalData.RData", verbose = T)

library(survival)
commonSamples <- intersect(rownames(stData), names(groups))
stData <- stData[commonSamples, ]
groups <- factor(groups[commonSamples])

par(mfrow=c(2,2))
aSurv <- survfit(stData ~ groups)
logrank <- survdiff(stData ~ groups)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groups))-1)), digits = 4))
plot(aSurv, col = c("#8D0000", "forestgreen", "blue2", "cyan2"), lwd = 2.5, font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=2.5, padj = 0.3, at=c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.5", "0.75", "1"))
axis(1,cex.axis=2.5, padj = 0.3)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("GPM", "MTC")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("#8D0000", "forestgreen"), lwd = 2.5, font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=2.5, padj = 0.3, at=c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.5", "0.75", "1"))
axis(1,cex.axis=2.5, padj = 0.3)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("MTC", "NEU")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("forestgreen", "blue2"), lwd = 2.5,  font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=2.5, padj = 0.3, at=c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.5", "0.75", "1"))
axis(1,cex.axis=2.5, padj = 0.3)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("MTC", "PPR")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("forestgreen", "cyan2"), lwd = 2.5, font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=2.5, padj = 0.3, at=c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.5", "0.75", "1"))
axis(1,cex.axis=2.5, padj = 0.3)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)


### ED Figure 6e
# upload the CGGA datasaset classification
library(openxlsx)
CGGAclassification <- read.xlsx("tables/Supplementary Table 10.xlsx", sheet = 2)
CGGAclassification <- CGGAclassification[2:nrow(CGGAclassification), ]
colnames(CGGAclassification) <- CGGAclassification[1, ]
CGGAclassification <- CGGAclassification[-1, ]
rownames(CGGAclassification) <- CGGAclassification$`CGGA ID`
CGGAclassification2 <- CGGAclassification$`GBM subtype`
names(CGGAclassification2) <- rownames(CGGAclassification)
table(CGGAclassification2)
CGGAclassification2 <- c(CGGAclassification2[CGGAclassification2 == "NEU"], CGGAclassification2[CGGAclassification2 == "PPR"], CGGAclassification2[CGGAclassification2 == "MTC"], CGGAclassification2[CGGAclassification2 == "GPM"])
groups <- CGGAclassification2

#upload survival data
library(survival)
load("RData/CGGAGBM_survData.RData")
stData <- survCGGA
commonSamples <- intersect(rownames(stData), names(groups))
stData <- stData[commonSamples, ]
groups <- factor(groups[commonSamples])

par(mfrow=c(2,2))
aSurv <- survfit(stData ~ groups)
logrank <- survdiff(stData ~ groups)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groups))-1)), digits = 4))
plot(aSurv, col = c("#8D0000", "forestgreen", "blue2", "cyan2"), lwd = 2.5, font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=2.5, padj = 0.3, at=c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.5", "0.75", "1"))
axis(1,cex.axis=2.5, padj = 0.3)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("GPM", "MTC")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("#8D0000", "forestgreen"), lwd = 2.5, font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=2.5, padj = 0.3, at=c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.5", "0.75", "1"))
axis(1,cex.axis=2.5, padj = 0.3)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("MTC", "NEU")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("forestgreen", "blue2"), lwd = 2.5,  font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=2.5, padj = 0.3, at=c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.5", "0.75", "1"))
axis(1,cex.axis=2.5, padj = 0.3)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("MTC", "PPR")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("forestgreen", "cyan2"), lwd = 2.5, font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=2.5, padj = 0.3, at=c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.5", "0.75", "1"))
axis(1,cex.axis=2.5, padj = 0.3)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)


### ED Figure 6f
# upload the Lee et al. datasaset classification
library(openxlsx)
Leeclassification <- read.xlsx("tables/Supplementary Table 10.xlsx", sheet = 3)
Leeclassification <- Leeclassification[2:nrow(Leeclassification), ]
colnames(Leeclassification) <- Leeclassification[1, ]
Leeclassification <- Leeclassification[-1, ]
rownames(Leeclassification) <- Leeclassification$`GEO accession number`
Leeclassification2 <- Leeclassification$`GBM subtype`
names(Leeclassification2) <- rownames(Leeclassification)
table(Leeclassification2)
Leeclassification2 <- c(Leeclassification2[Leeclassification2 == "NEU"], Leeclassification2[Leeclassification2 == "PPR"], Leeclassification2[Leeclassification2 == "MTC"], Leeclassification2[Leeclassification2 == "GPM"])
groups <- Leeclassification2

#upload survival data
library(survival)
load("RData/LeeEtAlGBM_survData.RData")
stData <- survLeeEtAl
commonSamples <- intersect(rownames(stData), names(groups))
stData <- stData[commonSamples, ]
groups <- factor(groups[commonSamples])

par(mfrow=c(2,2))
aSurv <- survfit(stData ~ groups)
logrank <- survdiff(stData ~ groups)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groups))-1)), digits = 4))
plot(aSurv, col = c("#8D0000", "forestgreen", "blue2", "cyan2"), lwd = 2.5, font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=2.5, padj = 0.3, at=c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.5", "0.75", "1"))
axis(1,cex.axis=2.5, padj = 0.3, at=c(0, 500, 1000, 1500, 2000), labels=c("0", "", "1000", "", "2000"))
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("GPM", "MTC")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("#8D0000", "forestgreen"), lwd = 2.5, font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=2.5, padj = 0.3, at=c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.5", "0.75", "1"))
axis(1,cex.axis=2.5, padj = 0.3, at=c(0, 500, 1000, 1500, 2000), labels=c("0", "", "1000", "", "2000"))
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("MTC", "NEU")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("forestgreen", "blue2"), lwd = 2.5,  font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=2.5, padj = 0.3, at=c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.5", "0.75", "1"))
axis(1,cex.axis=2.5, padj = 0.3, at=c(0, 500, 1000, 1500, 2000), labels=c("0", "", "1000", "", "2000"))
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("MTC", "PPR")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("forestgreen", "cyan2"), lwd = 2.5, font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=2.5, padj = 0.3, at=c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.5", "0.75", "1"))
axis(1,cex.axis=2.5, padj = 0.3, at=c(0, 500, 1000, 1500, 2000), labels=c("0", "", "1000", "", "2000"))
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)



### ED Figure 6g
# upload the TCGA-RNAseq classification
library(openxlsx)
TCGAclassification <- read.xlsx("tables/Supplementary Table 10.xlsx", sheet = 4)
TCGAclassification <- TCGAclassification[2:nrow(TCGAclassification), ]
colnames(TCGAclassification) <- TCGAclassification[1, ]
TCGAclassification <- TCGAclassification[-1, ]
rownames(TCGAclassification) <- TCGAclassification$`TCGA ID`
TCGAclassification2 <- TCGAclassification$`Wang et al. 2017 classification`
names(TCGAclassification2) <- rownames(TCGAclassification)
table(TCGAclassification2)
groups <- TCGAclassification2

#upload survival data
load("RData/TCGAGBM_clinicalData.RData")
library(survival)
commonSamples <- intersect(rownames(stData), names(groups))
stData <- stData[commonSamples, ]
groups <- factor(groups[commonSamples])

library(survival)
stData <- stData[!is.na(stData[, "time"]), ]
commonSamples <- intersect(rownames(stData), names(groups))
stData <- stData[commonSamples, ]
groups <- groups[commonSamples]

par(mfrow=c(2,2))
aSurv <- survfit(stData ~ groups)
logrank <- survdiff(stData ~ groups)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groups))-1)), digits = 4))
plot(aSurv, col = c("blue", "green3", "purple"), lwd = 2.5, font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=1.5)
axis(1,cex.axis=1.5)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("Classical", "Mesenchymal")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("blue", "green3"), lwd = 2.5, font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=1.5)
axis(1,cex.axis=1.5)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("Classical", "Proneural")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("blue", "purple"), lwd = 2.5,  font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=1.5)
axis(1,cex.axis=1.5)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("Mesenchymal", "Proneural")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("green3", "purple"), lwd = 2.5,  font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=1.5)
axis(1,cex.axis=1.5)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)



### ED Figure 6h
# upload the TCGA-RNAseq classification
library(openxlsx)
TCGAclassification <- read.xlsx("tables/Supplementary Table 10.xlsx", sheet = 1)
TCGAclassification <- TCGAclassification[2:nrow(TCGAclassification), ]
colnames(TCGAclassification) <- TCGAclassification[1, ]
TCGAclassification <- TCGAclassification[-1, ]
rownames(TCGAclassification) <- TCGAclassification$`TCGA ID`
TCGAclassification2 <- TCGAclassification$`Wang et al. 2017 classification`
names(TCGAclassification2) <- rownames(TCGAclassification)
table(TCGAclassification2)
groups <- TCGAclassification2

#upload survival data
load("RData/TCGAGBM_clinicalData.RData")
library(survival)
commonSamples <- intersect(rownames(stData), names(groups))
stData <- stData[commonSamples, ]
groups <- groups[commonSamples]

library(survival)
stData <- stData[!is.na(stData[, "time"]), ]
commonSamples <- intersect(rownames(stData), names(groups))
stData <- stData[commonSamples, ]
groups <- factor(groups[commonSamples])


par(mfrow=c(2,2))
aSurv <- survfit(stData ~ groups)
logrank <- survdiff(stData ~ groups)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groups))-1)), digits = 4))
plot(aSurv, col = c("blue", "green3", "purple"), lwd = 2.5, font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=1.5)
axis(1,cex.axis=1.5)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("Classical", "Mesenchymal")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("blue", "green3"), lwd = 2.5, font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=1.5)
axis(1,cex.axis=1.5)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("Classical", "Proneural")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("blue", "purple"), lwd = 2.5,  font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=1.5)
axis(1,cex.axis=1.5)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("Mesenchymal", "Proneural")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("green3", "purple"), lwd = 2.5,  font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=1.5)
axis(1,cex.axis=1.5)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)



### ED Figure 6i
# upload the TCGA-RNAseq classification
library(openxlsx)
TCGAclassification <- read.xlsx("tables/Supplementary Table 10.xlsx", sheet = 4)
TCGAclassification <- TCGAclassification[2:nrow(TCGAclassification), ]
colnames(TCGAclassification) <- TCGAclassification[1, ]
TCGAclassification <- TCGAclassification[-1, ]
rownames(TCGAclassification) <- TCGAclassification$`TCGA ID`
TCGAclassification2 <- TCGAclassification$`Phillips et al. 2006 classification`
names(TCGAclassification2) <- rownames(TCGAclassification)
table(TCGAclassification2)
groups <- TCGAclassification2

# upload survival data
load("RData/TCGAGBM_clinicalData.RData")
library(survival)
commonSamples <- intersect(rownames(stData), names(groups))
stData <- stData[commonSamples, ]
groups <- groups[commonSamples]

library(survival)
stData <- stData[!is.na(stData[, "time"]), ]
commonSamples <- intersect(rownames(stData), names(groups))
stData <- stData[commonSamples, ]
groups <- factor(groups[commonSamples])

par(mfrow=c(2,2))
aSurv <- survfit(stData ~ groups)
logrank <- survdiff(stData ~ groups)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groups))-1)), digits = 4))
plot(aSurv, col = c("firebrick3", "navy", "forestgreen"), lwd = 2.5, font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=1.5)
axis(1,cex.axis=1.5)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("Mesenchymal", "Proliferative")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("firebrick3", "navy"), lwd = 2.5, font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=1.5)
axis(1,cex.axis=1.5)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("Proliferative", "Proneural")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("navy", "forestgreen"), lwd = 2.5,  font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=1.5)
axis(1,cex.axis=1.5)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("Mesenchymal", "Proneural")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("firebrick3", "forestgreen"), lwd = 2.5,  font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=1.5)
axis(1,cex.axis=1.5)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)


### ED Figure 6j
# upload the TCGA-RNAseq classification
library(openxlsx)
TCGAclassification <- read.xlsx("tables/Supplementary Table 10.xlsx", sheet = 1)
TCGAclassification <- TCGAclassification[2:nrow(TCGAclassification), ]
colnames(TCGAclassification) <- TCGAclassification[1, ]
TCGAclassification <- TCGAclassification[-1, ]
rownames(TCGAclassification) <- TCGAclassification$`TCGA ID`
TCGAclassification2 <- TCGAclassification$`Phillips et al. 2006 classification`
names(TCGAclassification2) <- rownames(TCGAclassification)
table(TCGAclassification2)
groups <- TCGAclassification2

# upload survival data
load("RData/TCGAGBM_clinicalData.RData")
library(survival)
commonSamples <- intersect(rownames(stData), names(groups))
stData <- stData[commonSamples, ]
groups <- groups[commonSamples]

library(survival)
stData <- stData[!is.na(stData[, "time"]), ]
commonSamples <- intersect(rownames(stData), names(groups))
stData <- stData[commonSamples, ]
groups <- factor(groups[commonSamples])


par(mfrow=c(2,2))
aSurv <- survfit(stData ~ groups)
logrank <- survdiff(stData ~ groups)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groups))-1)), digits = 4))
plot(aSurv, col = c("firebrick3", "navy", "forestgreen"), lwd = 2.5, font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=1.5)
axis(1,cex.axis=1.5)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("Mesenchymal", "Proliferative")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("firebrick3", "navy"), lwd = 2.5, font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=1.5)
axis(1,cex.axis=1.5)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("Proliferative", "Proneural")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("navy", "forestgreen"), lwd = 2.5,  font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=1.5)
axis(1,cex.axis=1.5)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("Mesenchymal", "Proneural")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("firebrick3", "forestgreen"), lwd = 2.5,  font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=1.5)
axis(1,cex.axis=1.5)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)
