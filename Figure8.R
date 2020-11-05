### Figure 8e
#UPPER PANEL
# upload the classification of PDCs dataset
library(openxlsx)
PDCsclassification <- read.xlsx("tables/Supplementary Table 17.xlsx", sheet = 9)
PDCsclassification <- PDCsclassification[3:nrow(PDCsclassification), ]
colnames(PDCsclassification) <- PDCsclassification[1, ]
PDCsclassification <- PDCsclassification[-1, ]
rownames(PDCsclassification) <- PDCsclassification$`PDC ID`
PDCsclassification2 <- PDCsclassification$`Random Forest classification`
names(PDCsclassification2) <- rownames(PDCsclassification)

#upload NES activity GPM and MTC
NES_GPM <- as.numeric(PDCsclassification$`GPM NES`)
names(NES_GPM) <- rownames(PDCsclassification)
NES_MTC <- as.numeric(PDCsclassification$`MTC NES`)
names(NES_MTC) <- rownames(PDCsclassification)

# derive the difference of GPM activity versus MTC activity
diffGPM_MTCact <- NES_GPM - NES_MTC
diffGPM_MTCact_ord <- sort(diffGPM_MTCact, decreasing = F)

#derive the combined mitochondrial inhibitor sensitivity score
scoreMTCinh <- as.numeric(PDCsclassification$`Combined MTC inhibitor sensitivity score`)
names(scoreMTCinh) <- rownames(PDCsclassification)
scoreMTCinh <- scoreMTCinh[names(diffGPM_MTCact_ord)]

#Generate upper panel Fig. 8e 
par(mar=c(0,0,0,0))
plot(1:length(scoreMTCinh), scoreMTCinh,  type = "n", xlim = c(0.6, length(scoreMTCinh)+0.6), xaxs = "i", xaxt = "n", axes = F)
abline(lm(scoreMTCinh ~ order(sort(diffGPM_MTCact_ord, decreasing = F))), lwd = 2)
points(which(PDCsclassification2=="MTC"), scoreMTCinh[names(which(PDCsclassification2=="MTC"))], col = "blue4", pch = 16, cex=2)
points(which(PDCsclassification2=="GPM"), scoreMTCinh[names(which(PDCsclassification2=="GPM"))], col = "darkred", pch = 16, cex=2)

#Generate middle panel Fig. 8e 
plot(NES_MTC, type = "n", ylim=c(-3.5, 5), xaxt = "n", axes = F)
plot(NES_MTC, type = "n", ylim=c(-3.5, 5))
points(NES_MTC, col = "blue4", pch = 16, cex=1.3)
points(NES_GPM, col = "darkred", pch = 16, cex=1.3)
abline(lm(NES_MTC ~ order(sort(diffGPM_MTCact, decreasing = F))))
abline(lm(NES_GPM ~ order(sort(diffGPM_MTCact, decreasing = F))))
rect(0.39, min(NES_MTC)-0.63, length(NES_GPM)+0.50, max(NES_GPM)+0.11,
     col=NULL, border=par("fg"), lty=NULL, lwd=par("lwd"), xpd=FALSE)

#Generate lower panel Fig. 8e 

#upload genes identified as glycolytic and mitochondrial associated to GPM and MTC PDCs subtypes
load("RData/AMPDEL_GlycoMitoch_PDCs.RData")
load("RData/CNVfunc_PDCs.RData")

CNVfunc2AMP <- CNVfunc_PDCs[c(AMP_mitochgenes,AMP_glycolgenes), ]
table(CNVfunc2AMP)
CNVfunc2AMP[CNVfunc2AMP == -2] <- 0
CNVfunc2AMP[CNVfunc2AMP == -1 ] <- 0
CNVfunc2AMP[CNVfunc2AMP == 2] <- 1
CNVfunc2AMP <- CNVfunc2AMP[, names(diffGPM_MTCact_ord)]

CNVfunc2DEL <- CNVfunc_PDCs[c(DEL_mitochgenes,DEL_glycolgenes), ]
table(CNVfunc2DEL)
CNVfunc2DEL[CNVfunc2DEL == 2] <- 0
CNVfunc2DEL[CNVfunc2DEL == 1 ] <- 0
CNVfunc2DEL[CNVfunc2DEL == -2] <- 1
CNVfunc2DEL[CNVfunc2DEL == -1 ] <- 1
CNVfunc2DEL <- CNVfunc2DEL[, names(diffGPM_MTCact_ord)]

sumAMPmitoch <- colSums(CNVfunc2AMP[AMP_mitochgenes, ])
sumDELmitoch <- colSums(CNVfunc2DEL[DEL_glycolgenes, ])
sumAMPglyco <- colSums(CNVfunc2AMP[AMP_glycolgenes, ])
sumDELglyco <- colSums(CNVfunc2DEL[DEL_mitochgenes, ])

SupportPlot <- CNVfunc_PDCs[c("PTEN","RB1"), names(diffGPM_MTCact_ord)]
table(SupportPlot)
SupportPlot[SupportPlot==1] <- 0
SupportPlot[SupportPlot==-1] <- 1
SupportPlot[SupportPlot==-2] <- 1

toPlot <- SupportPlot

library(gplots)
tmp <- sumAMPmitoch
nCol <- 100
ccol <- colorRampPalette(c("white", "darkred"))(nCol)
ccut <- cut(tmp, nCol)
names(ccol) <- levels(ccut)
ccol <- ccol[as.character(ccut)]
bre <- seq(min(tmp), (max(tmp)+0.001), 0.001)
library(heatmap3)
ccol <- colByValue(tmp, col = colorRampPalette(c("white", "darkred"))(length(bre)-1), breaks = bre, cex.axis = 0.8)
ActivityAMPmitoch <- as.vector(ccol)

library(gplots)
tmp <- sumDELmitoch
nCol <- 100
ccol <- colorRampPalette(c("white", "blue4"))(nCol)
ccut <- cut(tmp, nCol)
names(ccol) <- levels(ccut)
ccol <- ccol[as.character(ccut)]
bre <- seq(min(tmp), (max(tmp)+0.001), 0.001)
library(heatmap3)
ccol <- colByValue(tmp, col = colorRampPalette(c("white", "blue4"))(length(bre)-1), breaks = bre, cex.axis = 0.8)
ActivityDELmitoch <- as.vector(ccol)

library(gplots)
tmp <- sumAMPglyco
nCol <- 100
ccol <- colorRampPalette(c("white", "darkred"))(nCol)
ccut <- cut(tmp, nCol)
names(ccol) <- levels(ccut)
ccol <- ccol[as.character(ccut)]
bre <- seq(min(tmp), (max(tmp)+0.001), 0.001)
library(heatmap3)
ccol <- colByValue(tmp, col = colorRampPalette(c("white", "darkred"))(length(bre)-1), breaks = bre, cex.axis = 0.8)
ActivityAMPglyco <- as.vector(ccol)

library(gplots)
tmp <- sumDELglyco
nCol <- 100
ccol <- colorRampPalette(c("white", "blue4"))(nCol)
ccut <- cut(tmp, nCol)
names(ccol) <- levels(ccut)
ccol <- ccol[as.character(ccut)]
bre <- seq(min(tmp), (max(tmp)+0.001), 0.001)
library(heatmap3)
ccol <- colByValue(tmp, col = colorRampPalette(c("white", "blue4"))(length(bre)-1), breaks = bre, cex.axis = 0.8)
ActivityDELglyco <- as.vector(ccol)

annot2 <- cbind(ActivityDELglyco, ActivityAMPglyco, ActivityDELmitoch, ActivityAMPmitoch)
colnames(annot2) <- c("DELglyco", "AMPglyco" ,"DELmito", "AMPmito")
ColSideColors = annot2
rownames(ColSideColors) <- NULL

NES_geneticAlt <- t(toPlot)
table(as.vector(NES_geneticAlt))

library(gplots)
library(heatmap3)
heatmap3(t(NES_geneticAlt), 
         ColSideColors = ColSideColors,  
         key = FALSE, 
         Colv = NA,
         Rowv = NA,
         scale = "none", 
         col = colorRampPalette(c("white","red", "blue"))(3), 
         dendrogram = "none",
         labRow = NA, labCol = NA,
         side.height.fraction = 2, keysize = 0.8, cexRow = 0.4)
