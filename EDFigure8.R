### ED Figure 8a
load("RData/TCGA_methData.RData", verbose = T)

library(openxlsx)
classTab <- read.xlsx("tables/Supplementary Table 6.xlsx", startRow = 3, sheet = 3)

groups <- classTab$Core.classification
names(groups) <- classTab$TCGA.ID

groupsMeth <- groups[colnames(methData)]

gg <- names(table(groupsMeth))
res <- vector("list", 4)
names(res) <- gg
for(i in 1:length(gg)){
  (whichOfInterest <- gg[i])
  
  ### degs methData
  gruppo1_meth <- names(groupsMeth)[groupsMeth == whichOfInterest]
  gruppo2_meth <- names(groupsMeth)[groupsMeth != whichOfInterest]
  ans <- apply(methData, 1, function(x) wilcox.test(x[gruppo1_meth], x[gruppo2_meth])$p.value)
  fc <- rowMeans(methData[, gruppo1_meth]) - rowMeans(methData[, gruppo2_meth])
  DEGS <- cbind(fc, ans)
  DEGS <- DEGS[order(DEGS[,"fc"], decreasing = TRUE),]
  fdr <- p.adjust(DEGS[, "ans"], method = "fdr")
  DEGS <- cbind(DEGS, fdr)
  
  DEGS_methData <- as.data.frame(DEGS)
  colnames(DEGS_methData)[2:3] <- c("pValue", "qValue")
  colnames(DEGS_methData) <- paste0("md_", colnames(DEGS_methData))
  DEGS_methData <- cbind(ProbeID = rownames(DEGS_methData), DEGS_methData, stringsAsFactors = F)
  
  
  signif <- DEGS_methData[DEGS_methData$md_fc > 0.58 & DEGS_methData$md_pValue < 0.01, ]
  signif <- signif[order(signif$md_fc, decreasing = T), ]
  
  res[[i]] <- signif
}

groupsMeth <- sort(groupsMeth)

res <- lapply(res, function(x) x$ProbeID[1:100])
allProbes <- unlist(res)
toPlot <- methData[allProbes, names(groupsMeth)]

classCol <- c("green3", "red", "blue", "cyan")
names(classCol) <- c("MTC", "GPM", "NEU", "PPR")

ColSideColors <- classCol[groupsMeth]

RowSideColors <- rep("red", nrow(toPlot))
RowSideColors[rownames(toPlot) %in% res$MTC] <- "green3"
RowSideColors[rownames(toPlot) %in% res$NEU] <- "blue"
RowSideColors[rownames(toPlot) %in% res$PPR] <- "cyan"

toPlot <- (toPlot - rowMeans(toPlot))/apply(toPlot, 1, sd)

summary(as.vector(toPlot))
quantile(toPlot, 0.05); quantile(toPlot, 0.9)
toPlot[toPlot <= quantile(toPlot, 0.05)] <- quantile(toPlot, 0.05)
toPlot[toPlot >= quantile(toPlot, 0.9)] <- quantile(toPlot, 0.9)

library(heatmap3)
heatmap3(toPlot[nrow(toPlot):1, ], showRowDendro = F, showColDendro = F,
         Rowv = NA,
         Colv = NA,
         ColSideColors = ColSideColors,
         ColSideLabs = NA,
         RowSideColors = RowSideColors[nrow(toPlot):1],
         RowSideLabs = NA,
         col = colorRampPalette(c("blue", "white", "red"))(100),
         labCol = NA,
         labRow = NA,
         cexCol = 0.5,
         cexRow = 1,
         scale = 'none', useRaster = F)


### ED Figure 8b
#upload the differential expression miRNA analysis in TCGA GBM subtypes
load("RData/miRNAdiff_TCGA.RData", verbose = T)

#upload the significant differential expressed miRNA in GPM subtype
library(openxlsx)
miRNAdiffexpr <- read.xlsx("tables/Supplementary Table 14.xlsx", sheet = 3)
miRNAdiffexprGPM <- miRNAdiffexpr[3:nrow(miRNAdiffexpr), 1:4]
colnames(miRNAdiffexprGPM) <- c("miRNA ID","FC","pValue", "FDR")
miRNAdiffexprGPM <- miRNAdiffexprGPM[!is.na(miRNAdiffexprGPM$`miRNA ID`),]
rownames(miRNAdiffexprGPM) <- miRNAdiffexprGPM$`miRNA ID`

#upload the significant differential expressed miRNA in MTC subtype
miRNAdiffexpr <- read.xlsx("tables/Supplementary Table 14.xlsx", sheet = 3)
miRNAdiffexprMTC <- miRNAdiffexpr[3:nrow(miRNAdiffexpr), 5:8]
colnames(miRNAdiffexprMTC) <- c("miRNA ID","FC","pValue", "FDR")
miRNAdiffexprMTC <- miRNAdiffexprMTC[!is.na(miRNAdiffexprMTC$`miRNA ID`),]
rownames(miRNAdiffexprMTC) <- miRNAdiffexprMTC$`miRNA ID`

#upload the significant differential expressed miRNA in NEU subtype
miRNAdiffexpr <- read.xlsx("tables/Supplementary Table 14.xlsx", sheet = 3)
miRNAdiffexprNEU <- miRNAdiffexpr[3:nrow(miRNAdiffexpr), 9:12]
colnames(miRNAdiffexprNEU) <- c("miRNA ID","FC","pValue", "FDR")
miRNAdiffexprNEU <- miRNAdiffexprNEU[!is.na(miRNAdiffexprNEU$`miRNA ID`),]
rownames(miRNAdiffexprNEU) <- miRNAdiffexprNEU$`miRNA ID`

#upload the significant differential expressed miRNA in PPR subtype
miRNAdiffexpr <- read.xlsx("tables/Supplementary Table 14.xlsx", sheet = 3)
miRNAdiffexprPPR <- miRNAdiffexpr[3:nrow(miRNAdiffexpr), 13:16]
colnames(miRNAdiffexprPPR) <- c("miRNA ID","FC","pValue", "FDR")
miRNAdiffexprPPR <- miRNAdiffexprPPR[!is.na(miRNAdiffexprPPR$`miRNA ID`),]
rownames(miRNAdiffexprPPR) <- miRNAdiffexprPPR$`miRNA ID`

res <- diffmiRNA_GPMvsOthers
vect <- rep("black", nrow(res))
names(vect) <- rownames(res)
vect[rownames(miRNAdiffexprGPM)] <- "red"
plot(res[, "FC"], -log10(res[, "pValue"]), pch=1, cex = 0.8, main="Volcano plot", xlim = c(-2.5, 3), col = vect, cex.axis = 1.3)
abline(h = -log10(0.0005), v = 0, col = "darkgray", lty = 3)

res <- diffmiRNA_MTCvsOthers
vect <- rep("black", nrow(res))
names(vect) <- rownames(res)
vect[rownames(miRNAdiffexprMTC)] <- "red"
plot(res[, "FC"], -log10(res[, "pValue"]), pch=1, cex = 0.8, main="Volcano plot", xlim = c(-2.5, 3), col = vect, cex.axis = 1.3)
abline(h = -log10(0.0005), v = 0, col = "darkgray", lty = 3)

res <- diffmiRNA_NEUvsOthers
vect <- rep("black", nrow(res))
names(vect) <- rownames(res)
vect[rownames(miRNAdiffexprNEU)] <- "red"
plot(res[, "FC"], -log10(res[, "pValue"]), pch=1, cex = 0.8, main="Volcano plot", xlim = c(-2.5, 3), col = vect, cex.axis = 1.3)
abline(h = -log10(0.0005), v = 0, col = "darkgray", lty = 3)

res <- diffmiRNA_PPRvsOthers
vect <- rep("black", nrow(res))
names(vect) <- rownames(res)
vect[rownames(miRNAdiffexprPPR)] <- "red"
plot(res[, "FC"], -log10(res[, "pValue"]), pch=1, cex = 0.8, main="Volcano plot", xlim = c(-2.5, 3), col = vect, cex.axis = 1.3)
abline(h = -log10(0.0005), v = 0, col = "darkgray", lty = 3)

