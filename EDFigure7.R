### ED Figure 7b
load("RData/percTumor_percNonTumor_allDatasets.RData", verbose = T)
#Removed tumors with < 5 normal cells: MGH129_D3 = 2 cells, MGH136_D3 = 2, MGH143_D3 = 1, MGH152_D3 = 1, MGH66_D3 = 1, MGH100_D3 = 0, MGH104_D3 = 0

ccor <- matrix(0, nrow = ncol(percTumor), ncol = ncol(percNonTumor))
rownames(ccor) <- colnames(percTumor)
colnames(ccor) <- colnames(percNonTumor)
ccorPval <- ccor

for(i in 1:ncol(percTumor))
  for(j in 1:ncol(percNonTumor)){
    tmp <- cor.test(percTumor[, i], percNonTumor[, j])
    ccor[i, j] <- tmp$estimate
    ccorPval[i, j] <- tmp$p.value
  }

colnames(ccor) <- colnames(ccorPval) <- gsub(".", "-", colnames(ccor), fixed = T)

ccor
ccorPval

tumorClassCol <- c("green3", "red", "blue", "cyan")
names(tumorClassCol) <- c("Mitochondrial", "Glycolytic", "Neuronal", "Proliferative")
RowSideColors <- tumorClassCol[rownames(ccor)]

nonTumorClassCol <- c("green3", "red", "yellow4", "magenta", "navy", "blue", "thistle1", "purple")
names(nonTumorClassCol) <- c("Endothelial", "Macrophages", "Microglia", "Monocytes", "Neutrophils", "Oligodendrocytes", "Perycites", "T-lymphocytes")
ColSideColors <- nonTumorClassCol[colnames(ccor)]


library(heatmap3)
library(gplots)
heatmap3(ccor[nrow(ccor):1, ], showRowDendro = F, showColDendro = T,
         Rowv = NA,
         ColSideColors = ColSideColors,
         RowSideColors = RowSideColors[nrow(ccor):1],
         labRow = c("PPR", "NEU", "MTC", "GPM"),
         RowSideLabs = NA,
         ColSideLabs = NA,
         col = colorRampPalette(c("blue4", "white", "darkred"))(75),
         cexCol = 1.5,
         cexRow = 1.5,
         scale = 'none', margins = c(20, 20))


### ED Figure 7c
load("RData/Dataset1_geData_nonTumor.RData", verbose = T)

load("RData/phenotypesNonTumor_allDatasets.RData", verbose = T)
pData <- phenotypesNonTumor[colnames(geData), ]
pData <- pData[pData$Sample %in% c("S1", "S5", "S4", "S12"), ]
pData$Phenotype[pData$Phenotype %in% c("Macrophages", "Microglia", "Monocytes", "Neutrophils")] <- "Myeloid-cells"
pData$Phenotype[pData$Phenotype %in% c("Endothelial", "Perycites")] <- "Blood-vessel"
pData <- pData[pData$Phenotype == "Myeloid-cells", ]
geData <- geData[, rownames(pData)]

whichOfInterest <- rownames(pData)[pData$Sample %in% c("S4", "S12")]
theOthers <- setdiff(rownames(pData), whichOfInterest)

load("RData/Venteicher2017_macrophages_microglia_signatures.RData", verbose = T)
geneSet <- lapply(geneSet, function(x) intersect(x, rownames(geData)))

o <- colMeans(geData[geneSet$VenteicherMacrophages, ]) - colMeans(geData[geneSet$VenteicherMicroglia, ])
o <- sort(o, decreasing = T)

toPlot <- geData[unlist(geneSet), names(o)]

tmp <- cor(o, t(toPlot), method = "spearman")
ccor <- as.numeric(tmp)
names(ccor) <- colnames(tmp)
ccor <- sort(ccor)
ccor <- c(head(ccor, 25), tail(ccor, 25))
toPlot <- toPlot[names(ccor), ]

ColSideColors <- rep("red", ncol(toPlot))
ColSideColors[colnames(toPlot) %in% theOthers] <- "green3"

toPlot <- (toPlot - rowMeans(toPlot))/apply(toPlot, 1, sd)

summary(as.vector(toPlot))
quantile(toPlot, 0.15); quantile(toPlot, 0.85)
toPlot[toPlot <= quantile(toPlot, 0.15)] <- quantile(toPlot, 0.15)
toPlot[toPlot >= quantile(toPlot, 0.85)] <- quantile(toPlot, 0.85)


library(heatmap3)
library(gplots)
heatmap3(toPlot, showRowDendro = F, showColDendro = F,
         Rowv = NA,
         Colv = NA,
         ColSideColors = ColSideColors, ColSideLabs = NA,
         RowSideLabs = NA,
         col = colorRampPalette(c("#053061", "#4393C3", "cornsilk", "#D6604D", "#67001F"))(75),
         labCol = NA,
         cexCol = 1,#3,
         cexRow = 0.6,
         scale = 'none', useRaster = F)
