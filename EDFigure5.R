### ED Figure 5a

#upload 192 pathways associated with GBM IDHwt patient survival
library(openxlsx)
survAssPath192 <- read.xlsx("tables/Supplementary Table 6.xlsx", sheet = 2)
survAssPath192 <- survAssPath192[2:nrow(survAssPath192), ]
colnames(survAssPath192) <- survAssPath192[1, ]
survAssPath192 <- survAssPath192[-1, ]
rownames(survAssPath192) <- survAssPath192$Pathway

#upload the NES (normalized enriched scores) activity of 5,032 biological pathways derived from all GBM IDHwt tumors using MWWGST
load("RData/TCGA_NESActivity_HallmarksGOKegg.RData", verbose = T)
ss <- intersect(rownames(survAssPath192), rownames(NES))

NES2 <- NES[ss,]
FDR2 <- FDR[ss,]

#Consensus clustering using 192 survival associated pathways
tmp <- matrix(0, nrow=ncol(NES2), ncol=ncol(NES2))
rownames(tmp) <- colnames(NES2)
colnames(tmp) <- colnames(NES2)

for (i in 1:nrow(NES2)){
  g1 <- names(which(NES2[i,] > 0 & FDR2[i,] < 0.01))
  g2 <- names(which(NES2[i,] < 0 & FDR2[i,] < 0.01))
  if (length(g1) == 0 | length(g2) == 0) next
  g3 <- setdiff(names(NES2[i,]), c(g1, g2))
  tmp[g1, g1] <- tmp[g1, g1]+1
  tmp[g2, g2] <- tmp[g2, g2]+1
  tmp[g3, g3] <- tmp[g3, g3]+1
}

tmp2 <- tmp/nrow(NES2)
consensusmatrix <- tmp2
tmp3 <- 1-tmp2
consensusdist <- tmp3

library(parallel)
library(doParallel)

source("functions.R")

ddist <- as.dist(consensusdist)

aMakeCluster <- makeCluster(4)
aConsensus <- consensusClust(ddist,
                             G = 5,
                             runs = 10000,
                             epsilon = .7,
                             clustMethod = "ward",
                             method = "monti",
                             bySample = TRUE,
                             cpuCluster = aMakeCluster)
stopCluster(aMakeCluster)

# generating ED Fig. 5a
obj <- aConsensus
G <- obj$G
consHc <- hclust(obj$consensus.diss, method = "complete")
consClust <- as.factor(cutree(consHc, k = G))
levels(consClust) <- palette()[1:G]
table(consClust)

tmp <- paste("g", c(1:G), sep = "")
names(tmp) <- palette()[1:G]
tmp <- tmp[consClust]
names(tmp) <- names(consClust)
consClust <- tmp
table(consClust)

col <- ifelse(obj$bySample, "darkred", "royalblue")
col <- colorRampPalette(c("white", col))(50)

RowSideColors <- as.character(consClust)[consHc$order]
RowSideColors <- cbind(RowSideColors)
colnames(RowSideColors) <- ""

colClusters <- as.character(consClust)
tmp <- palette()[1:G]
names(tmp) <- paste("g", c(1:G), sep = "")
colClusters <- as.character(tmp[consClust]) 
table(colClusters)
annot2 <- cbind(colClusters)
rownames(annot2) <- names(consClust)
ColSideColors = annot2

Colv  <- as.dendrogram(consHc)
consMatrix <- as.matrix(1 - obj$consensus.diss)
diag(consMatrix) <- 1
consMatrix <- consMatrix[, consHc$order]
library(heatmap3)
heatmap3(t(consMatrix), col = col, scale = "none", labCol = NA,
         labRow = NA, Colv = rev(Colv), Rowv = NA,
         ColSideColors = ColSideColors, 
         cexCol = 1.5, keep.dendro = T)


### ED Figure 5b

# upload the differential pathways among the 192 survival associated pathways in TCGA
library(openxlsx)
diff_survAssPath192 <- read.xlsx("tables/Supplementary Table 6.xlsx", sheet = 4)
diff_survAssPath192 <- diff_survAssPath192[2:nrow(diff_survAssPath192), ]
colnames(diff_survAssPath192) <- diff_survAssPath192[1, ]
diff_survAssPath192 <- diff_survAssPath192[-1, ]
rownames(diff_survAssPath192) <- diff_survAssPath192$Pathway

# order the differential pathways among the 192 survival associated pathways in each subtype by effect size
diff_survAssPath192_GPM <- diff_survAssPath192[diff_survAssPath192$`GBM subtype` == "GPM",]
diff_survAssPath192_MTC <- diff_survAssPath192[diff_survAssPath192$`GBM subtype` == "MTC",]
diff_survAssPath192_PPR <- diff_survAssPath192[diff_survAssPath192$`GBM subtype` == "PPR",]
diff_survAssPath192_NEU <- diff_survAssPath192[diff_survAssPath192$`GBM subtype` == "NEU",]
diff_survAssPath192_GPM <- diff_survAssPath192_GPM[order(diff_survAssPath192_GPM$`effect size subtype vs others`, decreasing = F),]
diff_survAssPath192_MTC <- diff_survAssPath192_MTC[order(diff_survAssPath192_MTC$`effect size subtype vs others`, decreasing = F),]
diff_survAssPath192_PPR <- diff_survAssPath192_PPR[order(diff_survAssPath192_PPR$`effect size subtype vs others`, decreasing = F),]
diff_survAssPath192_NEU <- diff_survAssPath192_NEU[order(diff_survAssPath192_NEU$`effect size subtype vs others`, decreasing = F),]

# create a list with the differential pathways among the 192 survival associated pathways in each subtype and color annotation
listOfDiffs <- vector("list", 4)
names(listOfDiffs) <- c("GPM", "MTC", "PPR", "NEU")
listOfDiffs$GPM <- rownames(diff_survAssPath192_GPM)
listOfDiffs$MTC <-  rownames(diff_survAssPath192_MTC)
listOfDiffs$PPR <-  rownames(diff_survAssPath192_PPR)
listOfDiffs$NEU <-  rownames(diff_survAssPath192_NEU)

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

colPath <- c(aa, bb, cc, dd)
table(colPath)
annot3 <- cbind(colPath)
colnames(annot3) <- c("colPATH")
rownames(annot3) <- as.character(unlist(listOfDiffs))
RowSideColors = annot3

# upload the TCGA core classification by consensus clustering
library(openxlsx)
TCGAclassification <- read.xlsx("tables/Supplementary Table 6.xlsx", sheet = 3)
TCGAclassification <- TCGAclassification[2:nrow(TCGAclassification), ]
colnames(TCGAclassification) <- TCGAclassification[1, ]
TCGAclassification <- TCGAclassification[-1, ]
rownames(TCGAclassification) <- TCGAclassification$`TCGA ID`
TCGAclassification2 <- TCGAclassification$`Core classification`
names(TCGAclassification2) <- rownames(TCGAclassification)
table(TCGAclassification2)
TCGAclassification2 <- c(TCGAclassification2[TCGAclassification2 == "NEU"], TCGAclassification2[TCGAclassification2 == "PPR"], TCGAclassification2[TCGAclassification2 == "MTC"], TCGAclassification2[TCGAclassification2 == "GPM"])

# upload the TCGA core classification color annotation
colClass <- c("red", "green3", "cyan", "blue")
names(colClass) <- c("GPM", "MTC", "PPR", "NEU")
colClass <- colClass[TCGAclassification2]
names(colClass) <- names(TCGAclassification2)
ColSideColors = colClass
colnames(ColSideColors) <- NULL

# plot the ED. Fig. 5b
load("RData/TCGA_NESActivity_HallmarksGOKegg.RData", verbose = T)
toPlot <- NES[as.character(unlist(listOfDiffs)), names(TCGAclassification2)]

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

PathActivity <- t(toPlot)

library(heatmap3)
heatmap3(t(PathActivity), showRowDendro = F, showColDendro = F,
         Rowv = NA,
         Colv = NA,
         ColSideColors = ColSideColors, ColSideLabs = NA,
         RowSideColors = RowSideColors, RowSideLabs = NA,
         col = colorRampPalette(c("#053061", "#4393C3", "cornsilk", "#D6604D", "#67001F"))(75),
         labCol = NA,
         cexCol = 0.5,
         cexRow = 0.3,
         scale = 'none', useRaster = F)


### ED Figure 5c
# upload the differential genes in TCGA subtypes
library(openxlsx)
diff_genesGPM <- read.xlsx("tables/Supplementary Table 6.xlsx", sheet = 6)
diff_genesGPM <- diff_genesGPM[2:nrow(diff_genesGPM), ]
colnames(diff_genesGPM) <- diff_genesGPM[1, ]
diff_genesGPM <- diff_genesGPM[-1, ]
rownames(diff_genesGPM) <- diff_genesGPM$Gene

diff_genesMTC <- read.xlsx("tables/Supplementary Table 6.xlsx", sheet = 7)
diff_genesMTC <- diff_genesMTC[2:nrow(diff_genesMTC), ]
colnames(diff_genesMTC) <- diff_genesMTC[1, ]
diff_genesMTC <- diff_genesMTC[-1, ]
rownames(diff_genesMTC) <- diff_genesMTC$Gene

diff_genesPPR <- read.xlsx("tables/Supplementary Table 6.xlsx", sheet = 9)
diff_genesPPR <- diff_genesPPR[2:nrow(diff_genesPPR), ]
colnames(diff_genesPPR) <- diff_genesPPR[1, ]
diff_genesPPR <- diff_genesPPR[-1, ]
rownames(diff_genesPPR) <- diff_genesPPR$Gene

diff_genesNEU <- read.xlsx("tables/Supplementary Table 6.xlsx", sheet = 8)
diff_genesNEU <- diff_genesNEU[2:nrow(diff_genesNEU), ]
colnames(diff_genesNEU) <- diff_genesNEU[1, ]
diff_genesNEU <- diff_genesNEU[-1, ]
rownames(diff_genesNEU) <- diff_genesNEU$Gene

# upload the TCGA core classification by consensus clustering
library(openxlsx)
TCGAclassification <- read.xlsx("tables/Supplementary Table 6.xlsx", sheet = 3)
TCGAclassification <- TCGAclassification[2:nrow(TCGAclassification), ]
colnames(TCGAclassification) <- TCGAclassification[1, ]
TCGAclassification <- TCGAclassification[-1, ]
rownames(TCGAclassification) <- TCGAclassification$`TCGA ID`
TCGAclassification2 <- TCGAclassification$`Core classification`
names(TCGAclassification2) <- rownames(TCGAclassification)
table(TCGAclassification2)
TCGAclassification2 <- c(TCGAclassification2[TCGAclassification2 == "NEU"], TCGAclassification2[TCGAclassification2 == "PPR"], TCGAclassification2[TCGAclassification2 == "MTC"], TCGAclassification2[TCGAclassification2 == "GPM"])

# upload the TCGA core classification color annotation
colClass <- c("red", "green3", "cyan", "blue")
names(colClass) <- c("GPM", "MTC", "PPR", "NEU")
colClass <- colClass[TCGAclassification2]
names(colClass) <- names(TCGAclassification2)
ColSideColors = colClass
colnames(ColSideColors) <- NULL

load("RData/TCGA_geData.RData", verbose = T)

# create a list with the differential genes in each subtype and color annotation
listOfDiffs <- vector("list", 4)
names(listOfDiffs) <- c("GPM", "MTC", "PPR", "NEU")
listOfDiffs$GPM <- rownames(diff_genesGPM)
listOfDiffs$MTC <-  rownames(diff_genesMTC)
listOfDiffs$PPR <-  rownames(diff_genesPPR)
listOfDiffs$NEU <-  rownames(diff_genesNEU)

gg <- intersect(as.character(unlist(listOfDiffs)), rownames(geData))
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

geData <- geData[gg, names(TCGAclassification2)]
toPlot <- geData 

#plot the ED. Fig. 5c
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


### ED Figure 5d

load("RData/sc_dataset1_NEUvsOthers_DEGs.RData")
dataset1_FC <- DEGs_NEUvsOth_dataset1[, "FC"]
load("RData/sc_dataset2_NEUvsOthers_DEGs.RData")
dataset2_FC <- DEGs_NEUvsOth_dataset2[, "FC"]
load("RData/sc_dataset3_NEUvsOthers_DEGs.RData")
dataset3_FC <- DEGs_NEUvsOth_dataset3[, "FC"]
load("RData/DEGs_NEUvsOth_TCGA.RData")
dataset_TCGA <- DEGs_NEUvsOth_TCGA[, "FC"]

# upload the dataframe relative significantly differential expressed genes (see methods for thresholds used) between the subtype of interest compared to the others
# the rownames of dataframe should be the gene names of significantly differential expressed genes (see supplementary table 8a)
library(openxlsx)
SuppTab3 <- read.xlsx("tables/Supplementary Table 8.xlsx", sheet = 1)

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

# generate plots ED. Fig 5d
plot(x = 1:length(dataset1_FCord), y = dataset1_FCord, col="gray", pch=16, type = 'n', axes = F, ylim = c(-4, 4))
points(which(!names(dataset1_FCord)%in%rownames(sigDEGs_dataset1)), dataset1_FCord[which(!names(dataset1_FCord)%in%rownames(sigDEGs_dataset1))], pch = 16, col="gray78", cex =0.2)
points(which(!names(dataset2_FCord)%in%rownames(sigDEGs_dataset2)), dataset2_FCord[which(!names(dataset2_FCord)%in%rownames(sigDEGs_dataset2))], pch = 16, col="gray38", cex = 0.2)
points(which(!names(dataset3_FCord)%in%rownames(sigDEGs_dataset3)), dataset3_FCord[which(!names(dataset3_FCord)%in%rownames(sigDEGs_dataset3))], pch = 16, col="gray58", cex = 0.2)
points(which(!names(dataset_TCGA)%in%rownames(sigDEGs_datasetTCGA)), dataset_TCGA[which(!names(dataset_TCGA)%in%rownames(sigDEGs_datasetTCGA))], pch = 16, col="black", cex = 0.2)
abline(h = 0.3,  col = "darkgray", lty = 3)
points(which(names(dataset1_FCord)%in%rownames(sigDEGs_dataset1)), dataset1_FCord[which(names(dataset1_FCord)%in%rownames(sigDEGs_dataset1))], pch = 16, col="red", cex = 0.45)
points(which(names(dataset2_FCord)%in%rownames(sigDEGs_dataset2)), dataset2_FCord[which(names(dataset2_FCord)%in%rownames(sigDEGs_dataset2))], pch = 16, col="red", cex = 0.45)
points(which(names(dataset3_FCord)%in%rownames(sigDEGs_dataset3)), dataset3_FCord[which(names(dataset3_FCord)%in%rownames(sigDEGs_dataset3))], pch = 16, col="red", cex = 0.45)
points(which(names(dataset_TCGA)%in%rownames(sigDEGs_datasetTCGA)), dataset_TCGA[which(names(dataset_TCGA)%in%rownames(sigDEGs_datasetTCGA))], pch = 16, col="red", cex = 0.45)
axis(side = 2, at= -4:4, labels = -4:4)


### ED Figure 5e

load("RData/sc_dataset1_PPRvsOthers_DEGs.RData")
dataset1_FC <- DEGs_PPRvsOth_dataset1[, "FC"]
load("RData/sc_dataset2_PPRvsOthers_DEGs.RData")
dataset2_FC <- DEGs_PPRvsOth_dataset2[, "FC"]
load("RData/sc_dataset3_PPRvsOthers_DEGs.RData")
dataset3_FC <- DEGs_PPRvsOth_dataset3[, "FC"]
load("RData/DEGs_PPRvsOth_TCGA.RData")
dataset_TCGA <- DEGs_PPRvsOth_TCGA[, "FC"]

# upload the dataframe relative significantly differential expressed genes (see methods for thresholds used) between the subtype of interest compared to the others
# the rownames of dataframe should be the gene names of significantly differential expressed genes (see supplementary table 8a)
library(openxlsx)
SuppTab3 <- read.xlsx("tables/Supplementary Table 8.xlsx", sheet = 2)

sigDEGs_datasetTCGA <- SuppTab3[SuppTab3$X5 == "TCGA",]
colnames(sigDEGs_datasetTCGA) <- c("Gene",  "fold change", 	"p-value", "FDR", "Dataset")
sigDEGs_datasetTCGA <- sigDEGs_datasetTCGA[-1, ]
rownames(sigDEGs_datasetTCGA) <- sigDEGs_datasetTCGA$Gene

sigDEGs_dataset1 <- SuppTab3[SuppTab3$X5 == "Dataset 1",]
colnames(sigDEGs_dataset1) <- c("Gene", "fold change", 	"p-value", "FDR", "Dataset")
sigDEGs_dataset1 <- sigDEGs_dataset1[-1, ]
rownames(sigDEGs_dataset1) <- sigDEGs_dataset1$Gene

sigDEGs_dataset2 <-SuppTab3[SuppTab3$X5 == "Dataset 2",]
colnames(sigDEGs_dataset2) <- c("Gene", "fold change", 	"p-value", "FDR", "Dataset")
sigDEGs_dataset2 <- sigDEGs_dataset2[-1, ]
rownames(sigDEGs_dataset2) <- sigDEGs_dataset2$Gene

sigDEGs_dataset3 <-SuppTab3[SuppTab3$X5 == "Dataset 3",]
colnames(sigDEGs_dataset3) <- c("Gene",  "fold change", 	"p-value", "FDR", "Dataset")
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

# generate plots ED Fig 5e

plot(x = 1:length(dataset1_FCord), y = dataset1_FCord, col="gray", pch=16, type = 'n', axes = F, ylim = c(-5, 6))
points(which(!names(dataset1_FCord)%in%rownames(sigDEGs_dataset1)), dataset1_FCord[which(!names(dataset1_FCord)%in%rownames(sigDEGs_dataset1))], pch = 16, col="gray78", cex =0.2)
points(which(!names(dataset2_FCord)%in%rownames(sigDEGs_dataset2)), dataset2_FCord[which(!names(dataset2_FCord)%in%rownames(sigDEGs_dataset2))], pch = 16, col="gray38", cex = 0.2)
points(which(!names(dataset3_FCord)%in%rownames(sigDEGs_dataset3)), dataset3_FCord[which(!names(dataset3_FCord)%in%rownames(sigDEGs_dataset3))], pch = 16, col="gray58", cex = 0.2)
points(which(!names(dataset_TCGA)%in%rownames(sigDEGs_datasetTCGA)), dataset_TCGA[which(!names(dataset_TCGA)%in%rownames(sigDEGs_datasetTCGA))], pch = 16, col="black", cex = 0.2)
abline(h = 0.3,  col = "darkgray", lty = 3)
points(which(names(dataset1_FCord)%in%rownames(sigDEGs_dataset1)), dataset1_FCord[which(names(dataset1_FCord)%in%rownames(sigDEGs_dataset1))], pch = 16, col="red", cex = 0.45)
points(which(names(dataset2_FCord)%in%rownames(sigDEGs_dataset2)), dataset2_FCord[which(names(dataset2_FCord)%in%rownames(sigDEGs_dataset2))], pch = 16, col="red", cex = 0.45)
points(which(names(dataset3_FCord)%in%rownames(sigDEGs_dataset3)), dataset3_FCord[which(names(dataset3_FCord)%in%rownames(sigDEGs_dataset3))], pch = 16, col="red", cex = 0.45)
points(which(names(dataset_TCGA)%in%rownames(sigDEGs_datasetTCGA)), dataset_TCGA[which(names(dataset_TCGA)%in%rownames(sigDEGs_datasetTCGA))], pch = 16, col="red", cex = 0.45)
axis(side = 2, at= -5:5, labels = -5:5)


### ED Figure 5f
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

# upload the TCGA core classification by consensus clustering
library(openxlsx)
TCGAclassification <- read.xlsx("tables/Supplementary Table 6.xlsx", sheet = 3)
TCGAclassification <- TCGAclassification[2:nrow(TCGAclassification), ]
colnames(TCGAclassification) <- TCGAclassification[1, ]
TCGAclassification <- TCGAclassification[-1, ]
rownames(TCGAclassification) <- TCGAclassification$`TCGA ID`
TCGAclassification2 <- TCGAclassification$`Core classification`
names(TCGAclassification2) <- rownames(TCGAclassification)
table(TCGAclassification2)
TCGAclassification2 <- c(TCGAclassification2[TCGAclassification2 == "NEU"], TCGAclassification2[TCGAclassification2 == "PPR"], TCGAclassification2[TCGAclassification2 == "MTC"], TCGAclassification2[TCGAclassification2 == "GPM"])

# upload the TCGA core classification color annotation
colClass <- c("red", "green3", "cyan", "blue")
names(colClass) <- c("GPM", "MTC", "PPR", "NEU")
colClass <- colClass[TCGAclassification2]
names(colClass) <- names(TCGAclassification2)
ColSideColors = colClass
colnames(ColSideColors) <- NULL

load("RData/TCGA_geData.RData", verbose = T)

# create a list with the gene signatures in each subtype and color annotation
listOfDiffs <- vector("list", 4)
names(listOfDiffs) <- c("GPM", "MTC", "PPR", "NEU")
listOfDiffs$GPM <- geneSigGPM
listOfDiffs$MTC <-  geneSigMTC
listOfDiffs$PPR <-  geneSigPPR
listOfDiffs$NEU <-  geneSigNEU

gg <- intersect(as.character(unlist(listOfDiffs)), rownames(geData))
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

geData <- geData[gg, names(TCGAclassification2)]
toPlot <- geData 

#plot the ED. Fig. 5f
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


### ED Figure 5g
load("RData/TCGA_NES_class.RData", verbose = T)

library(openxlsx)
tmp <- read.xlsx("tables/Supplementary Table 12.xlsx", startRow = 4, sheet = 1)
class <- tmp$Core.classification
names(class) <- tmp$TCGA.ID
class <- class[class != "Undefined"]

NES <- NES[, names(class)]

deltaGM <- NES["GPM", ] - NES["MTC", ]
deltaNP <- NES["NEU", ] - NES["PPR", ]

deltaPM <- NES["PPR", ] - NES["MTC", ]
deltaNG <- NES["NEU", ] - NES["GPM", ]

xx <- yy <- rep(NA, length(class))
names(xx) <- names(yy) <- names(class)
#MTC
idx <- names(class[class == "MTC"])
xx[idx] <- deltaGM[idx]
yy[idx] <- deltaPM[idx]
#GPM
idx <- names(class[class == "GPM"])
xx[idx] <- deltaGM[idx]
yy[idx] <- deltaNG[idx]
#NEU
idx <- names(class[class == "NEU"])
xx[idx] <- deltaNP[idx]
yy[idx] <- deltaNG[idx]
#PPR
idx <- names(class[class == "PPR"])
xx[idx] <- deltaNP[idx]
yy[idx] <- deltaPM[idx]

(xRange <- range(xx))
(yRange <- range(yy))

semplicityScore <- rep(NA, length(class))
names(semplicityScore) <- names(class)
#MTC
whichOfInterest <- "MTC"  
idx <- names(class[class == whichOfInterest])
semplicityScore[idx] <- NES[whichOfInterest, idx] - colMeans(NES[rownames(NES) != whichOfInterest, idx])
#GPM
whichOfInterest <- "GPM"  
idx <- names(class[class == whichOfInterest])
semplicityScore[idx] <- NES[whichOfInterest, idx] - colMeans(NES[rownames(NES) != whichOfInterest, idx])
#NEU
whichOfInterest <- "NEU"  
idx <- names(class[class == whichOfInterest])
semplicityScore[idx] <- NES[whichOfInterest, idx] - colMeans(NES[rownames(NES) != whichOfInterest, idx])
#PPR
whichOfInterest <- "PPR"  
idx <- names(class[class == whichOfInterest])
semplicityScore[idx] <- NES[whichOfInterest, idx] - colMeans(NES[rownames(NES) != whichOfInterest, idx])

toPlot <- semplicityScore
toPlot <- (toPlot - mean(toPlot))/sd(toPlot)

bre <- seq(min(toPlot), (max(toPlot)+0.001), 0.001)
library(gplots)
library(heatmap3)
bre <- colByValue(toPlot, col = colorRampPalette(c("gray80", "gold", "red"))(length(bre)-1), breaks = bre, cex.axis = 0.8)
ccol <- as.character(bre)
names(ccol) <- rownames(bre)

plot(xx, yy, pch = 16, col = ccol[names(xx)], xlab = "", ylab = "", xlim = c(-8, 12), ylim = c(-9, 10), cex = 2, cex.axis = 1.5)
abline(h = 0, lty = 2, col = "gray")
abline(v = 0, lty = 2, col = "gray")

arrows(x0 = 0, y0 = -8, x1 = xRange[2]-0.7, y1 = -8, length = 0.1, lwd = 2)
arrows(x0 = 0, y0 = -8, x1 = xRange[1]+1, y1 = -8, length = 0.1, lwd = 2)

arrows(x0 = 0, y0 = 9, x1 = xRange[2]-0.7, y1 = 9, length = 0.1, lwd = 2)
arrows(x0 = 0, y0 = 9, x1 = xRange[1]+1, y1 = 9, length = 0.1, lwd = 2)

arrows(x0 = -7, y0 = 0, x1 = -7, y1 = yRange[2], length = 0.1, lwd = 2)
arrows(x0 = -7, y0 = 0, x1 = -7, y1 = yRange[1], length = 0.1, lwd = 2)

arrows(x0 = 11, y0 = 0, x1 = 11, y1 = yRange[2], length = 0.1, lwd = 2)
arrows(x0 = 11, y0 = 0, x1 = 11, y1 = yRange[1], length = 0.1, lwd = 2)
