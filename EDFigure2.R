### ED Figure 2a

# upload the dataframe relative to the differential expression analysis between MTC compared to the others
# the dataframe should include a column called "FC" in which log2FC between two conditions is reported

load("RData/sc_dataset1_MTCvsOthers_DEGs.RData")
dataset1_FC <- DEGs_MTCvsOth_dataset1[, "FC"]
load("RData/sc_dataset2_MTCvsOthers_DEGs.RData")
dataset2_FC <- DEGs_MTCvsOth_dataset2[, "FC"]
load("RData/sc_dataset3_MTCvsOthers_DEGs.RData")
dataset3_FC <- DEGs_MTCvsOth_dataset3[, "FC"]

# upload the dataframe relative significantly differential expressed genes (see methods for thresholds used) between the subtype of interest compared to the others
# the rownames of dataframe should be the gene names of significantly differential expressed genes (see supplementary table 3a)
library(openxlsx)
SuppTab3 <- read.xlsx("tables/Supplementary Table 3.xlsx", sheet = 1)
sigDEGs_dataset1 <- SuppTab3[SuppTab3$X6 == "Dataset1",]
colnames(sigDEGs_dataset1) <- c("Gene", "fold change", 	"p-value", "FDR", "Dataset")
sigDEGs_dataset1 <- sigDEGs_dataset1[-1, ]
rownames(sigDEGs_dataset1) <- sigDEGs_dataset1$Gene

sigDEGs_dataset2 <-SuppTab3[SuppTab3$X6 == "Dataset 2",]
colnames(sigDEGs_dataset2) <- c("Gene", "fold change", 	"p-value", "FDR", "Dataset")
sigDEGs_dataset2 <- sigDEGs_dataset2[-1, ]
rownames(sigDEGs_dataset2) <- sigDEGs_dataset2$Gene

sigDEGs_dataset3 <-SuppTab3[SuppTab3$X6 == "Dataset 3",]
colnames(sigDEGs_dataset3) <- c("Gene", "fold change", 	"p-value", "FDR", "Dataset")
sigDEGs_dataset3 <- sigDEGs_dataset3[-1, ]
rownames(sigDEGs_dataset3) <- sigDEGs_dataset3$Gene

# intersect genes which expression has been detected in all three single cell datasets
gg <- intersect(intersect(names(dataset1_FC), names(dataset2_FC)), names(dataset3_FC))
dataset1_FC <- dataset1_FC[gg]
dataset2_FC <- dataset2_FC[gg]
dataset3_FC <- dataset3_FC[gg]

# order the vector of the FC derived for all genes in decreasing order
dataset1_FCord <- dataset1_FC[order(dataset1_FC, decreasing = F)]
dataset2_FCord <- dataset2_FC[order(dataset2_FC, decreasing = F)]
dataset3_FCord <- dataset3_FC[order(dataset3_FC, decreasing = F)]

# generate plots ED Fig. 2a
plot(x = 1:length(dataset1_FCord), y = dataset1_FCord, col="gray", pch=16, type = 'n', axes = F, ylim = c(-4, 4))
points(which(!names(dataset1_FCord)%in%rownames(sigDEGs_dataset1)), dataset1_FCord[which(!names(dataset1_FCord)%in%rownames(sigDEGs_dataset1))], pch = 16, col="gray78", cex =0.2)
points(which(!names(dataset2_FCord)%in%rownames(sigDEGs_dataset2)), dataset2_FCord[which(!names(dataset2_FCord)%in%rownames(sigDEGs_dataset2))], pch = 16, col="gray38", cex = 0.2)
points(which(!names(dataset3_FCord)%in%rownames(sigDEGs_dataset3)), dataset3_FCord[which(!names(dataset3_FCord)%in%rownames(sigDEGs_dataset3))], pch = 16, col="gray58", cex = 0.2)
abline(h = 0.3,  col = "darkgray", lty = 3)
points(which(names(dataset1_FCord)%in%rownames(sigDEGs_dataset1)), dataset1_FCord[which(names(dataset1_FCord)%in%rownames(sigDEGs_dataset1))], pch = 16, col="red", cex = 0.45)
points(which(names(dataset2_FCord)%in%rownames(sigDEGs_dataset2)), dataset2_FCord[which(names(dataset2_FCord)%in%rownames(sigDEGs_dataset2))], pch = 16, col="red", cex = 0.45)
points(which(names(dataset3_FCord)%in%rownames(sigDEGs_dataset3)), dataset3_FCord[which(names(dataset3_FCord)%in%rownames(sigDEGs_dataset3))], pch = 16, col="red", cex = 0.45)
axis(side = 2, at= -4:4, labels = -4:4)


### ED Figure 2b

# upload the dataframe relative to the differential expression analysis between NEU compared to the others
# the dataframe should include a column called "log2FC" in which log2FC between two conditions is reported
load("RData/sc_dataset1_NEUvsOthers_DEGs.RData")
dataset1_FC <- DEGs_NEUvsOth_dataset1[, "FC"]
load("RData/sc_dataset2_NEUvsOthers_DEGs.RData")
dataset2_FC <- DEGs_NEUvsOth_dataset2[, "FC"]
load("RData/sc_dataset3_NEUvsOthers_DEGs.RData")
dataset3_FC <- DEGs_NEUvsOth_dataset3[, "FC"]

# upload the dataframe relative significantly differential expressed genes (see methods for thresholds used) between the subtype of interest compared to the others
# the rownames of dataframe should be the gene names of significantly differential expressed genes (see supplementary table 3b)
library(openxlsx)
SuppTab3 <- read.xlsx("tables/Supplementary Table 3.xlsx", sheet = 2)
sigDEGs_dataset1 <- SuppTab3[SuppTab3$X6 == "Dataset 1",]
colnames(sigDEGs_dataset1) <- c("Gene", "fold change", 	"p-value", "FDR", "Dataset")
sigDEGs_dataset1 <- sigDEGs_dataset1[-1, ]
rownames(sigDEGs_dataset1) <- sigDEGs_dataset1$Gene

sigDEGs_dataset2 <-SuppTab3[SuppTab3$X6 == "Dataset 2",]
colnames(sigDEGs_dataset2) <- c("Gene", "fold change", 	"p-value", "FDR", "Dataset")
sigDEGs_dataset2 <- sigDEGs_dataset2[-1, ]
rownames(sigDEGs_dataset2) <- sigDEGs_dataset2$Gene

sigDEGs_dataset3 <-SuppTab3[SuppTab3$X6 == "Dataset 3",]
colnames(sigDEGs_dataset3) <- c("Gene", "fold change", 	"p-value", "FDR", "Dataset")
sigDEGs_dataset3 <- sigDEGs_dataset3[-1, ]
rownames(sigDEGs_dataset3) <- sigDEGs_dataset3$Gene

# intersect genes which expression has been detected in all three single cell datasets
gg <- intersect(intersect(names(dataset1_FC), names(dataset2_FC)), names(dataset3_FC))
dataset1_FC <- dataset1_FC[gg]
dataset2_FC <- dataset2_FC[gg]
dataset3_FC <- dataset3_FC[gg]

# order the vector of the FC derived for all genes in decreasing order
dataset1_FCord <- dataset1_FC[order(dataset1_FC, decreasing = F)]
dataset2_FCord <- dataset2_FC[order(dataset2_FC, decreasing = F)]
dataset3_FCord <- dataset3_FC[order(dataset3_FC, decreasing = F)]

# generate plots ED Fig. 2b

plot(x = 1:length(dataset1_FCord), y = dataset1_FCord, col="gray", pch=16, type = 'n', axes = F, ylim = c(-4, 4))
points(which(!names(dataset1_FCord)%in%rownames(sigDEGs_dataset1)), dataset1_FCord[which(!names(dataset1_FCord)%in%rownames(sigDEGs_dataset1))], pch = 16, col="gray78", cex =0.2)
points(which(!names(dataset2_FCord)%in%rownames(sigDEGs_dataset2)), dataset2_FCord[which(!names(dataset2_FCord)%in%rownames(sigDEGs_dataset2))], pch = 16, col="gray38", cex = 0.2)
points(which(!names(dataset3_FCord)%in%rownames(sigDEGs_dataset3)), dataset3_FCord[which(!names(dataset3_FCord)%in%rownames(sigDEGs_dataset3))], pch = 16, col="gray58", cex = 0.2)
abline(h = 0.3,  col = "darkgray", lty = 3)
points(which(names(dataset1_FCord)%in%rownames(sigDEGs_dataset1)), dataset1_FCord[which(names(dataset1_FCord)%in%rownames(sigDEGs_dataset1))], pch = 16, col="red", cex = 0.45)
points(which(names(dataset2_FCord)%in%rownames(sigDEGs_dataset2)), dataset2_FCord[which(names(dataset2_FCord)%in%rownames(sigDEGs_dataset2))], pch = 16, col="red", cex = 0.45)
points(which(names(dataset3_FCord)%in%rownames(sigDEGs_dataset3)), dataset3_FCord[which(names(dataset3_FCord)%in%rownames(sigDEGs_dataset3))], pch = 16, col="red", cex = 0.45)
axis(side = 2, at= -4:4, labels = -4:4)


### ED Figure 2c

# upload the dataframe relative to the differential expression analysis between PPR compared to the others
# the dataframe should include a column called "log2FC" in which log2FC between two conditions is reported

load("RData/sc_dataset1_PPRvsOthers_DEGs.RData")
dataset1_FC <- DEGs_PPRvsOth_dataset1[, "FC"]
load("RData/sc_dataset2_PPRvsOthers_DEGs.RData")
dataset2_FC <- DEGs_PPRvsOth_dataset2[, "FC"]
load("RData/sc_dataset3_PPRvsOthers_DEGs.RData")
dataset3_FC <- DEGs_PPRvsOth_dataset3[, "FC"]

# upload the dataframe relative significantly differential expressed genes (see methods for thresholds used) between the subtype of interest compared to the others
# the rownames of dataframe should be the gene names of significantly differential expressed genes (see supplementary table 3c-d)
library(openxlsx)
SuppTab3 <- read.xlsx("tables/Supplementary Table 3.xlsx", sheet = 4)
sigDEGs_dataset1 <- SuppTab3[SuppTab3$X5 == "Dataset 1",]
colnames(sigDEGs_dataset1) <- c("Gene", "fold change", 	"p-value", "FDR", "Dataset")
sigDEGs_dataset1 <- sigDEGs_dataset1[-1, ]
rownames(sigDEGs_dataset1) <- sigDEGs_dataset1$Gene

sigDEGs_dataset2 <-SuppTab3[SuppTab3$X5 == "Dataset 2",]
colnames(sigDEGs_dataset2) <- c("Gene", "fold change", 	"p-value", "FDR", "Dataset")
sigDEGs_dataset2 <- sigDEGs_dataset2[-1, ]
rownames(sigDEGs_dataset2) <- sigDEGs_dataset2$Gene

sigDEGs_dataset3 <-SuppTab3[SuppTab3$X5 == "Dataset 3",]
colnames(sigDEGs_dataset3) <- c("Gene", "fold change", 	"p-value", "FDR", "Dataset")
sigDEGs_dataset3 <- sigDEGs_dataset3[-1, ]
rownames(sigDEGs_dataset3) <- sigDEGs_dataset3$Gene

# intersect genes which expression has been detected in all three single cell datasets
gg <- intersect(intersect(names(dataset1_FC), names(dataset2_FC)), names(dataset3_FC))
dataset1_FC <- dataset1_FC[gg]
dataset2_FC <- dataset2_FC[gg]
dataset3_FC <- dataset3_FC[gg]

# order the vector of the FC derived for all genes in decreasing order
dataset1_FCord <- dataset1_FC[order(dataset1_FC, decreasing = F)]
dataset2_FCord <- dataset2_FC[order(dataset2_FC, decreasing = F)]
dataset3_FCord <- dataset3_FC[order(dataset3_FC, decreasing = F)]

# generate plots ED Fig. 2c
plot(x = 1:length(dataset1_FCord), y = dataset1_FCord, col="gray", pch=16, type = 'n', axes = F, ylim = c(-4, 4))
points(which(!names(dataset1_FCord)%in%rownames(sigDEGs_dataset1)), dataset1_FCord[which(!names(dataset1_FCord)%in%rownames(sigDEGs_dataset1))], pch = 16, col="gray78", cex =0.2)
points(which(!names(dataset2_FCord)%in%rownames(sigDEGs_dataset2)), dataset2_FCord[which(!names(dataset2_FCord)%in%rownames(sigDEGs_dataset2))], pch = 16, col="gray38", cex = 0.2)
points(which(!names(dataset3_FCord)%in%rownames(sigDEGs_dataset3)), dataset3_FCord[which(!names(dataset3_FCord)%in%rownames(sigDEGs_dataset3))], pch = 16, col="gray58", cex = 0.2)
abline(h = 0.3,  col = "darkgray", lty = 3)
points(which(names(dataset1_FCord)%in%rownames(sigDEGs_dataset1)), dataset1_FCord[which(names(dataset1_FCord)%in%rownames(sigDEGs_dataset1))], pch = 16, col="red", cex = 0.45)
points(which(names(dataset2_FCord)%in%rownames(sigDEGs_dataset2)), dataset2_FCord[which(names(dataset2_FCord)%in%rownames(sigDEGs_dataset2))], pch = 16, col="red", cex = 0.45)
points(which(names(dataset3_FCord)%in%rownames(sigDEGs_dataset3)), dataset3_FCord[which(names(dataset3_FCord)%in%rownames(sigDEGs_dataset3))], pch = 16, col="red", cex = 0.45)
axis(side = 2, at= -4:4, labels = -4:4)


### ED Figure 2d

library(openxlsx)
classTab <- read.xlsx("tables/Supplementary Table 4.xlsx", startRow = 4)

class <- classTab$`Pathway-based`
neftelClass <- classTab$Neftel.et.al.

table(class, neftelClass)
chisq.test(table(class, neftelClass))

transMat <- as.data.frame(table(neftelClass, class), stringsAsFactors = F)
transMat$class <- factor(transMat$class, levels = c("GPM", "MTC", "NEU", "PPR"))
transMat$neftelClass <- factor(transMat$neftelClass, levels = c("AC", "MES1", "MES2", "NPC1", "NPC2", "OPC"))

library(flipPlots)
SankeyDiagram(transMat[, -3],
              link.color = "Source", 
              label.show.varname = FALSE,
              font.size = 0,
              weights = transMat$Freq,
              node.width = 100,
              colors = c("aquamarine3", "chocolate1", "slateblue1", "violet", "chartreuse3", "yellow2", "red", "green3", "blue", "cyan")
) 

table(class)/length(class)
table(neftelClass)/length(neftelClass)


###ED Figure 2e
load("RData/Dataset1_class.RData", verbose = T)
load("RData/Dataset1_pData.RData", verbose = T)
pData <- cbind(pData, class = class, stringsAsFactors = F)
dataset1 <- pData
dataset1$TumorID <- paste0(dataset1$TumorID, "_D1")

load("RData/Dataset2_class.RData", verbose = T)
load("RData/Dataset2_pData.RData", verbose = T)
pData <- cbind(pData, class = class, stringsAsFactors = F)
dataset2 <- pData
dataset2$TumorID <- paste0(dataset2$TumorID, "_D2")

load("RData/Dataset3_class.RData", verbose = T)
load("RData/Dataset3_pData.RData", verbose = T)
pData <- cbind(pData, class = class, stringsAsFactors = F)
dataset3 <- pData
dataset3$TumorID <- paste0(dataset3$TumorID, "_D3")

pData <- data.frame(TumorID = c(dataset1$TumorID, dataset2$TumorID, dataset3$TumorID),
                    class = c(dataset1$class, dataset2$class, dataset3$class), stringsAsFactors = F)
ccol <- c("red", "green3", "blue", "cyan")
names(ccol) <- c("Glycolytic", "Mitochondrial", "Neuronal", "Proliferative")

percTumor <- table(pData$TumorID, pData$class)
percTumor <- as.matrix(percTumor)

tmp <- rowSums(percTumor != 0)
(tmp <- table(tmp))
toAdd <- c(0, 0)
names(toAdd) <- c("1", "2")
tmp <- c(toAdd, tmp)

barplot(tmp, col = "gray80", xlab = "Number of GBM subtypes", ylab = "Number of tumors", ylim = c(0, 30), )


###ED Figure 2g-h
#see Figure 2d, e for code