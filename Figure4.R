### Figure 4a
library(openxlsx)
GPMmut <- read.xlsx("tables/Supplementary Table 12.xlsx", startRow = 4, sheet = 2)[, 1]
GPMamp <- read.xlsx("tables/Supplementary Table 12.xlsx", startRow = 4, sheet = 3)
GPMamp <- GPMamp[GPMamp$`p-value` < 0.01 & GPMamp$`p-value.GPM` < 0.01, 1]
GPMdel <- read.xlsx("tables/Supplementary Table 12.xlsx", startRow = 4, sheet = 4)
GPMdel <- GPMdel[GPMdel$`p-value` < 0.01 & GPMdel$`p-value.GPM` < 0.01, 1]
GPMmutamp <- read.xlsx("tables/Supplementary Table 12.xlsx", startRow = 4, sheet = 5)[, 1]
GPMmutdel <- read.xlsx("tables/Supplementary Table 12.xlsx", startRow = 4, sheet = 6)[, 1]

MTCmut <- read.xlsx("tables/Supplementary Table 12.xlsx", startRow = 4, sheet = 7)[, 1]
MTCamp <- read.xlsx("tables/Supplementary Table 12.xlsx", startRow = 4, sheet = 8)
MTCamp <- MTCamp[MTCamp$`p-value` < 0.01 & MTCamp$`p-value.MTC` < 0.01, 1]
MTCdel <- read.xlsx("tables/Supplementary Table 12.xlsx", startRow = 4, sheet = 9)
MTCdel <- MTCdel[MTCdel$`p-value` < 0.01 & MTCdel$`p-value.MTC` < 0.01, 1]
MTCmutamp <- read.xlsx("tables/Supplementary Table 12.xlsx", startRow = 4, sheet = 10)[, 1]
MTCmutdel <- read.xlsx("tables/Supplementary Table 12.xlsx", startRow = 4, sheet = 11)[, 1]

NEUmut <- read.xlsx("tables/Supplementary Table 12.xlsx", startRow = 4, sheet = 12)[, 1]
NEUamp <- read.xlsx("tables/Supplementary Table 12.xlsx", startRow = 4, sheet = 13)
NEUamp <- NEUamp[NEUamp$`p-value` < 0.01 & NEUamp$`p-value.NEU` < 0.01, 1]
NEUdel <- read.xlsx("tables/Supplementary Table 12.xlsx", startRow = 4, sheet = 14)
NEUdel <- NEUdel[NEUdel$`p-value` < 0.01 & NEUdel$`p-value.NEU` < 0.01, 1]
NEUmutamp <- read.xlsx("tables/Supplementary Table 12.xlsx", startRow = 4, sheet = 15)[, 1]
NEUmutdel <- read.xlsx("tables/Supplementary Table 12.xlsx", startRow = 4, sheet = 16)[, 1]

PPRmut <- read.xlsx("tables/Supplementary Table 12.xlsx", startRow = 4, sheet = 17)[, 1]
PPRamp <- read.xlsx("tables/Supplementary Table 12.xlsx", startRow = 4, sheet = 18)
PPRamp <- PPRamp[PPRamp$`p-value` < 0.01 & PPRamp$`p-value.PPR` < 0.01, 1]
PPRdel <- read.xlsx("tables/Supplementary Table 12.xlsx", startRow = 4, sheet = 19)
PPRdel <- PPRdel[PPRdel$`p-value` < 0.01 & PPRdel$`p-value.PPR` < 0.01, 1]
PPRmutamp <- read.xlsx("tables/Supplementary Table 12.xlsx", startRow = 4, sheet = 20)[, 1]
PPRmutdel <- read.xlsx("tables/Supplementary Table 12.xlsx", startRow = 4, sheet = 21)[, 1]

#Select amp/del genes 2-fold frequency in group
tmp <- read.xlsx("tables/Supplementary Table 12.xlsx", startRow = 4, sheet = 1)
groups <- tmp$kNN.classification
names(groups) <- tmp$TCGA.ID
groups <- groups[groups != "Undefined"]

load("RData/TCGA_CNVfunc.RData", verbose = T)
cs <- intersect(colnames(CNVfunc), names(groups))
CNVfunc <- CNVfunc[, cs]
groups <- groups[cs]

CNVDEL <- CNVfunc
CNVDEL[CNVDEL > 0] <- 0
CNVDEL[CNVDEL < 0] <- 1

toTake <- rowMeans(CNVDEL[GPMdel, names(groups)[groups == "GPM"]])/rowMeans(CNVDEL[GPMdel, names(groups)[groups != "GPM"]])
toTake <- toTake[toTake > 2]
GPMdel <- names(sort(toTake))

toTake <- rowMeans(CNVDEL[MTCdel, names(groups)[groups == "MTC"]])/rowMeans(CNVDEL[MTCdel, names(groups)[groups != "MTC"]])
toTake <- toTake[toTake > 2]
MTCdel <- names(sort(toTake))

toTake <- rowMeans(CNVDEL[NEUdel, names(groups)[groups == "NEU"]])/rowMeans(CNVDEL[NEUdel, names(groups)[groups != "NEU"]])
toTake <- toTake[toTake > 2]
NEUdel <- names(sort(toTake))

toTake <- rowMeans(CNVDEL[PPRdel, names(groups)[groups == "PPR"]])/rowMeans(CNVDEL[PPRdel, names(groups)[groups != "PPR"]])
toTake <- toTake[toTake > 2]
PPRdel <- names(sort(toTake))

CNVAMP <- CNVfunc
CNVAMP[CNVAMP < 0] <- 0
CNVAMP[CNVAMP > 0] <- 1

toTake <- rowMeans(CNVAMP[GPMamp, names(groups)[groups == "GPM"]])/rowMeans(CNVAMP[GPMamp, names(groups)[groups != "GPM"]])
toTake <- toTake[toTake > 2]
GPMamp <- names(sort(toTake))

toTake <- rowMeans(CNVAMP[MTCamp, names(groups)[groups == "MTC"]])/rowMeans(CNVAMP[MTCamp, names(groups)[groups != "MTC"]])
toTake <- toTake[toTake > 2]
MTCamp <- names(sort(toTake))

toTake <- rowMeans(CNVAMP[NEUamp, names(groups)[groups == "NEU"]])/rowMeans(CNVAMP[NEUamp, names(groups)[groups != "NEU"]])
toTake <- toTake[toTake > 2]
NEUamp <- names(sort(toTake))

toTake <- rowMeans(CNVAMP[PPRamp, names(groups)[groups == "PPR"]])/rowMeans(CNVAMP[PPRamp, names(groups)[groups != "PPR"]])
toTake <- toTake[toTake > 2]
PPRamp <- names(sort(toTake))

###remove mut or amp/del when the gene is mut/amp or mut/del
allMUT <- unique(c(GPMmut, MTCmut, NEUmut, PPRmut))
allAMP <- unique(c(GPMamp, MTCamp, NEUamp, PPRamp))
allDEL <- unique(c(GPMdel, MTCdel, NEUdel, PPRdel))
allMUTAMP <- unique(c(GPMmutamp, MTCmutamp, NEUmutamp, PPRmutamp))
allMUTDEL <- unique(c(GPMmutdel, MTCmutdel, NEUmutdel, PPRmutdel))

toRemove_mut <- intersect(allMUT, c(allMUTAMP, allMUTDEL))
toRemove_amp <- intersect(allAMP, allMUTAMP)
toRemove_del <- intersect(allDEL, allMUTDEL)

allMUT <- allMUT[!allMUT %in% toRemove_mut]
allAMP <- allAMP[!allAMP %in% toRemove_amp]
allDEL <- allDEL[!allDEL %in% toRemove_del]

GPMmut <- GPMmut[GPMmut %in% allMUT]
MTCmut <- MTCmut[MTCmut %in% allMUT]
NEUmut <- NEUmut[NEUmut %in% allMUT]
PPRmut <- PPRmut[PPRmut %in% allMUT]

GPMamp <- GPMamp[GPMamp %in% allAMP]
MTCamp <- MTCamp[MTCamp %in% allAMP]
NEUamp <- NEUamp[NEUamp %in% allAMP]
PPRamp <- PPRamp[PPRamp %in% allAMP]

GPMdel <- GPMdel[GPMdel %in% allDEL]
MTCdel <- MTCdel[MTCdel %in% allDEL]
NEUdel <- NEUdel[NEUdel %in% allDEL]
PPRdel <- PPRdel[PPRdel %in% allDEL]


###create matrix for each alt type
tmp <- read.xlsx("tables/Supplementary Table 12.xlsx", startRow = 4, sheet = 1)
groups <- tmp$kNN.classification
names(groups) <- tmp$TCGA.ID
groups <- groups[groups != "Undefined"]

load("RData/TCGA_CNVfunc.RData", verbose = T)

load("RData/TCGA_mutMat.RData", verbose = T)
load("RData/TCGA_samplesHyperMut.RData", verbose = T)
mutMat <- mutMat[, !colnames(mutMat)%in%samplesHiper]

allSamples <- names(groups)
allSamples <- unique(c(allSamples[allSamples %in% colnames(mutMat)], allSamples[allSamples %in% colnames(CNVfunc)]))
groups <- groups[allSamples]

noMut <- setdiff(allSamples, colnames(mutMat))
noCNV <- setdiff(allSamples, colnames(CNVfunc))

#AMP
cs <- intersect(names(groups), colnames(CNVfunc))
CNVfunc <- CNVfunc[, cs]

CNVAMP <- CNVfunc[c(GPMamp, MTCamp, NEUamp, PPRamp), ]
CNVAMP[CNVAMP < 0] <- 0
CNVAMP[CNVAMP > 0] <- 1
rownames(CNVAMP) <- paste(rownames(CNVAMP), "amp", sep = "")
GPMamp <- paste(GPMamp, "amp", sep = "")
MTCamp <- paste(MTCamp, "amp", sep = "")
NEUamp <- paste(NEUamp, "amp", sep = "")
PPRamp <- paste(PPRamp, "amp", sep = "")

colAMP <- c(rep("red", length(GPMamp)), rep("green3", length(MTCamp)), rep("blue", length(NEUamp)), rep("cyan", length(PPRamp)))
names(colAMP) <- rownames(CNVAMP)

toAdd <- matrix(0, nrow=nrow(CNVAMP), ncol=length(setdiff(allSamples, colnames(CNVfunc))))
colnames(toAdd) <- setdiff(allSamples, colnames(CNVfunc))
CNVAMP <- cbind(CNVAMP, toAdd)
CNVAMP <- CNVAMP[, allSamples]

#DEL
CNVDEL <- CNVfunc[c(GPMdel, MTCdel, NEUdel, PPRdel), ]
CNVDEL[CNVDEL > 0] <- 0
CNVDEL[CNVDEL < 0] <- -1
rownames(CNVDEL) <- paste(rownames(CNVDEL), "del", sep = "")
GPMdel <- paste(GPMdel, "del", sep = "")
MTCdel <- paste(MTCdel, "del", sep = "")
NEUdel <- paste(NEUdel, "del", sep = "")
PPRdel <- paste(PPRdel, "del", sep = "")

colDEL <- c(rep("red", length(GPMdel)), rep("green3", length(MTCdel)), rep("blue", length(NEUdel)), rep("cyan", length(PPRdel)))
names(colDEL) <- rownames(CNVDEL)

toAdd <- matrix(0, nrow=nrow(CNVDEL), ncol=length(setdiff(allSamples, colnames(CNVfunc))))
colnames(toAdd) <- setdiff(allSamples, colnames(CNVfunc))
CNVDEL <- cbind(CNVDEL, toAdd)
CNVDEL <- CNVDEL[, allSamples]

#MUT
cs <- intersect(names(groups), colnames(mutMat))
mutMat <- mutMat[, cs]

MUT <- mutMat[c(GPMmut, MTCmut, NEUmut, PPRmut), ]
MUT[MUT != 0] <- 3
rownames(MUT) <- paste(rownames(MUT), "mut", sep = "")
GPMmut <- paste(GPMmut, "mut", sep = "")
MTCmut <- paste(MTCmut, "mut", sep = "")
NEUmut <- paste(NEUmut, "mut", sep = "")
PPRmut <- paste(PPRmut, "mut", sep = "")

colMUT <- c(rep("red", length(GPMmut)), rep("green3", length(MTCmut)), rep("blue", length(NEUmut)), rep("cyan", length(PPRmut)))
names(colMUT) <- rownames(MUT)

toAdd <- matrix(0, nrow=nrow(MUT), ncol=length(setdiff(allSamples, colnames(mutMat))))
colnames(toAdd) <- setdiff(allSamples, colnames(MUT))
MUT <- cbind(MUT, toAdd)
MUT <- MUT[, allSamples]


#MUT AMP
toAdd <- matrix(0, nrow=nrow(CNVfunc), ncol=length(setdiff(allSamples, colnames(CNVfunc))))
colnames(toAdd) <- setdiff(allSamples, colnames(CNVfunc))
CNVfunc <- cbind(CNVfunc, toAdd)
CNVfunc <- CNVfunc[, allSamples]

toAdd <- matrix(0, nrow=nrow(mutMat), ncol=length(setdiff(allSamples, colnames(mutMat))))
colnames(toAdd) <- setdiff(allSamples, colnames(mutMat))
mutMat <- cbind(mutMat, toAdd)
mutMat <- mutMat[, allSamples]

tmpMut <- mutMat[c(GPMmutamp, MTCmutamp, NEUmutamp, PPRmutamp), ]
tmpMut[tmpMut != 0] <- 3
tmpAmp <- CNVfunc[c(GPMmutamp, MTCmutamp, NEUmutamp, PPRmutamp), ]
tmpAmp[tmpAmp < 0] <- 0
tmpAmp[tmpAmp > 0] <- 1
MUTAMP <- tmpMut + tmpAmp
MUTAMP[MUTAMP == 4] <- 2

rownames(MUTAMP) <- paste(rownames(MUTAMP), "mutamp", sep = "")
GPMmutamp <- paste(GPMmutamp, "mutamp", sep = "")
MTCmutamp <- paste(MTCmutamp, "mutamp", sep = "")
NEUmutamp <- paste(NEUmutamp, "mutamp", sep = "")
PPRmutamp <- paste(PPRmutamp, "mutamp", sep = "")

colMUTAMP <- c(rep("red", length(GPMmutamp)), rep("green3", length(MTCmutamp)), rep("blue", length(NEUmutamp)), rep("cyan", length(PPRmutamp)))
names(colMUTAMP) <- rownames(MUTAMP)


#MUT DEL
tmpMut <- mutMat[c(GPMmutdel, MTCmutdel, NEUmutdel, PPRmutdel), ]
tmpMut[tmpMut != 0] <- 3
tmpDel <- CNVfunc[c(GPMmutdel, MTCmutdel, NEUmutdel, PPRmutdel), ]
tmpDel[tmpDel > 0] <- 0
tmpDel[tmpDel < 0] <- 1
MUTDEL <- tmpMut + tmpDel
MUTDEL[MUTDEL == 1] <- -1
MUTDEL[MUTDEL == 4] <- -2

rownames(MUTDEL) <- paste(rownames(MUTDEL), "mutdel", sep = "")
GPMmutdel <- paste(GPMmutdel, "mutdel", sep = "")
MTCmutdel <- paste(MTCmutdel, "mutdel", sep = "")
NEUmutdel <- paste(NEUmutdel, "mutdel", sep = "")
PPRmutdel <- paste(PPRmutdel, "mutdel", sep = "")

colMUTDEL <- c(rep("red", length(GPMmutdel)), rep("green3", length(MTCmutdel)), rep("blue", length(NEUmutdel)), rep("cyan", length(PPRmutdel)))
names(colMUTDEL) <- rownames(MUTDEL)

#integrated matrix
allGenes <- c(rownames(MUT), rownames(CNVAMP), rownames(CNVDEL), rownames(MUTAMP), rownames(MUTDEL))

bigHM <- matrix(0, nrow=length(allGenes), ncol=length(allSamples))
rownames(bigHM) <- allGenes
colnames(bigHM) <- allSamples

bigHM[rownames(MUT), ] <- MUT
bigHM[rownames(CNVAMP),] <- CNVAMP
bigHM[rownames(CNVDEL),] <- CNVDEL
bigHM[rownames(MUTAMP),] <- MUTAMP
bigHM[rownames(MUTDEL),] <- MUTDEL
table(bigHM)

groups <- sort(groups)

RowSideColors <- c(colMUT, colAMP, colDEL, colMUTAMP, colMUTDEL)
RowSideColors <- c(RowSideColors[RowSideColors == "red"], RowSideColors[RowSideColors == "green3"], RowSideColors[RowSideColors == "blue"], RowSideColors[RowSideColors == "cyan"])

toPlot <- bigHM[names(RowSideColors), names(groups)]

RowSideColors <- rbind(RowSideColors)
colnames(RowSideColors) <- NULL

classCol <- c("green3", "red", "blue", "cyan")
names(classCol) <- c("MTC", "GPM", "NEU", "PPR")
classCol <- classCol[groups]

NOmut <- rep("white", length(colnames(toPlot)))
names(NOmut) <- colnames(toPlot)
NOmut[noMut] <- "gray23"

NOcn <- rep("white", length(colnames(toPlot)))
names(NOcn) <- colnames(toPlot)
NOcn[noCNV] <- "gray23"

ColSideColors <- cbind(NOcn, NOmut, classCol)
rownames(ColSideColors) <- NULL

source("functions.R")
library(gplots)

heatmap.3(toPlot, 
          ColSideColors = ColSideColors,  
          RowSideColors = RowSideColors,
          key = FALSE, 
          Colv = NA,
          Rowv = NA,
          scale = "none", 
          col = colorRampPalette(c("turquoise1", "blue4", "snow", "firebrick", "darkorange", "green2"))(6), 
          dendrogram = "none",
          labRow = NA, labCol = NA,
          side.height.fraction = 0.3, keysize = 0.1, cexRow = 0.3,
)


### Figure 4b
# AMPLIFICATIONS
load("RData/listOfGO_enrich_mutAMP.RData")

path <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION","GO_MITOCHONDRION","HALLMARK_GLYCOLYSIS","GO_REGULATION_OF_GLUCOSE_TRANSPORT", "GO_REGULATION_OF_CARBOHYDRATE_METABOLIC_PROCESS", "GO_REGULATION_OF_CARBOHYDRATE_BIOSYNTHETIC_PROCESS","HALLMARK_HYPOXIA", "GO_REGULATION_OF_LIPID_METABOLIC_PROCESS","GO_LIPID_STORAGE", "GO_LIPID_OXIDATION","GO_POSITIVE_REGULATION_OF_FATTY_ACID_METABOLIC_PROCESS","GO_FATTY_ACID_CATABOLIC_PROCESS", "GO_CELLULAR_RESPONSE_TO_AMINO_ACID_STIMULUS", "GO_REGULATION_OF_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS")

pval_GPM <- listOfGO_enrich_mutAMP$GPM[path, "pval"]
pval_MTC <- listOfGO_enrich_mutAMP$MTC[path, "pval"]
pval_NEU <- listOfGO_enrich_mutAMP$NEU[path, "pval"]
pval_PPR <- listOfGO_enrich_mutAMP$PPR[path, "pval"]

estim_GPM <- listOfGO_enrich_mutAMP$GPM[path, "estim"]
estim_MTC <- listOfGO_enrich_mutAMP$MTC[path, "estim"]
estim_NEU <- listOfGO_enrich_mutAMP$NEU[path, "estim"]
estim_PPR <- listOfGO_enrich_mutAMP$PPR[path, "estim"]

M <- cbind(pval_GPM, pval_MTC, pval_NEU, pval_PPR)
Mstat <- cbind(estim_GPM, estim_MTC, estim_NEU, estim_PPR)

library(RColorBrewer)
M2 <- -log10(M)
M2[Mstat < 1] <- -M2[Mstat < 1]

colnames(M2) <- c("GPM", "MTC", "NEU", "PPR")
x <- brewer.pal(n=9, name="YlGnBu")[1:9]
rownames(M2)[grep(rownames(M2), pattern = "GO_")] <- unlist(lapply(strsplit(rownames(M2)[grep(rownames(M2), pattern = "GO_")], split = "GO_"), function(x) x[[2]]))
rownames(M2)[grep(rownames(M2), pattern = "REGULATION_OF")] <- unlist(lapply(strsplit(rownames(M2)[grep(rownames(M2), pattern = "REGULATION_OF")], split = "REGULATION_OF_"), function(x) x[[2]]))
rownames(M2)[grep(rownames(M2), pattern = "HALLMARK_")] <- unlist(lapply(strsplit(rownames(M2)[grep(rownames(M2), pattern = "HALLMARK_")], split = "HALLMARK_"), function(x) x[[2]]))

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
rownames(M2) <- firstup(tolower(gsub(rownames(M2), pattern = "_", replacement = " ")))

library(corrplot)
library(RColorBrewer)

corrplot(M2, method="color", is.corr=F, p.mat = M, sig.level = c(0.01, 0.05, 0.10), insig = "label_sig", pch.cex = 2,
         col=colorRampPalette(c("blue4", "white", "darkred"))(51), tl.col = "black", cl.pos = "b", tl.srt = 35, tl.pos='n', 
         cl.lim=c(-14,14))

# DELETIONS
load("RData/listOfGO_enrich_mutDEL.RData")

path <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION","GO_MITOCHONDRION","GO_POSITIVE_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS","GO_POSITIVE_REGULATION_OF_GLUCOSE_TRANSPORT", "GO_REGULATION_OF_CARBOHYDRATE_METABOLIC_PROCESS", "GO_REGULATION_OF_CARBOHYDRATE_BIOSYNTHETIC_PROCESS", "HALLMARK_HYPOXIA", "GO_REGULATION_OF_LIPID_METABOLIC_PROCESS","GO_LIPID_STORAGE", "GO_LIPID_OXIDATION", "GO_POSITIVE_REGULATION_OF_FATTY_ACID_METABOLIC_PROCESS","GO_FATTY_ACID_CATABOLIC_PROCESS", "GO_CELLULAR_RESPONSE_TO_AMINO_ACID_STIMULUS", "GO_REGULATION_OF_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS")

pval_GPM <- listOfGO_enrich_mutDEL$GPM[path, "pval"]
pval_MTC <- listOfGO_enrich_mutDEL$MTC[path, "pval"]
pval_NEU <- listOfGO_enrich_mutDEL$NEU[path, "pval"]
pval_PPR <- listOfGO_enrich_mutDEL$PPR[path, "pval"]

estim_GPM <- listOfGO_enrich_mutDEL$GPM[path, "estim"]
estim_MTC <- listOfGO_enrich_mutDEL$MTC[path, "estim"]
estim_NEU <- listOfGO_enrich_mutDEL$NEU[path, "estim"]
estim_PPR <- listOfGO_enrich_mutDEL$PPR[path, "estim"]

M <- cbind(pval_GPM, pval_MTC, pval_NEU, pval_PPR)
Mstat <- cbind(estim_GPM, estim_MTC, estim_NEU, estim_PPR)

library(RColorBrewer)
M2 <- -log10(M)
M2[Mstat < 1] <- -M2[Mstat < 1]

colnames(M2) <- c("GPM", "MTC", "NEU", "PPR")
x <- brewer.pal(n=9, name="YlGnBu")[1:9]
rownames(M2)[grep(rownames(M2), pattern = "GO_")] <- unlist(lapply(strsplit(rownames(M2)[grep(rownames(M2), pattern = "GO_")], split = "GO_"), function(x) x[[2]]))
rownames(M2)[grep(rownames(M2), pattern = "REGULATION_OF")] <- unlist(lapply(strsplit(rownames(M2)[grep(rownames(M2), pattern = "REGULATION_OF")], split = "REGULATION_OF_"), function(x) x[[2]]))
rownames(M2)[grep(rownames(M2), pattern = "HALLMARK_")] <- unlist(lapply(strsplit(rownames(M2)[grep(rownames(M2), pattern = "HALLMARK_")], split = "HALLMARK_"), function(x) x[[2]]))

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
rownames(M2) <- firstup(tolower(gsub(rownames(M2), pattern = "_", replacement = " ")))

library(corrplot)
library(RColorBrewer)


corrplot(M2, method="color", is.corr=F, p.mat = M, sig.level = c(0.01, 0.05, 0.10), insig = "label_sig", pch.cex = 2,
         col=colorRampPalette(c("blue4", "white", "darkred"))(51), tl.col = "black", cl.pos = "b", tl.srt = 35, tl.pos='n', 
         cl.lim=c(-14,14))


### Figure 4c
#UPPER PANEL
# upload the TCGA survival data
library(survival)
load("RData/TCGAGBM_clinicalData.RData")

# upload the TCGA kNN classification relative to GPM and MTC subtypes
library(openxlsx)
TCGAclassification <- read.xlsx("tables/Supplementary Table 12.xlsx", sheet = 1)
TCGAclassification <- TCGAclassification[2:nrow(TCGAclassification), ]
colnames(TCGAclassification) <- TCGAclassification[1, ]
TCGAclassification <- TCGAclassification[-1, ]
rownames(TCGAclassification) <- TCGAclassification$`TCGA ID`
TCGAclassification2 <- TCGAclassification$`kNN classification`
names(TCGAclassification2) <- rownames(TCGAclassification)
table(TCGAclassification2)
TCGAclassification2 <- c(TCGAclassification2[TCGAclassification2 == "MTC"], TCGAclassification2[TCGAclassification2 == "GPM"])

#upload NES activity GPM and MTC
library(openxlsx)
NES <- read.xlsx("tables/Supplementary Table 13.xlsx", sheet = 1)
NES <- NES[3:nrow(NES), ]
colnames(NES) <- NES[1, ]
NES <- NES[-1, ]
rownames(NES) <- NES[,1]

NES_GPM <- as.numeric(NES[, "GPM NES"])
names(NES_GPM) <- rownames(NES)
NES_MTC <-  as.numeric(NES[, "MTC NES"])
names(NES_MTC) <- rownames(NES)

#upload functional CN calls in TCGA tumors
load("RData/TCGA_CNVfunc.RData", verbose = T)
CNVfuncTCGA <- CNVfunc
ss <- intersect(intersect(intersect(names(TCGAclassification2), names(NES_GPM)), rownames(stData)), colnames(CNVfuncTCGA))

TCGAclassification2 <- TCGAclassification2[ss]
NES_GPM <- NES_GPM[ss]
NES_MTC <- NES_MTC[ss]
stData <- stData[ss, ]
CNVfuncTCGA <- CNVfuncTCGA[,ss]

# derive the difference of GPM activity versus MTC activity
diffGPM_MTCact <- NES_GPM - NES_MTC
diffGPM_MTCact_ord <- sort(diffGPM_MTCact, decreasing = F)

# derive the cox model
aCox <- coxph(stData[names(diffGPM_MTCact_ord)] ~ diffGPM_MTCact_ord)
aPred <- predict(aCox, type = "risk", se.fit = T)

hr <- aPred$fit
high <- hr + 2*aPred$se.fit
low <- hr - 2*aPred$se.fit
x <- 1:length(diffGPM_MTCact)

# Fig. 4c UPPER PANEL
par(mar=c(0,0,0,0))
plot(x, hr, type = "n", ylim = c(0.5, 1.65), xlim=c(0.5, 273.5), xaxs="i", yaxs="i", xaxt = "n", axes = F)
polygon(c(x, rev(x)), c(low, rev(high)), col = "#DEEBF7",border=NA)
lines(x, hr, lwd = 2)

#MIDDLE PANEL
# upload the TCGA survival data
library(survival)
load("RData/TCGAGBM_clinicalData.RData", verbose = T)

# upload the TCGA kNN classification relative to GPM and MTC subtypes
library(openxlsx)
TCGAclassification <- read.xlsx("tables/Supplementary Table 12.xlsx", sheet = 1)
TCGAclassification <- TCGAclassification[2:nrow(TCGAclassification), ]
colnames(TCGAclassification) <- TCGAclassification[1, ]
TCGAclassification <- TCGAclassification[-1, ]
rownames(TCGAclassification) <- TCGAclassification$`TCGA ID`
TCGAclassification2 <- TCGAclassification$`kNN classification`
names(TCGAclassification2) <- rownames(TCGAclassification)
table(TCGAclassification2)
TCGAclassification2 <- c(TCGAclassification2[TCGAclassification2 == "MTC"], TCGAclassification2[TCGAclassification2 == "GPM"])

#upload NES activity GPM and MTC
library(openxlsx)
NES <- read.xlsx("tables/Supplementary Table 13.xlsx", sheet = 1)
NES <- NES[3:nrow(NES), ]
colnames(NES) <- NES[1, ]
NES <- NES[-1, ]
rownames(NES) <- NES[,1]

NES_GPM <- as.numeric(NES[, "GPM NES"])
names(NES_GPM) <- rownames(NES)
NES_MTC <-  as.numeric(NES[, "MTC NES"])
names(NES_MTC) <- rownames(NES)

#upload functional CN calls in TCGA tumors
load("RData/TCGA_CNVfunc.RData", verbose = T)
CNVfuncTCGA <- CNVfunc
ss <- intersect(intersect(intersect(names(TCGAclassification2), names(NES_GPM)), rownames(stData)), colnames(CNVfuncTCGA))

TCGAclassification2 <- TCGAclassification2[ss]
NES_GPM <- NES_GPM[ss]
NES_MTC <- NES_MTC[ss]
stData <- stData[ss, ]
CNVfuncTCGA <- CNVfuncTCGA[,ss]

# derive the difference of GPM activity versus MTC activity
diffGPM_MTCact <- NES_GPM - NES_MTC
diffGPM_MTCact_ord <- sort(diffGPM_MTCact, decreasing = F)

NES_GPM <- NES_GPM[names(diffGPM_MTCact_ord)]
NES_MTC <- NES_MTC[names(diffGPM_MTCact_ord)]

par(mar=c(0,0,0,0))
range(NES_MTC)
range(NES_GPM)
plot(1:length(NES_MTC), NES_MTC, type = "n", xlim = c(0, length(NES_MTC)+1), ylim = c(-3.5, 5), xaxs = "i", xaxt = "n", axes = F)
plot(1:length(NES_MTC), NES_MTC, type = "n", xlim = c(0, length(NES_MTC)+1), ylim = c(-3.5, 5), xaxs = "i", xaxt = "n")
abline(lm(NES_MTC ~ order(sort(diffGPM_MTCact, decreasing = F))), lwd = 2)
abline(lm(NES_GPM ~ order(sort(diffGPM_MTCact, decreasing = F))), lwd = 2)
points(NES_MTC, col = "blue4", pch = 16, cex=1.5)
points(NES_GPM, col = "darkred", pch = 16, cex=1.5)

cor.test(NES_GPM, NES_MTC,  method = "spearman", alternative = "two.sided")

#LOWER PANEL
# upload the TCGA survival data
library(survival)
load("RData/TCGAGBM_clinicalData.RData")

# upload the TCGA kNN classification relative to GPM and MTC subtypes
library(openxlsx)
TCGAclassification <- read.xlsx("tables/Supplementary Table 12.xlsx", sheet = 1)
TCGAclassification <- TCGAclassification[2:nrow(TCGAclassification), ]
colnames(TCGAclassification) <- TCGAclassification[1, ]
TCGAclassification <- TCGAclassification[-1, ]
rownames(TCGAclassification) <- TCGAclassification$`TCGA ID`
TCGAclassification2 <- TCGAclassification$`kNN classification`
names(TCGAclassification2) <- rownames(TCGAclassification)
table(TCGAclassification2)
TCGAclassification2 <- c(TCGAclassification2[TCGAclassification2 == "MTC"], TCGAclassification2[TCGAclassification2 == "GPM"])

#upload NES activity GPM and MTC
library(openxlsx)
NES <- read.xlsx("tables/Supplementary Table 13.xlsx", sheet = 1)
NES <- NES[3:nrow(NES), ]
colnames(NES) <- NES[1, ]
NES <- NES[-1, ]
rownames(NES) <- NES[,1]

NES_GPM <- as.numeric(NES[, "GPM NES"])
names(NES_GPM) <- rownames(NES)
NES_MTC <-  as.numeric(NES[, "MTC NES"])
names(NES_MTC) <- rownames(NES)

#upload functional CN calls in TCGA tumors
load("RData/TCGA_CNVfunc.RData", verbose = T)
CNVfuncTCGA <- CNVfunc
ss <- intersect(intersect(intersect(names(TCGAclassification2), names(NES_GPM)), rownames(stData)), colnames(CNVfuncTCGA))

TCGAclassification2 <- TCGAclassification2[ss]
NES_GPM <- NES_GPM[ss]
NES_MTC <- NES_MTC[ss]
stData <- stData[ss, ]
CNVfuncTCGA <- CNVfuncTCGA[,ss]

# derive the difference of GPM activity versus MTC activity
diffGPM_MTCact <- NES_GPM - NES_MTC
diffGPM_MTCact_ord <- sort(diffGPM_MTCact, decreasing = F)

#upload genes identified as glycolytic and mitochondrial associated to GPM and MTC subtypes

library(openxlsx)
genetAlt_MTCGPM <- read.xlsx("tables/Supplementary Table 13.xlsx", sheet = 2)
genetAlt_MTCGPM <- genetAlt_MTCGPM[2:nrow(genetAlt_MTCGPM), ]
colnames(genetAlt_MTCGPM) <- genetAlt_MTCGPM[1, ]
genetAlt_MTCGPM <- genetAlt_MTCGPM[-1, ]

AMP_mitochgenes <- genetAlt_MTCGPM$`Amplified Mitochondrial genes`
DEL_glycolgenes <- genetAlt_MTCGPM$`Deleted Glycolytic genes`[!is.na(genetAlt_MTCGPM$`Deleted Glycolytic genes`)]
AMP_glycolgenes <- genetAlt_MTCGPM$`Amplified Glycolytic genes`[!is.na(genetAlt_MTCGPM$`Amplified Glycolytic genes`)]
DEL_mitochgenes <- genetAlt_MTCGPM$`Deleted Mitochondrial genes`[!is.na(genetAlt_MTCGPM$`Deleted Mitochondrial genes`)]

CNVfunc2AMP <- CNVfuncTCGA[c(AMP_mitochgenes,AMP_glycolgenes), ]
table(CNVfunc2AMP)
CNVfunc2AMP[CNVfunc2AMP == -2] <- 0
CNVfunc2AMP[CNVfunc2AMP == -1 ] <- 0
CNVfunc2AMP[CNVfunc2AMP == 2] <- 1
CNVfunc2AMP <- CNVfunc2AMP[, names(diffGPM_MTCact_ord)]

CNVfunc2DEL <- CNVfuncTCGA[c(DEL_mitochgenes,DEL_glycolgenes), ]
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

SupportPlot <- CNVfuncTCGA[c("PTEN","RB1"), names(diffGPM_MTCact_ord)]
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


### Figure 4d
source("functions.R")

load("RData/TCGA_geData.RData", verbose = T)
load("RData/TCGA_methData.RData", verbose = T)

library(openxlsx)
classTab <- read.xlsx("tables/Supplementary Table 6.xlsx", startRow = 3, sheet = 3)

groups <- classTab$Core.classification
names(groups) <- classTab$TCGA.ID

groupsMeth <- groups[colnames(methData)]

gg <- names(table(groups))
res <- vector("list", length(gg))
names(res) <- gg
for(i in 1:length(gg)){
  (whichOfInterest <- gg[i])
  
  ### degs geData
  gruppo1 <- names(groups)[groups == whichOfInterest]
  gruppo2 <- names(groups)[groups != whichOfInterest]
  ans <- apply(geData, 1, function(x) wilcox.test(x[gruppo1], x[gruppo2])$p.value)
  fc <- rowMeans(geData[, gruppo1]) - rowMeans(geData[, gruppo2])
  DEGS <- cbind(fc, ans)
  DEGS <- DEGS[order(DEGS[,"fc"], decreasing = TRUE),]
  fdr <- p.adjust(DEGS[, "ans"], method = "fdr")
  DEGS <- cbind(DEGS, fdr)
  
  DEGS_geData <- as.data.frame(DEGS)
  colnames(DEGS_geData)[2:3] <- c("pValue", "qValue")
  colnames(DEGS_geData) <- paste0("ge_", colnames(DEGS_geData))
  DEGS_geData <- cbind(GeneSymbol = rownames(DEGS_geData), DEGS_geData, stringsAsFactors = F)
  
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
  
  probesInfo <- probesInfo[probesInfo$Gene_Name %in% rownames(DEGS_geData), ]
  DEGS_methData <- DEGS_methData[probesInfo$ProbID, ]
  DEGS_geData <- DEGS_geData[probesInfo$Gene_Name, ]
  
  DEGS <- data.frame(ProbeID = DEGS_methData$ProbeID, GeneSymbol = DEGS_geData$GeneSymbol, DEGS_methData[, -1], DEGS_geData[, -1], stringsAsFactors = F)
  rownames(DEGS) <- NULL
  
  res[[i]] <- DEGS
}

par(mfrow = c(2, 2))
for(i in 1:length(res)){
  DEGsSB <- starburstPlot(res[[i]], cut.off.p = 0.05, cut.off.fc.ge = .4, cut.off.fc.meth = .3, porq = "p")
}



