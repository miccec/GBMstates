### ED Figure 9a

#UPPER PANEL
load("RData/NES_PDCs_GBMsubtypes.RData")

# upload the classification of PDCs dataset
load("RData/PDCsclassification.RData")

# generate the upper panel 
left <- c(0, 0.22, 0.47, 0.61)+0.10
rigth <- c(0.22, 0.47, 0.61, 0.78)+0.10
bottom <- rep(0.1, 4)
top <- rep(1, 4)

m <- cbind(left, rigth, bottom, top)
split.screen(m)

screen(1)
par(mar = c(2, 2, 5, 0))
gg <- "GPM"
x <- list(Prolif = NES_PDCs_PPR[names(PDCsclassification)[PDCsclassification == gg]], 
          Neuro = NES_PDCs_NEU[names(PDCsclassification)[PDCsclassification == gg]], 
          Mito = NES_PDCs_MTC[names(PDCsclassification)[PDCsclassification == gg]], 
          Glyco = NES_PDCs_GPM[names(PDCsclassification)[PDCsclassification == gg]])
mmeans <- sapply(x, FUN = mean)
ssd <- sapply(x, FUN = sd)
ssd <- ssd/sqrt(sapply(x, length))
li <- mmeans - ssd #* 2
ui <- mmeans + ssd #* 2
barCenters <- barplot(height = mmeans,
                      beside = T, las = 2,
                      xlim=c(-1, 1.5),
                      cex.names = 0.75, yaxt = "n",cex.axis = 1.3,
                      border = "black", axes = TRUE, 
                      col = c("cyan", "blue", "green3", "red"), width = c(0.3, 0.3, 0.3, 0.3), space = 0.1, horiz = T)
segments(ui,barCenters, li, barCenters, lwd = 1.5)

arrows(ui, barCenters, li, barCenters,
       lwd = 1.5, angle = 90,
       code = 3, length = 0.05)

screen(2)
par(mar = c(2, 2, 5, 0))
gg <- "MTC"
x <- list(Prolif = NES_PDCs_PPR[names(PDCsclassification)[PDCsclassification == gg]], 
          Neuro = NES_PDCs_NEU[names(PDCsclassification)[PDCsclassification == gg]], 
          Mito = NES_PDCs_MTC[names(PDCsclassification)[PDCsclassification == gg]], 
          Glyco = NES_PDCs_GPM[names(PDCsclassification)[PDCsclassification == gg]])
mmeans <- sapply(x, FUN = mean)
ssd <- sapply(x, FUN = sd)
ssd <- ssd/sqrt(sapply(x, length))
li <- mmeans - ssd
ui <- mmeans + ssd
barCenters <- barplot(height = mmeans,
                      beside = T, las = 2,
                      xlim=c(-0.6, 0.6), 
                      cex.names = 0.75, yaxt = "n",cex.axis = 1.3,
                      border = "black", axes = TRUE,
                      col = c("cyan", "blue", "green3", "red"), width = c(0.3, 0.3, 0.3, 0.3), space = 0.1, horiz = T)
segments(ui, barCenters, li, barCenters, lwd = 1.5)
arrows(ui, barCenters, li, barCenters, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)

screen(3)
par(mar = c(2, 2, 5, 0))
gg <- "NEU"
x <- list(Prolif = NES_PDCs_PPR[names(PDCsclassification)[PDCsclassification == gg]], 
          Neuro = NES_PDCs_NEU[names(PDCsclassification)[PDCsclassification == gg]], 
          Mito = NES_PDCs_MTC[names(PDCsclassification)[PDCsclassification == gg]], 
          Glyco = NES_PDCs_GPM[names(PDCsclassification)[PDCsclassification == gg]])
mmeans <- sapply(x, FUN = mean)
ssd <- sapply(x, FUN = sd)
ssd <- ssd/sqrt(sapply(x, length))
li <- mmeans - ssd
ui <- mmeans + ssd
barCenters <- barplot(height = mmeans,
                      beside = T, las = 2,cex.axis = 1.3,
                      xlim=c(-1.1, 1.5),
                      cex.names = 0.75, yaxt = "n",
                      border = "black", axes = TRUE,
                      col = c("cyan", "blue", "green3", "red"), width = c(0.3, 0.3, 0.3, 0.3), space = 0.1, horiz=T)
segments(ui, barCenters, li, barCenters, lwd = 1.5)
arrows(ui, barCenters, li, barCenters, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)

screen(4)
par(mar = c(2, 2, 5, 0))
gg <- "PPR"
x <- list(Prolif = NES_PDCs_PPR[names(PDCsclassification)[PDCsclassification == gg]], 
          Neuro = NES_PDCs_NEU[names(PDCsclassification)[PDCsclassification == gg]], 
          Mito = NES_PDCs_MTC[names(PDCsclassification)[PDCsclassification == gg]], 
          Glyco = NES_PDCs_GPM[names(PDCsclassification)[PDCsclassification == gg]])
mmeans <- sapply(x, FUN = mean)
ssd <- sapply(x, FUN = sd)
ssd <- ssd/sqrt(sapply(x, length))
li <- mmeans - ssd
ui <- mmeans + ssd
barCenters <- barplot(height = mmeans,
                      beside = T, las = 2,cex.axis = 1.3,
                      xlim=c(-1, 2),
                      cex.names = 0.75, yaxt = "n",
                      border = "black", axes = TRUE,
                      col = c("cyan", "blue", "green3", "red"), width = c(0.3, 0.3, 0.3, 0.3), space = 0.1, horiz = T)
segments(ui, barCenters, li, barCenters, lwd = 1.5)
arrows(ui, barCenters, li, barCenters, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)

close.screen(all.screens = TRUE)

#MIDDLE PANEL 

#upload the NES of pathway in PDCs subtype
load("RData/NES_pathwaysPDCs.RData")

#upload the resulted significantly differential expressed pathways in the GBM subtypes
listOfDiffs <- vector("list", 4)
names(listOfDiffs) <- c("GPM", "MTC","NEU", "PPR")
listOfDiffs$GPM <- c("HALLMARK_GLYCOLYSIS", "GO_HEXOSE_METABOLIC_PROCESS", "HALLMARK_HYPOXIA", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "GO_LIPID_LOCALIZATION")
listOfDiffs$MTC <- c("GO_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_OXIDATIVE_PHOSPHORYLATION","Oxidative phosphorylation", "GO_ELECTRON_TRANSPORT_CHAIN", 
                     "GO_MITOCHONDRIAL_TRANSLATION", "GO_MITOCHONDRIAL_MEMBRANE_PART",
                     "GO_INNER_MITOCHONDRIAL_MEMBRANE_PROTEIN_COMPLEX")
listOfDiffs$NEU <-  c("GO_REGULATION_OF_NEUROTRANSMITTER_LEVELS","GO_NEUROTRANSMITTER_TRANSPORT", 
                      "GO_AXON_PART", 
                      "GO_PRESYNAPTIC_ACTIVE_ZONE","GO_PRESYNAPTIC_MEMBRANE", "GO_CELL_MORPHOGENESIS_INVOLVED_IN_NEURON_DIFFERENTIATION")
listOfDiffs$PPR <-  c("GO_CELL_CYCLE_PROCESS","GO_MISMATCH_REPAIR", 
                      "GO_MITOTIC_SPINDLE_ORGANIZATION",  "GO_CHROMATIN_ORGANIZATION")

#upload the PDCs classification
load("RData/PDCsclassification.RData")

ss <- intersect(colnames(NESpathways_PDCs), names(PDCsclassification))
PDCsclassification <- PDCsclassification[ss]
PDCsclassification <- PDCsclassification[c(which(PDCsclassification == "GPM"), which(PDCsclassification == "MTC"), 
                                           which(PDCsclassification == "NEU"), which(PDCsclassification == "PPR"))]
NESpathways_PDCs <- NESpathways_PDCs[, names(PDCsclassification)]
FDRpathways_PDCs <- FDRpathways_PDCs[, names(PDCsclassification)]
ss <- intersect(rownames(NESpathways_PDCs), as.character(unlist(listOfDiffs)))
listOfDiffs <- lapply(listOfDiffs, function(x) x[x%in%ss])
NESpathways_PDCs<- NESpathways_PDCs[as.character(unlist(listOfDiffs)), ]
FDRpathways_PDCs <- FDRpathways_PDCs[as.character(unlist(listOfDiffs)), ]

# generate the middle panel
toPlot <- NESpathways_PDCs

colRow <- rep(0, nrow(toPlot))
names(colRow) <- rownames(toPlot)
colRow[names(colRow)%in%listOfDiffs$GPM] <- "red"
colRow[names(colRow)%in%listOfDiffs$MTC] <- "green3"
colRow[names(colRow)%in%listOfDiffs$PPR] <- "cyan"
colRow[names(colRow)%in%listOfDiffs$NEU] <- "blue"
table(colRow)
row_annotation <- cbind(colRow)
colnames(row_annotation) <- NULL
rownames(row_annotation) <- NULL
RowSideColors = row_annotation

colCLUSTERS2 <- PDCsclassification[colnames(toPlot)]
tmp <- c("red","green3", "blue","cyan")
names(tmp) <- c("GPM", "MTC", "NEU", "PPR")
aa <- tmp[colCLUSTERS2]
names(aa) <- names(colCLUSTERS2)
colCLUSTERS2 <- aa
names(colCLUSTERS2) <- colnames(toPlot)
colCLUSTERS2[is.na(colCLUSTERS2)] <- "white"
table(colCLUSTERS2, useNA = "ifany")
annot2 <- cbind(colCLUSTERS2)
colnames(annot2) <- c("Clusters")
ColSideColors = annot2[, "Clusters"]

annotGenes <- c(rep("GPM", length(listOfDiffs$GPM)), rep("MTC", length(listOfDiffs$MTC)), 
                rep("NEU", length(listOfDiffs$NEU)), rep("PPR", length(listOfDiffs$PPR)))
names(annotGenes) <- as.character(unlist(listOfDiffs))

matExprs <- matrix(0, nrow=length(as.character(unlist(listOfDiffs))), ncol=ncol(toPlot))
rownames(matExprs) <- as.character(unlist(listOfDiffs))
colnames(matExprs) <- colnames(toPlot)

nesTresh <- 0.3
FDRTresh <- 0.05

for (i in 1:nrow(matExprs)){
  gene <- rownames(matExprs)[i]
  groupGene <- as.character(annotGenes[gene])
  concordSS <- names(which(NESpathways_PDCs[gene, names(which(PDCsclassification == groupGene))] > nesTresh & FDRpathways_PDCs[gene, names(which(PDCsclassification == groupGene))] < FDRTresh))
  DisconcordSS <- names(which(NESpathways_PDCs[gene, names(which(PDCsclassification != groupGene))] > nesTresh & FDRpathways_PDCs[gene, names(which(PDCsclassification != groupGene))] < FDRTresh))
  matExprs[gene, concordSS] <- 1
  matExprs[gene, DisconcordSS] <- 2
}

matExprs[names(which(annotGenes == "MTC")), names(which(PDCsclassification == "MTC"))][matExprs[names(which(annotGenes == "MTC")), names(which(PDCsclassification == "MTC"))] == 1] <- 3
matExprs[names(which(annotGenes == "NEU")), names(which(PDCsclassification == "NEU"))][matExprs[names(which(annotGenes == "NEU")), names(which(PDCsclassification == "NEU"))] == 1] <- 4
matExprs[names(which(annotGenes == "PPR")), names(which(PDCsclassification == "PPR"))][matExprs[names(which(annotGenes == "PPR")), names(which(PDCsclassification == "PPR"))] == 1] <- 5

GPM <- names(PDCsclassification[PDCsclassification == "GPM"])
MTC <- names(PDCsclassification[PDCsclassification == "MTC"])
NEU <- names(PDCsclassification[PDCsclassification == "NEU"])
PPR <- names(PDCsclassification[PDCsclassification == "PPR"])

CombinedConcord <- names(which(colSums(matExprs[c("HALLMARK_GLYCOLYSIS", "GO_HEXOSE_METABOLIC_PROCESS"), GPM]) > 0))
others <- setdiff(colnames(NES2), GPM)
CombinedDisConcord <- names(which(colSums(matExprs[c("HALLMARK_GLYCOLYSIS", "GO_HEXOSE_METABOLIC_PROCESS"), others]) > 0))
matExprs2 <- matExprs[!rownames(matExprs)%in%c("HALLMARK_GLYCOLYSIS"), ]
matExprs2["GO_HEXOSE_METABOLIC_PROCESS", ] <- 0
matExprs2["GO_HEXOSE_METABOLIC_PROCESS", CombinedConcord] <- 1
matExprs2["GO_HEXOSE_METABOLIC_PROCESS", CombinedDisConcord] <- 2

CombinedConcord <- names(which(colSums(matExprs2[c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "GO_OXIDATIVE_PHOSPHORYLATION","Oxidative phosphorylation"), MTC]) > 0))
others <- setdiff(colnames(NESpathways_PDCs), MTC)
CombinedDisConcord <- names(which(colSums(matExprs2[c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "GO_OXIDATIVE_PHOSPHORYLATION","Oxidative phosphorylation"), others]) > 0))
matExprs2 <- matExprs2[!rownames(matExprs2)%in%c("HALLMARK_OXIDATIVE_PHOSPHORYLATION","GO_OXIDATIVE_PHOSPHORYLATION"), ]
matExprs2["Oxidative phosphorylation", ] <- 0
matExprs2["Oxidative phosphorylation", CombinedConcord] <- 3
matExprs2["Oxidative phosphorylation", CombinedDisConcord] <- 2

CombinedConcord <- names(which(colSums(matExprs2[c("GO_MITOCHONDRIAL_MEMBRANE_PART", "GO_INNER_MITOCHONDRIAL_MEMBRANE_PROTEIN_COMPLEX"), Mitoch]) > 0))
others <- setdiff(colnames(NESpathways_PDCs), Mitoch)
CombinedDisConcord <- names(which(colSums(matExprs2[c("GO_MITOCHONDRIAL_MEMBRANE_PART", "GO_INNER_MITOCHONDRIAL_MEMBRANE_PROTEIN_COMPLEX"), others]) > 0))
matExprs2 <- matExprs2[!rownames(matExprs2)%in%c("GO_MITOCHONDRIAL_MEMBRANE_PART"), ]
matExprs2["GO_INNER_MITOCHONDRIAL_MEMBRANE_PROTEIN_COMPLEX", ] <- 0
matExprs2["GO_INNER_MITOCHONDRIAL_MEMBRANE_PROTEIN_COMPLEX", CombinedConcord] <- 3
matExprs2["GO_INNER_MITOCHONDRIAL_MEMBRANE_PROTEIN_COMPLEX", CombinedDisConcord] <- 2

CombinedConcord <- names(which(colSums(matExprs2[c("GO_PRESYNAPTIC_MEMBRANE", "GO_PRESYNAPTIC_ACTIVE_ZONE"), NEU]) > 0))
others <- setdiff(colnames(NESpathways_PDCs), NEU)
CombinedDisConcord <- names(which(colSums(matExprs2[c("GO_PRESYNAPTIC_MEMBRANE", "GO_PRESYNAPTIC_ACTIVE_ZONE"), others]) > 0))
matExprs2 <- matExprs2[!rownames(matExprs2)%in%c("GO_PRESYNAPTIC_MEMBRANE"), ]
matExprs2["GO_PRESYNAPTIC_ACTIVE_ZONE", ] <- 0
matExprs2["GO_PRESYNAPTIC_ACTIVE_ZONE", CombinedConcord] <- 4
matExprs2["GO_PRESYNAPTIC_ACTIVE_ZONE", CombinedDisConcord] <- 2

CombinedConcord <- names(which(colSums(matExprs2[c("GO_REGULATION_OF_NEUROTRANSMITTER_LEVELS", "GO_NEUROTRANSMITTER_TRANSPORT"), NEU]) > 0))
others <- setdiff(colnames(NESpathways_PDCs), NEU)
CombinedDisConcord <- names(which(colSums(matExprs2[c("GO_REGULATION_OF_NEUROTRANSMITTER_LEVELS", "GO_NEUROTRANSMITTER_TRANSPORT"), others]) > 0))
matExprs2 <- matExprs2[!rownames(matExprs2)%in%c("GO_REGULATION_OF_NEUROTRANSMITTER_LEVELS"), ]
matExprs2["GO_NEUROTRANSMITTER_TRANSPORT", ] <- 0
matExprs2["GO_NEUROTRANSMITTER_TRANSPORT", CombinedConcord] <- 4
matExprs2["GO_NEUROTRANSMITTER_TRANSPORT", CombinedDisConcord] <- 2

library(gplots)
library(heatmap3)
heatmap3(matExprs2[rev(rownames(matExprs2)), ], 
         ColSideColors = ColSideColors,  
         key = FALSE, 
         Colv = NA,
         Rowv = NA,
         scale = "none", 
         col = colorRampPalette(c("white","red","gainsboro","green3","blue","cyan"))(6),
         dendrogram = "none",
         labRow = NA, labCol = NA,
         side.height.fraction = 0.3, keysize = 0.1, cexRow = 0.8)


#BOTTOM PANEL

#upload the expression matrix relative to the PDCs dataset
load("RData/geData_exprs_PDCs.RData")

#upload the resulted significantly differential expressed genes in the GBM subtypes
listOfDiffs <- vector("list", 4)
names(listOfDiffs) <- c("GPM", "MTC", "NEU", "PPR")
listOfDiffs$GPM <- c("SLC16A3", "MET", "HK1", "H6PD","G6PD", "PFKFB3",  "CD44")
listOfDiffs$MTC <-  c("NDUFS3","NDUFS4","TOMM7", "CYCS", "COX16", "TIMM9" , "SDHAF1")
listOfDiffs$PPR <-  c("EZH2", "CCNA2","MKI67", "RAD51", "MCM6", "CCNB2", "PRC1")
listOfDiffs$NEU <-  c("GRIA2","GRIN2D","SEMA6C","SATB1", "WNT4", "SYN3", "DCX")

#upload the PDCs classification
load("RData/PDCsclassification.RData")

ss <- intersect(colnames(geData_exprs_PDCs), names(PDCsclassification))
PDCsclassification <- PDCsclassification[ss]
PDCsclassification <- PDCsclassification[c(which(PDCsclassification == "GPM"), which(PDCsclassification == "MTC"), 
                                           which(PDCsclassification == "NEU"), which(PDCsclassification == "PPR"))]
geData_exprs_PDCs <- geData_exprs_PDCs[, names(PDCsclassification)]
ss <- intersect(rownames(geData_exprs_PDCs), as.character(unlist(listOfDiffs)))
listOfDiffs <- lapply(listOfDiffs, function(x) x[x%in%ss])
geData_exprs_PDCs <- geData_exprs_PDCs[as.character(unlist(listOfDiffs)), ]

#generate the bottom panel of the ED Fig.10a
toPlot <- geData_exprs_PDCs[, names(PDCsclassification)]
colRow <- rep(0, nrow(toPlot))
names(colRow) <- rownames(toPlot)
colRow[names(colRow)%in%listOfDiffs$GPM] <- "red"
colRow[names(colRow)%in%listOfDiffs$MTC] <- "green3"
colRow[names(colRow)%in%listOfDiffs$PPR] <- "cyan"
colRow[names(colRow)%in%listOfDiffs$NEU] <- "blue"
row_annotation <- cbind(colRow)
colnames(row_annotation) <- NULL
rownames(row_annotation) <- NULL
RowSideColors = row_annotation

colCLUSTERS2 <- PDCsclassification[colnames(toPlot)]
tmp <- c("red","green3", "blue","cyan")
names(tmp) <- c("GPM", "MTC", "NEU", "PPR")
aa <- tmp[colCLUSTERS2]
names(aa) <- names(colCLUSTERS2)
colCLUSTERS2 <- aa
names(colCLUSTERS2) <- colnames(toPlot)
colCLUSTERS2[is.na(colCLUSTERS2)] <- "white"
table(colCLUSTERS2, useNA = "ifany")

annot2 <- cbind(colCLUSTERS2)
colnames(annot2) <- c("Clusters")
ColSideColors = annot2[, "Clusters"]

GPM <- names(PDCsclassification[PDCsclassification == "GPM"])
MTC <- names(PDCsclassification[PDCsclassification == "MTC"])
NEU <- names(PDCsclassification[PDCsclassification == "NEU"])
PPR <- names(PDCsclassification[PDCsclassification == "PPR"])

annotGenes <- c(rep("GPM", length(listOfDiffs$GPM)), rep("MTC", length(listOfDiffs$MTC)), rep("NEU", length(listOfDiffs$NEU)), rep("PPR", length(listOfDiffs$PPR)))
names(annotGenes) <- as.character(unlist(listOfDiffs))

matExprs <- matrix(0, nrow=length(as.character(unlist(listOfDiffs))), ncol=ncol(toPlot))
rownames(matExprs) <- as.character(unlist(listOfDiffs))
colnames(matExprs) <- colnames(toPlot)

for (i in 1:nrow(matExprs)){
  gene <- rownames(matExprs)[i]
  xx <- quantile(toPlot[gene, ], probs = .75)
  groupGene <- as.character(annotGenes[gene])
  concordSS <- names(which(toPlot[gene, names(which(PDCsclassification == groupGene))] > xx))
  DisconcordSS <- names(which(toPlot[gene, names(which(PDCsclassification != groupGene))] > xx))
  matExprs[gene, concordSS] <- 1
  matExprs[gene, DisconcordSS] <- 2
}

matExprs[names(which(annotGenes == "MTC")), names(which(PDCsclassification == "MTC"))][matExprs[names(which(annotGenes == "MTC")), names(which(PDCsclassification == "MTC"))] == 1] <- 3
matExprs[names(which(annotGenes == "NEU")), names(which(PDCsclassification == "NEU"))][matExprs[names(which(annotGenes == "NEU")), names(which(PDCsclassification == "NEU"))] == 1] <- 4
matExprs[names(which(annotGenes == "PPR")), names(which(PDCsclassification == "PPR"))][matExprs[names(which(annotGenes == "PPR")), names(which(PDCsclassification == "PPR"))] == 1] <- 5

library(gplots)
library(heatmap3)
heatmap3(matExprs[rev(rownames(matExprs)), ], 
         ColSideColors = ColSideColors,  
         key = FALSE, 
         Colv = NA,
         Rowv = NA,
         scale = "none", 
         col = colorRampPalette(c("white","red","gainsboro","green3","blue","cyan"))(6),
         dendrogram = "none",
         labRow = rownames(matExprs), labCol = NA,
         side.height.fraction = 0.3, keysize = 0.1, cexRow = 1.3)

