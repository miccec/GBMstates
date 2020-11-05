### Figure 3a

# upload the differential pathways among the all 5032 biological pathways in TCGA subtypes
library(openxlsx)
diff_survAssPathAll <- read.xlsx("tables/Supplementary Table 6.xlsx", sheet = 5)
diff_survAssPathAll <- diff_survAssPathAll[3:nrow(diff_survAssPathAll), ]

diff_survAssPathAll_GPM <- diff_survAssPathAll[, 1:4]
colnames(diff_survAssPathAll_GPM) <- diff_survAssPathAll_GPM[1, ]
diff_survAssPathAll_GPM <- diff_survAssPathAll_GPM[-1, ]
rownames(diff_survAssPathAll_GPM) <- diff_survAssPathAll_GPM$Pathway

diff_survAssPathAll_MTC <- diff_survAssPathAll[, 5:8]
colnames(diff_survAssPathAll_MTC) <- diff_survAssPathAll_MTC[1, ]
diff_survAssPathAll_MTC <- diff_survAssPathAll_MTC[-1, ]
diff_survAssPathAll_MTC <- diff_survAssPathAll_MTC[!is.na(diff_survAssPathAll_MTC$Pathway),]
rownames(diff_survAssPathAll_MTC) <- diff_survAssPathAll_MTC$Pathway

diff_survAssPathAll_PPR <- diff_survAssPathAll[, 13:16]
colnames(diff_survAssPathAll_PPR) <- diff_survAssPathAll_PPR[1, ]
diff_survAssPathAll_PPR <- diff_survAssPathAll_PPR[-1, ]
diff_survAssPathAll_PPR <- diff_survAssPathAll_PPR[!is.na(diff_survAssPathAll_PPR$Pathway),]
rownames(diff_survAssPathAll_PPR) <- diff_survAssPathAll_PPR$Pathway

diff_survAssPathAll_NEU <- diff_survAssPathAll[, 9:12]
colnames(diff_survAssPathAll_NEU) <- diff_survAssPathAll_NEU[1, ]
diff_survAssPathAll_NEU <- diff_survAssPathAll_NEU[-1, ]
diff_survAssPathAll_NEU <- diff_survAssPathAll_NEU[!is.na(diff_survAssPathAll_NEU$Pathway),]
rownames(diff_survAssPathAll_NEU) <- diff_survAssPathAll_NEU$Pathway

# order the differential pathways among the all 5032 biological pathways in each subtype by effect size
diff_survAssPathAll_GPM <- diff_survAssPathAll_GPM[order(diff_survAssPathAll_GPM$`effect size`, decreasing = F),]
diff_survAssPathAll_MTC <- diff_survAssPathAll_MTC[order(diff_survAssPathAll_MTC$`effect size`, decreasing = F),]
diff_survAssPathAll_PPR <- diff_survAssPathAll_PPR[order(diff_survAssPathAll_PPR$`effect size`, decreasing = F),]
diff_survAssPathAll_NEU <- diff_survAssPathAll_NEU[order(diff_survAssPathAll_NEU$`effect size`, decreasing = F),]

# create a list with the differential pathwaysamong the all 5032 biological pathways in each subtype and color annotation
listOfDiffs <- vector("list", 4)
names(listOfDiffs) <- c("GPM", "MTC", "PPR", "NEU")
listOfDiffs$GPM <- rownames(diff_survAssPathAll_GPM)
listOfDiffs$MTC <-  rownames(diff_survAssPathAll_MTC)
listOfDiffs$PPR <-  rownames(diff_survAssPathAll_PPR)
listOfDiffs$NEU <-  rownames(diff_survAssPathAll_NEU)

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

# plot the Fig. 3a
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
         labRow = NA,
         cexCol = 0.5,#3,
         cexRow = 0.3,
         scale = 'none', useRaster = F)


### Figure 3b
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

metabolicClassCol <- c("green3", "red", "blue", "cyan")
names(metabolicClassCol) <- c("MTC", "GPM", "NEU", "PPR")

(xRange <- range(xx))
(yRange <- range(yy))


ggray <- "gray90"

library(gplots)
library(heatmap3)
ccol <- NULL
for(i in 1:nrow(NES)){
  whichOfInterest <- rownames(NES)[i]
  tmp <- NES[whichOfInterest, names(class[class == whichOfInterest])]
  tmp <- (tmp - mean(tmp))/sd(tmp)
  bre <- seq(min(tmp), (max(tmp)+0.001), 0.001)
  bre <- colByValue(tmp, col = colorRampPalette(c(ggray, metabolicClassCol[whichOfInterest]))(length(bre)-1), breaks = bre, cex.axis = 0.8)
  tmp <- as.character(bre)
  names(tmp) <- rownames(bre)
  ccol <- c(ccol, tmp)
}


idx <- names(class[class == "MTC"])
ccol[idx[xx[idx] > 0]] <- ggray
ccol[idx[yy[idx] > 0]] <- ggray

idx <- names(class[class == "GPM"])
ccol[idx[yy[idx] > 0]] <- ggray

idx <- names(class[class == "NEU"])
ccol[idx[xx[idx] < 0]] <- ggray
ccol[idx[yy[idx] < 0]] <- ggray

idx <- names(class[class == "PPR"])
ccol[idx[xx[idx] > 0]] <- ggray
ccol[idx[yy[idx] < 0]] <- ggray

ccol <- ccol[names(xx)]

plot(xx, yy, type = "n", xlab = "", ylab = "", xlim = c(-8, 12), ylim = c(-9, 10), cex.axis = 1.5)
points(xx[ccol == ggray], yy[ccol == ggray], col = ggray, pch = 16, cex = 2)
points(xx[ccol != ggray], yy[ccol != ggray], col = ccol[ccol != ggray], pch = 16, cex = 2)
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


### Figure 3c
# upload the TCGA core classification by consensus clustering
library(openxlsx)
TCGAclassification <- read.xlsx("tables/Supplementary Table 6.xlsx", sheet = 3)
TCGAclassification <- TCGAclassification[2:nrow(TCGAclassification), ]
colnames(TCGAclassification) <- TCGAclassification[1, ]
TCGAclassification <- TCGAclassification[-1, ]
rownames(TCGAclassification) <- TCGAclassification$`TCGA ID`
TCGAclassification2 <- TCGAclassification$`Core classification`
names(TCGAclassification2) <- rownames(TCGAclassification)

groups <- factor(TCGAclassification2)
table(groups)

#upload survival data
load("RData/TCGAGBM_clinicalData.RData")
library(survival)
commonSamples <- intersect(rownames(stData), names(groups))
stData <- stData[commonSamples, ]
groups <- groups[commonSamples]

#plot Fig. 3c
par(mfrow = c(2,2),
    oma = c(5, 4, 0, 0) + 0.1,
    mar = c(2.2, 1.8, 1.5, 0) + 1)
aSurv <- survfit(stData ~ groups)
logrank <- survdiff(stData ~ groups)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groups))-1)), digits = 4))
plot(aSurv, col = c("#8D0000", "forestgreen", "blue2", "cyan2"), lwd = 2, font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=1.2)
axis(1,cex.axis=1.2)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("GPM", "MTC")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("#8D0000", "forestgreen"), lwd = 2, font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=1.2)
axis(1,cex.axis=1.2)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("MTC", "NEU")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("forestgreen", "blue2"), lwd = 2,  font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=1.2)
axis(1,cex.axis=1.2)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)

toTake <- names(groups)[groups %in% c("MTC", "PPR")]
groupsGG <- factor(groups[toTake])
stDataGG <- stData[names(groupsGG), ]
aSurv <- survfit(stDataGG ~ groupsGG)
logrank <- survdiff(stDataGG ~ groupsGG)
(pv <- format(1 - pchisq(logrank$chisq, df = (length(levels(groupsGG))-1)), digits = 4))
plot(aSurv, col = c("forestgreen", "cyan2"), lwd = 2,  font.lab=2, bty = "L", yaxt="n", xaxt="n")
axis(2,cex.axis=1.2)
axis(1,cex.axis=1.2)
points(aSurv$time[aSurv$n.censor == 1], aSurv$surv[aSurv$n.censor == 1], pch = 3)


### Figure 3d
### HR calculated from TCGA activity
load("RData/TCGA_NES_class.RData")

NES_GPM <- NES["GPM", ]
NES_MTC <- NES["MTC", ]
NES_NEU <- NES["NEU", ]
NES_PPR <- NES["PPR", ]

library(survival)
# upload survival data
load("RData/TCGAGBM_clinicalData.RData")

# upload the TCGA core classification by consensus clustering
library(openxlsx)
TCGAclassification <- read.xlsx("tables/Supplementary Table 6.xlsx", sheet = 3)
TCGAclassification <- TCGAclassification[2:nrow(TCGAclassification), ]
colnames(TCGAclassification) <- TCGAclassification[1, ]
TCGAclassification <- TCGAclassification[-1, ]
rownames(TCGAclassification) <- TCGAclassification$`TCGA ID`
TCGAclassification2 <- TCGAclassification$`Core classification`
names(TCGAclassification2) <- rownames(TCGAclassification)
groups <- factor(TCGAclassification2)
table(groups)

ss <- intersect(names(groups), names(stData))
NES_MTC <- NES_MTC[ss]
NES_GPM <- NES_GPM[ss]
NES_NEU <- NES_NEU[ss]
NES_PPR <- NES_PPR[ss]
stData <- stData[ss]
groups <- groups[ss]

NES_MTC <- sort(NES_MTC, decreasing = F)
NES_GPM <- sort(NES_GPM, decreasing = F)
NES_NEU <- sort(NES_NEU, decreasing = F)
NES_PPR <- sort(NES_PPR, decreasing = F)

aCoxMTC <- coxph(stData[names(NES_MTC)] ~ NES_MTC)
aPredMTC <- predict(aCoxMTC, type = "risk", se.fit = T)

xMTC <- 1:length(NES_MTC)
hrMTC <- aPredMTC$fit
highMTC <- hrMTC + 2*aPredMTC$se.fit
lowMTC <- hrMTC - 2*aPredMTC$se.fit

aCoxGPM <- coxph(stData[names(NES_GPM)] ~ NES_GPM)
aPredGPM <- predict(aCoxGPM, type = "risk", se.fit = T)

xGPM <- 1:length(NES_GPM)
hrGPM <- aPredGPM$fit
highGPM <- hrGPM + 2*aPredGPM$se.fit
lowGPM <- hrGPM - 2*aPredGPM$se.fit

aCoxNEU <- coxph(stData[names(NES_NEU)] ~ NES_NEU)
aPredNEU <- predict(aCoxNEU, type = "risk", se.fit = T)

xNEU <- 1:length(NES_NEU)
hrNEU <- aPredNEU$fit
highNEU <- hrNEU + 2*aPredNEU$se.fit
lowNEU <- hrNEU - 2*aPredNEU$se.fit

aCoxPPR <- coxph(stData[names(NES_PPR)] ~ NES_PPR)
aPredPPR <- predict(aCoxPPR, type = "risk", se.fit = T)

xPPR <- 1:length(NES_PPR)
hrPPR <- aPredPPR$fit
highPPR <- hrPPR + 2*aPredPPR$se.fit
lowPPR <- hrPPR - 2*aPredPPR$se.fit

par(mfrow = c(2,2),
    oma = c(5, 4, 0, 0) + 0.1,
    mar = c(2.2, 1.8, 1.5, 0) + 1)
plot(xMTC, hrMTC, type = "n", ylim = c(0.5, 1.65), xlim=c(0.7, 302.5),  xaxs="i", xaxt = "n", cex.axis=1.5)
polygon(c(xMTC, rev(xMTC)), c(lowMTC, rev(highMTC)), col = "#87b9e8",border=NA)
lines(xMTC, hrMTC, lwd = 2)

plot(xGPM, hrGPM, type = "n", ylim = c(0.5, 1.65), xlim=c(0.7, 302.5),  xaxs="i", xaxt = "n", cex.axis=1.5)
polygon(c(xGPM, rev(xGPM)), c(lowGPM, rev(highGPM)), col = "#87b9e8",border=NA)
lines(xGPM, hrGPM, lwd = 2)

plot(xNEU, hrNEU, type = "n", ylim = c(0.5, 1.65), xlim=c(0.7, 302.5),  xaxs="i", xaxt = "n", cex.axis=1.5)
polygon(c(xNEU, rev(xNEU)), c(lowNEU, rev(highNEU)), col = "#87b9e8",border=NA)
lines(xNEU, hrNEU, lwd = 2)

plot(xPPR, hrPPR, type = "n", ylim = c(0.5, 1.65), xlim=c(0.7, 302.5),  xaxs="i", xaxt = "n", cex.axis=1.5)
polygon(c(xPPR, rev(xPPR)), c(lowPPR, rev(highPPR)), col = "#87b9e8",border=NA)
lines(xPPR, hrPPR, lwd = 2)


### Figure 3e
library(openxlsx)
classTab <- read.xlsx("tables/Supplementary Table 10.xlsx", startRow = 3, sheet = 4)

class <- classTab$`Pathway-based.classification`
wangClass <- classTab$Wang.et.al..2017.classification
phillipsClass <- classTab$Phillips.et.al..2006.classification

transMat <- as.data.frame(table(class, wangClass), stringsAsFactors = F)
transMat$class <- factor(transMat$class, levels = c("GPM", "MTC", "NEU", "PPR"))
transMat$wangClass <- factor(transMat$wangClass, levels = c("Classical", "Mesenchymal", "Proneural"))

library(flipPlots)
SankeyDiagram(transMat[, -3],
              link.color = "Source", 
              label.show.varname = FALSE,
              font.size = 0,
              weights = transMat$Freq,
              node.width = 100,
              colors = c("red", "green3", "blue", "cyan", "dodgerblue2", "green", "purple")
) 

table(class)/length(class)
table(wangClass)/length(wangClass)

transMat <- as.data.frame(table(class, phillipsClass), stringsAsFactors = F)
transMat$class <- factor(transMat$class, levels = c("GPM", "MTC", "NEU", "PPR"))
transMat$phillipsClass <- factor(transMat$phillipsClass, levels = c("Mesenchymal", "Proliferative", "Proneural"))

library(flipPlots)
SankeyDiagram(transMat[, -3],
              link.color = "Source", 
              label.show.varname = FALSE,
              font.size = 0,
              weights = transMat$Freq,
              node.width = 100,
              colors = c("red", "green3", "blue", "cyan", "red3", "blue4", "forestgreen")
) 

table(class)/length(class)
table(phillipsClass)/length(phillipsClass)


### Figure 3f
source("functions.R")

library(openxlsx)
classTab <- read.xlsx("tables/Supplementary Table 11.xlsx", startRow = 3)

TumorP <- paste0(classTab$Tumor.ID, "_P")
TumorR <- paste0(classTab$Tumor.ID, "_R")
classP <- classTab$`Primary.Tumor.-.GBM.subtype`
classR <- classTab$`Recurrent.Tumor.-.GBM.subtype`
classR[classR == "Unclassified"] <- "nc"

table(classP)
table(classR)

mat <- table(classP, classR)
mat <- as.matrix(as.data.frame.matrix(mat))
rownames(mat) <- paste0(rownames(mat), "_P")
colnames(mat) <- paste0(colnames(mat), "_R")

mat <- mat[, c(1:2, 4:5, 3)]
mat <- rbind(mat, 0)
rownames(mat)[5] <- "nc_P"

library(Gmisc)
library(grid)

transitionPlot2(mmat, box_txt = rep("", 5), tot_spacing = 0.1, box_width = 0.2,
                fill_start_box = c("#E32419", "#45A048", "#005EDB", "#69E3DC", "gray80"), 
                fill_end_box = c("#E32419", "#45A048", "#005EDB", "#69E3DC", "gray80"), 
                type_of_arrow = "gradient",
                lwd_prop_total = F,
                abs_arrow_width = T
)

