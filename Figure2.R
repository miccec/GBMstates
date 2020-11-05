### Figure 2a
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


round_percent <- function(x) { 
  x <- x/sum(x)*100
  res <- floor(x)
  rsum <- sum(res)
  if(rsum<100) {
    o <- order(x%%1, sample(length(x)), decreasing=TRUE) 
    res[o[1:(100-rsum)]] <- res[o[1:(100-rsum)]]+1
  } 
  res 
}

percTumor <- table(pData$TumorID, pData$class)
percTumor <- t(apply(percTumor, 1, round_percent))

ccor <- matrix(0, nrow = ncol(percTumor), ncol = ncol(percTumor))
rownames(ccor) <- colnames(percTumor)
colnames(ccor) <- colnames(percTumor)

for(i in 1:ncol(percTumor))
  for(j in 1:ncol(percTumor)){
    tmp <- cor(percTumor[, i], percTumor[, j], method = "spearman")
    ccor[i, j] <- tmp
  }
colnames(ccor) <- rownames(ccor) <- c("GPM", "MTC", "NEU", "PPR")


library(heatmap3)
library(gplots)
heatmap3(ccor, showRowDendro = F, showColDendro = T,
         ColSideColors = c("red", "green3", "blue", "cyan"), ColSideLabs = NA,
         RowSideColors = c("red", "green3", "blue", "cyan"), RowSideLabs = NA,
         col = greenred(75),
         cexCol = 1.5,
         cexRow = 1.5,
         scale = 'none', margins = c(20, 15))


### Figure 2b
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


round_percent <- function(x) { 
  x <- x/sum(x)*100
  res <- floor(x)
  rsum <- sum(res)
  if(rsum<100) {
    o <- order(x%%1, sample(length(x)), decreasing=TRUE) 
    res[o[1:(100-rsum)]] <- res[o[1:(100-rsum)]]+1
  } 
  res 
}

percTumor <- table(pData$TumorID, pData$class)
percTumor <- t(apply(percTumor, 1, round_percent))

library(magrittr)
library(dplyr)

# Cmpute MDS
mds <- percTumor %>%
  dist() %>%          
  cmdscale() %>%
  as_tibble()
colnames(mds) <- c("Dim.1", "Dim.2")

# K-means clustering
clust <- kmeans(mds, 2)$cluster %>%
  as.factor()
levels(clust) <- c("NEU - PPR", "GPM - MTC")
mds <- mds %>%
  mutate(groups = clust)
load("RData/mds.RData")
library(ggpubr)
ggscatter(mds, x = "Dim.1", y = "Dim.2", 
          label = "",
          color = "groups",
          # palette = "jco",
          palette = c("royalblue2", "darkorange2"),
          size = 1.5,
          font.label = 14,
          ellipse = TRUE,
          ellipse.type = "convex",
          repel = TRUE)


### Figure 2c
source("functions.R")

load("RData/Dataset1_class.RData", verbose = T)
load("RData/Dataset1_pData.RData", verbose = T)
pData <- cbind(pData, class = class, stringsAsFactors = F)
dataset1 <- pData
dataset1$TumorID <- paste0(dataset1$TumorID, "_D1")

round_percent <- function(x) { 
  x <- x/sum(x)*100
  res <- floor(x)
  rsum <- sum(res)
  if(rsum<100) { 
    o <- order(x%%1, sample(length(x)), decreasing=TRUE) 
    res[o[1:(100-rsum)]] <- res[o[1:(100-rsum)]]+1
  } 
  res 
}

tmp <- table(dataset1$Sample.type, dataset1$class)
tmp <- t(apply(tmp, 1, round_percent))
tmp

mmat <- diag(0, ncol(tmp))
rownames(mmat) <- colnames(mmat) <- colnames(tmp)

ccond <- which(tmp[1, ] < tmp[2, ])
for(i in ccond){
  mmat[colnames(tmp)[i], colnames(tmp)[i]] <- tmp[1, i]
  tmp[, colnames(tmp)[i]] <- tmp[, colnames(tmp)[i]] - tmp[1, colnames(tmp)[i]]
} 

ccond <- which(tmp[1, ] > tmp[2, ])
for(i in ccond){
  mmat[colnames(tmp)[i], colnames(tmp)[i]] <- tmp[2, i]
  tmp[, colnames(tmp)[i]] <- tmp[, colnames(tmp)[i]] - tmp[2, colnames(tmp)[i]]
}

ccond <- which(tmp[1, ] != 0)
for(i in ccond){
  idx <- which(tmp[2, ] != 0)
  ttmp <- tmp[1, i] * tmp[2, idx]/sum(tmp[2, idx])
  ttmp <- round(ttmp)
  ttmp[length(ttmp)] <- ttmp[length(ttmp)] + (tmp[1, i] - sum(ttmp))
  
  mmat[colnames(tmp)[i], colnames(tmp)[idx]] <- ttmp
  tmp[1, i] <- tmp[1, i] - tmp[1, i]
  tmp[2, idx] <- tmp[2, idx] - ttmp
}
rownames(mmat) <- paste0(rownames(mmat), "_C")
colnames(mmat) <- paste0(colnames(mmat), "_R")


library(Gmisc)
library(grid)

transitionPlot2(mmat, box_txt = rep("", 4), tot_spacing = 0.1, box_width = 0.2,
                fill_start_box = c("#E32419", "#45A048", "#005EDB", "#69E3DC"), 
                fill_end_box = c("#E32419", "#45A048", "#005EDB", "#69E3DC"), 
                type_of_arrow = "gradient",
                lwd_prop_total = F,
                # min_lwd = unit(0.5, "mm"), max_lwd = unit(10, "mm"),
                abs_arrow_width = T
)


### Figure 2d, e (plus ED Figure 2g-i)
#Tu run in a Jupyter environment

#upload sc expression matrix filtered by most variable genes or selected features
adata=st.read(file_name='./sc_dataset1_exprData_mostVarGenes.tsv')

#upload single cells subtype annotation and subtype colors OR single cells type annotation (rim core) and colors
st.add_cell_labels(adata,file_name='./sc_dataset1_TypeGBMsubtypes.tsv')
st.add_cell_colors(adata,file_name='./sc_dataset1_colGBMsubtypes.tsv')

st.add_cell_labels(adata,file_name='./sc_dataset1_TypeGBMsubtypes.tsv')
st.add_cell_colors(adata,file_name='./sc_dataset1_colRimCore.tsv')

#run STREAM
st.dimension_reduction(adata, method='mlle', n_components=2, feature = 'all', nb_pct=0.1, n_jobs=4)
st.seed_elastic_principal_graph(adata)
st.elastic_principal_graph(adata)
st.optimize_branching(adata)
st.extend_elastic_principal_graph(adata)

# figure relative to rim and core OR GBM subtypes (ED Figure 2g-i)
st.stream_plot(adata, root='S0', preference=['S1', 'S2'], fig_size=(8,8), factor_min_win=1, factor_width=3, save_fig=True)

st.stream_plot_gene(adata, root='S0', fig_size=(8,8), genes=['CCNE2','CDK1', 'CDK2', 'EOMES', 'EMX1', 'SSTR2', 'TBR1', 'SATB1', 'LRRC4', 'GABRB3', 'CHRNA4'], save_fig=True)