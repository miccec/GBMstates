library(yaGST)
library(doMC)

library(Rtsne)
library(RSpectra)

library(doParallel)
library(heatmap3)
library(cluster)
library(ConsensusClusterPlus)

### single-sample pathway activity
ssMwwGST <- function(geData, geneSet, nCores = 8){
  library(yaGST)
  library(doMC)
  
  means <- rowMeans(geData)
  sds <- apply(geData, 1, sd)
  
  registerDoMC(nCores)
  ans <- foreach(ss = 1:ncol(geData)) %dopar% {
    currentSample <- (geData[, ss] - means)/sds
    rankedList <- sort(currentSample, decreasing = T)
    
    aMwwGST <- lapply(geneSet, function(x) mwwGST(rankedList, geneSet = x, minLenGeneSet = 15, alternative = "two.sided", verbose = F))
    aMwwGST <- aMwwGST[sapply(aMwwGST, length) != 0]
    tmp_NES <- sapply(aMwwGST, function(x) x$log.pu)
    tmp_pValue <- sapply(aMwwGST, function(x) x$p.value)
    
    ans <- list(tmp_NES = tmp_NES, tmp_pValue = tmp_pValue)
    print(ss)
    return(ans)
  }
  NES <- sapply(ans, function(x) x$tmp_NES)
  pValue <- sapply(ans, function(x) x$tmp_pValue)
  colnames(NES) <- colnames(pValue) <- colnames(geData)
  
  FDR <- t(apply(pValue, 1, function(x) p.adjust(x, method = "fdr")))
  res <- list(NES = NES, pValue = pValue, FDR = FDR)
  return(res)
}

### scBiPaD
scBiPaD <- function(dataList, GO, ppath = "scBipaD_output"){
  dir.create(ppath)
  
  for(ddat in 1:length(dataList)){
    geData <- dataList[[ddat]]$geData
    pData <- dataList[[ddat]]$pData
    NES <- dataList[[ddat]]$NES
    
    cs <- Reduce(intersect, list(colnames(geData), rownames(pData), colnames(NES)))
    geData <- geData[, cs]
    pData <- pData[cs, ]
    NNES <- NES[, cs]
    
    pp <- sort(unique(pData$TumorID))
    
    ansMwwGST <- NULL
    
    for(i in 1:length(pp)){
      (tumorID <- pp[i])
      pdi <- pData[pData$TumorID == tumorID, ]
      NESi <- NNES[, rownames(pdi)]
      
      #STEP 1:
      calinsky <- function(hhc, dist = NULL, gMax = round(1 + 3.3 * log(length(hhc$order), 10))) {
        msg <- ""
        if (is.null(dist)) {
          require(clue)
          dist <- sqrt(as.cl_ultrametric(hhc))
          # message(msg <- "The distance matrix not is provided, using the cophenetic matrix")
        } else if (attr(dist, "method") != "euclidean") {
          require(clue)
          dist <- sqrt(as.cl_ultrametric(hhc))
          # message(msg <- "The distance matrix is not euclidean, using the cophenetic matrix")
        }
        dist <- as.matrix(dist)^2
        A <- -dist/2
        A_bar <- apply(A, 1, mean)
        totalSum <- sum(diag(A) - 2 * A_bar + mean(A))
        n <- length(hhc$order)
        ans <- rep(0, gMax)
        for (g in 2:gMax) {
          cclust <- cutree(hhc, k = g)
          withinSum <- 0
          for (k in 1:g) {
            if (sum(cclust == k) == 1) 
              next
            A <- as.matrix(-dist/2)[cclust == k, cclust == k]
            A_bar <- apply(A, 1, mean)
            withinSum <- withinSum + sum(diag(A) - 2 * A_bar + 
                                           mean(A))
          }
          betweenSum <- totalSum - withinSum
          betweenSum <- betweenSum/(g - 1)
          withinSum <- withinSum/(n - g)
          ans[g] <- betweenSum/withinSum
        }
        class(ans) <- "calinsky"
        attr(ans, "message") <- msg
        return(ans)
      }
      
      wData <- t(NESi)
      ddist <- dist(wData, method = "euclidean")
      sHc <- hclust(ddist, method = "ward.D2")
      aCalinsky <- calinsky(sHc, gMax = 5)
      # plot(as.vector(aCalinsky), type = "l", col = "grey", main = tumorID, xlab = "# of groups", ylab = "")
      # text(1:length(aCalinsky), aCalinsky, paste(1:length(aCalinsky)))
      nCls <- which.max(aCalinsky)
      
      library(ConsensusClusterPlus)
      aConsensus <- ConsensusClusterPlus(d = ddist, maxK = ifelse(nCls == 2, nCls + 1, nCls), reps = 10000, pItem = 0.7,
                                         innerLinkage = "ward.D2", finalLinkage = "ward.D2", distance = "euclidean")
      aConsensus <- aConsensus[[nCls]]
      
      obj <- aConsensus
      G <- nCls
      
      consHc <- obj$consensusTree
      consClust <- as.factor(cutree(consHc, k = G))
      names(consClust) <- colnames(NESi)
      levels(consClust) <- palette()[1:G]
      
      ppath_ddat <- paste0(ppath, "/", names(dataList)[ddat])
      dir.create(ppath_ddat)
      # ffile <- paste0(ppath_ddat, "/aConsensus_", tumorID, ".RData")
      # save(aConsensus, consHc, consClust, file = ffile)
      
      col <- colorRampPalette(c("lemonchiffon2", "#053061"))(51)
      
      Colv  <- as.dendrogram(consHc)
      
      consMatrix <- obj$consensusMatrix
      colnames(consMatrix) <- colnames(NESi)
      consMatrix <- consMatrix[, consHc$order]
      
      ColSideColors <- cbind(tumor = pdi$TumorCol)
      
      RowSideColors <- as.character(consClust[colnames(consMatrix)])
      
      # library(heatmap3)
      # ffile <- paste0(ppath_ddat, "/heatmap_", tumorID, ".pdf")
      # pdf(file = ffile, width = 12, height = 8)
      # heatmap3(t(consMatrix), col = col, scale = "none",
      #         Colv = Colv, Rowv = NA,
      #         labRow = NA, labCol = NA,
      #         ColSideColors = ColSideColors, RowSideColors = RowSideColors)
      # dev.off()
      
      #STEP 2:
      groups <- consClust
      levels(groups) <- paste0("c", 1:length(levels(groups)))
      
      consMatrix <- 1 - aConsensus$consensusMatrix
      rownames(consMatrix) <- colnames(consMatrix) <- colnames(NESi)
      consMatrix <- as.dist(consMatrix)
      tmp <- as.numeric(consClust)
      names(tmp) <- names(consClust)
      
      library(cluster)
      aSil <- silhouette(x = tmp, dist = consMatrix)
      # plot(aSil, main = tumorID, col = levels(consClust))
      # abline(v=0.5)
      
      groups <- groups[aSil[, 3] > 0.5]
      
      nGG <- table(groups)
      toTake <- names(which(nGG > 10))
      if(length(toTake) == 0) next
      
      rankedLists <- vector("list", length(toTake))
      names(rankedLists) <- toTake
      for(gg in 1:length(rankedLists)){
        whichOfInterest <- names(groups)[groups == names(rankedLists)[gg]]
        theOthers <- setdiff(names(groups), whichOfInterest)
        
        ans <- apply(geData, 1, function(x) wilcox.test(x[whichOfInterest], x[theOthers]))
        rankedList <- unlist(sapply(ans, function(x) x$statistic))/length(whichOfInterest)/length(theOthers) 
        names(rankedList) <- gsub("\\.W", "", names(rankedList))
        rankedList <- log2(rankedList/(1-rankedList))
        rankedList <- sort(rankedList, decreasing = TRUE)
        
        
        rankedLists[[gg]] <- rankedList
        # print(gg)
      }
      
      library(yaGST)
      aMwwGSTs <- vector("list", length(rankedLists))
      names(aMwwGSTs) <- names(rankedLists)
      for(gg in 1:length(rankedLists)){
        rankedList <- rankedLists[[gg]]
        
        aMwwGST <- lapply(GO, function(x) mwwGST(rankedList, geneSet = x, minLenGeneSet = 15, alternative = "two.sided", verbose = F))
        aMwwGST <- aMwwGST[sapply(aMwwGST, length) != 0]
        
        originalGeneSetCount <- sapply(aMwwGST, function(x) x$originalGeneSetCount)
        actualGeneSetCount <- sapply(aMwwGST, function(x) x$actualGeneSetCount)
        NES <- sapply(aMwwGST, function(x) x$nes)
        odd_NES <- sapply(aMwwGST, function(x) x$pu)
        logit2_NES <- sapply(aMwwGST, function(x) x$log.pu)
        pValue <- sapply(aMwwGST, function(x) x$p.value)
        qValue <- p.adjust(pValue, method = "BH")
        
        aMwwGST <- cbind(originalGeneSetCount, actualGeneSetCount, NES, odd_NES, logit2_NES, pValue, qValue)
        aMwwGST <- as.data.frame(aMwwGST)
        
        aMwwGSTs[[gg]] <- aMwwGST
        # print(gg)
      }
      ffile <- paste0(ppath_ddat, "/rankedLists_aMwwGSTs_", tumorID, ".RData")
      save(rankedLists, aMwwGSTs, file = ffile)
      
      names(aMwwGSTs) <- paste0(tumorID, "_", names(aMwwGSTs))
      ansMwwGST <- c(ansMwwGST, aMwwGSTs)
      
      print(paste0(names(dataList)[ddat], " - ", i, ": TumorID ", tumorID))
    }
  }
  
  #STEP 3:
  rowNames <- sort(unique(unlist(lapply(ansMwwGST, rownames))))
  colNames <- names(ansMwwGST)
  
  NES_allClusters <- matrix(NA, nrow = length(rowNames), ncol = length(colNames))
  rownames(NES_allClusters) <- rowNames
  colnames(NES_allClusters) <- colNames
  pValue_allClusters <- FDR_allClusters <- NES_allClusters
  for(j in 1:length(ansMwwGST)){
    NES_allClusters[rownames(ansMwwGST[[j]]), j] <- ansMwwGST[[j]]$logit2_NES
    pValue_allClusters[rownames(ansMwwGST[[j]]), j] <- ansMwwGST[[j]]$pValue
    FDR_allClusters[rownames(ansMwwGST[[j]]), j] <- ansMwwGST[[j]]$qValue
  }
  
  pValue_allClusters <- pValue_allClusters[rowSums(is.na(NES_allClusters)) == 0, ]
  FDR_allClusters <- FDR_allClusters[rowSums(is.na(NES_allClusters)) == 0, ]
  NES_allClusters <- NES_allClusters[rowSums(is.na(NES_allClusters)) == 0, ]
  
  ffile <- paste0(ppath, "/aMwwGSTs_allClusters.RData")
  save(NES_allClusters, pValue_allClusters, FDR_allClusters, file = ffile)
  
  NES <- NES_allClusters
  NES[NES_allClusters < 0.58] <- 0
  # NES[FDR_allClusters > 0.01] <- 0
  NES[pValue_allClusters > 0.01] <- 0
  NES <- NES[which(rowSums(NES > 0) != 0), ]
  NES <- NES[, which(colSums(NES > 0) != 0)]
  
  jaccard <- matrix(0, nrow = ncol(NES), ncol = ncol(NES))
  rownames(jaccard) <- colnames(jaccard) <- colnames(NES)
  for(i in 1:ncol(NES))
    for(j in 1:ncol(NES)){
      a <- rownames(NES[NES[, i] != 0, i, drop = F])
      b <- rownames(NES[NES[, j] != 0, j, drop = F])
      jaccard[i, j] <- length(intersect(a, b))/length(union(a, b))
      #print(i)
    }
  jaccard <- 1 - jaccard
  
  sHc <- hclust(as.dist(jaccard), method = "average")
  aCalinsky <- calinsky(sHc, gMax = 10)
  nCls <- which.max(aCalinsky)
  
  aConsensus <- ConsensusClusterPlus(d = as.dist(jaccard), maxK = ifelse(nCls == 2, nCls + 1, nCls), reps = 10000, pItem = 0.7,
                                     innerLinkage = "ward.D2", finalLinkage = "ward.D2", distance = "euclidean")
  aConsensus <- aConsensus[[nCls]]
  
  obj <- aConsensus
  G <- nCls
  
  consHc <- obj$consensusTree
  consClust <- as.factor(cutree(consHc, k = G))
  names(consClust) <- colnames(NES)
  levels(consClust) <- palette()[1:G]
  
  col <- colorRampPalette(c("lemonchiffon2", "#053061"))(51)
  
  Colv  <- as.dendrogram(consHc)
  
  consMatrix <- obj$consensusMatrix
  colnames(consMatrix) <- colnames(NES)
  consMatrix <- consMatrix[, consHc$order]
  
  RowSideColors <- as.character(consClust[colnames(consMatrix)])
  
  library(heatmap3)
  llabRow <- colnames(consMatrix)
  ffile <- paste0(ppath, "/heatmap_allClusters.pdf")
  pdf(file = ffile, width = 12, height = 8)
  heatmap3(t(consMatrix), col = col, scale = "none",
           labCol = NA, labRow = llabRow, cexRow = 0.5,
           Colv = Colv, Rowv = NA,
           RowSideColors = RowSideColors)
  dev.off()
  
  ffile <- paste0(ppath, "/aConsensus_allClusters.RData")
  save(aConsensus, consHc, consClust, file = ffile)
}

### Differential expression/pathway analysis more than 2 groups
DEAgroups <- function(ddata, groups, method = c("MWW", "t-test")){
  gg <- sort(unique(groups))
  
  ans <- vector("list", length(gg))
  names(ans) <- gg
  for(i in 1:length(gg)){
    whichOfInterest <- names(groups)[groups == gg[i]]
    theOthers <- setdiff(colnames(ddata), whichOfInterest)
    diffActivity <- apply(ddata, 1, function(x){
      if(method == "MWW") suppressWarnings(a <- wilcox.test(x[whichOfInterest], x[theOthers]))
      if(method == "t-test") suppressWarnings(a <- t.test(x[whichOfInterest], x[theOthers]))
      a <- c(as.numeric(a$statistic), a$p.value)
      return(a)
    })
    diffActivity <- t(diffActivity)
    colnames(diffActivity) <- c("wTest", "pValue")
    fc <- rowMeans(ddata[, whichOfInterest]) - rowMeans(ddata[, theOthers])
    qValue <- p.adjust(diffActivity[, "pValue"], method = "fdr")
    
    diffActivity <- data.frame(statistic = diffActivity[, "wTest"], dm = fc, p.value = diffActivity[, "pValue"], fdr = qValue)
    diffActivity <- diffActivity[, -1]
    colnames(diffActivity) <- c("logFC", "pValue", "qValue")
    # diffActivity <- diffActivity[order(diffActivity$p.value), ]
    # diffActivity <- diffActivity[order(diffActivity$dm, decreasing = T), ]
    ans[[i]] <- diffActivity
  }
  return(ans)
}

### transitionPlot2
transitionPlot2 <- function (transition_flow, type_of_arrow = c("grid", "simple", 
                                                                "gradient"), box_txt = rownames(transition_flow), tot_spacing = 0.2, 
                             box_width = 1/4, fill_start_box = "darkgreen", txt_start_clr = "white", 
                             fill_end_box = fill_start_box, txt_end_clr = txt_start_clr, 
                             cex = 2, min_lwd = if (type_of_arrow == "grid") 1 else unit(0.1, 
                                                                                         "mm"), max_lwd = if (type_of_arrow == "grid") 6 else unit(5, 
                                                                                                                                                   "mm"), lwd_prop_total = TRUE, arrow_clr = "#000000", 
                             abs_arrow_width = FALSE, overlap_bg_clr = "#FFFFFF", overlap_order = 1:nrow(transition_flow), 
                             overlap_add_width = if (type_of_arrow == "grid") 1.5 else unit(1, 
                                                                                            "mm"), box_prop, mar = unit(rep(3, times = 4), "mm"), 
                             main = NULL, box_label = NULL, box_label_pos = "top", box_label_cex = cex, 
                             color_bar = TRUE, color_bar_cex = cex * 0.33, color_bar_labels, 
                             color_bar_subspace = NULL, new_page = FALSE) 
{
  no_boxes <- nrow(transition_flow)
  if (length(dim(transition_flow)) > 2) {
    if (length(dim(transition_flow)) > 3) 
      stop("Your transition matrix should be created through:", 
           " table(var_a, var_b, var_c) providing a 3D-matrix", 
           " you have provided a ", length(dim(transition_flow)), 
           "D matrix.")
    if (!missing(box_prop)) 
      stop("You can't have both box_prop and a three dimensional matrix as input")
    if (dim(transition_flow)[3] != 2) 
      stop("Your third dimension should be a proportion,", 
           " i.e. a variable with two alternatives.", " You have provided a variable with ", 
           dim(transition_flow)[3], " alternatives")
    prop_fn <- function(x) {
      if (x[1] == 0) 
        return(0)
      if (x[2] == 0) 
        return(1)
      return(x[1]/x[2])
    }
    no_1_start <- rowSums(transition_flow[, , 1])
    no_tot_start <- rowSums(transition_flow)
    no_1_end <- colSums(transition_flow[, , 1])
    no_tot_end <- rowSums(colSums(transition_flow[, , 1:2]))
    box_prop <- cbind(apply(cbind(no_1_start, no_tot_start), 
                            1, prop_fn), apply(cbind(no_1_end, no_tot_end), 1, 
                                               prop_fn))
    transition_arrow_props <- transition_flow[, , 1]/(transition_flow[, 
                                                                      , 1] + transition_flow[, , 2])
    if (color_bar == FALSE) {
      color_bar <- "none"
    }
    else if (!is.character(color_bar)) {
      color_bar <- "bottom"
    }
    if (missing(color_bar_labels) && !is.null(dimnames(transition_flow))) {
      color_bar_labels <- dimnames(transition_flow)[[3]]
    }
    transition_flow <- transition_flow[, , 1] + transition_flow[, 
                                                                , 2]
  }
  else if (!missing(box_prop)) {
    transition_arrow_props <- t(sapply(box_prop[, 1], function(x) rep(x, 
                                                                      no_boxes)))
    color_bar <- "none"
  }
  else {
    transition_arrow_props <- matrix(1, ncol = no_boxes, 
                                     nrow = no_boxes)
    color_bar <- "none"
  }
  if (length(arrow_clr) == no_boxes) {
    arrow_clr <- t(sapply(arrow_clr, FUN = function(x) {
      rep(x, ncol(transition_flow))
    }))
  }
  else if (length(arrow_clr) == 1) {
    arrow_clr <- rep(arrow_clr, no_boxes * ncol(transition_flow))
  }
  if (length(arrow_clr) != no_boxes * ncol(transition_flow)) 
    stop("You have provided an invalid number of arrow colors,", 
         " you have ", length(arrow_clr), " colors, while you should provide either 1, ", 
         no_boxes, ", or ", no_boxes * ncol(transition_flow), 
         " colors")
  if (length(overlap_order) != no_boxes) 
    stop("You have the wrong number of overlap orders, you provided ", 
         length(overlap_order), " while it should be ", no_boxes)
  else if (all(overlap_order %in% 1:no_boxes) == FALSE) 
    stop("Your overlap numbers contain numbers outside the rowrange of", 
         " the transition rows, i.e. not between 1 and ", 
         no_boxes)
  type_of_arrow <- match.arg(type_of_arrow)
  if (type_of_arrow != "grid") {
    if (!"unit" %in% class(min_lwd) || !"unit" %in% class(max_lwd)) 
      stop("Your line widths must be in units when you specify the alternative arrows, e.g. unit(10, 'pt')")
    min_lwd <- convertUnit(min_lwd, unitTo = "npc", valueOnly = TRUE)
    max_lwd <- convertUnit(max_lwd, unitTo = "npc", valueOnly = TRUE)
  }
  if (tot_spacing < 0 || tot_spacing > 1) 
    stop("Total spacing, the tot_spacing param,", " must be a fraction between 0-1,", 
         " you provided ", tot_spacing)
  if (box_width < 0 || box_width > 1) 
    stop("Box width, the box_width param,", " must be a fraction between 0-1,", 
         " you provided ", box_width)
  if (is.null(box_txt)) 
    box_txt = matrix("", ncol = 2, nrow = no_boxes)
  if (is.null(dim(box_txt)) && is.vector(box_txt)) 
    if (length(box_txt) != no_boxes) 
      stop("You have an invalid length of text description, the box_txt param,", 
           " it should have the same length as the boxes, ", 
           no_boxes, ",", " but you provided a length of ", 
           length(box_txt))
  else box_txt <- cbind(box_txt, box_txt)
  else if (nrow(box_txt) != no_boxes || ncol(box_txt) != 2) 
    stop("Your box text matrix doesn't have the right dimension, ", 
         no_boxes, " x 2, it has: ", paste(dim(box_txt), collapse = " x "))
  if (missing(box_prop)) {
    fill_start_box <- rep(fill_start_box, length.out = no_boxes)
    txt_start_clr <- rep(txt_start_clr, length.out = no_boxes)
    fill_end_box <- rep(fill_end_box, length.out = no_boxes)
    txt_end_clr <- rep(txt_end_clr, length.out = no_boxes)
  }
  else {
    fill_start_box <- prTpGetBoxPropClr(fill_start_box, no_boxes = no_boxes)
    fill_end_box <- prTpGetBoxPropClr(fill_end_box, no_boxes = no_boxes)
    txt_start_clr <- prTpGetBoxPropClr(txt_start_clr, no_boxes = no_boxes, 
                                       lengthOneOK = TRUE)
    txt_end_clr <- prTpGetBoxPropClr(txt_end_clr, no_boxes = no_boxes, 
                                     lengthOneOK = TRUE)
    if (is.matrix(box_prop) == FALSE) 
      stop("You have to provide the box_prop as a matrix corresponding to the boxes")
    else if (nrow(box_prop) != no_boxes || ncol(box_prop) != 
             2) 
      stop("Your box_prop matrix must have ", no_boxes, 
           "x", 2, " dimensions, your matrix is currently of ", 
           nrow(box_prop), "x", ncol(box_prop), " dimensions")
    else if (any(box_prop > 1 | box_prop < 0)) 
      stop("You have provided in box_prop invalid quantiles outside the 0-1 range")
    else if (length(fill_start_box) == 0) 
      stop("You have provided invalid number of fill colors (fill_start_box) when used together with box_prop")
    else if (length(fill_end_box) == 0) 
      stop("You have provided invalid number of fill colors (fill_end_box) when used together with box_prop")
    else if (length(txt_start_clr) == 0) 
      stop("You have provided invalid number of text colors (txt_start_clr) when used together with box_prop")
    else if (length(txt_end_clr) == 0) 
      stop("You have provided invalid number of text colors (txt_end_clr) when used together with box_prop")
  }
  if (nrow(transition_flow) != ncol(transition_flow)) 
    stop("Invalid input array, the matrix is not square but ", 
         nrow(transition_flow), " x ", ncol(transition_flow))
  prop_start_sizes <- rowSums(transition_flow)/sum(transition_flow)
  prop_end_sizes <- colSums(transition_flow)/sum(transition_flow)
  if (sum(prop_end_sizes) == 0) 
    stop("You can't have all empty boxes after the transition")
  if (new_page) 
    grid.newpage()
  vp_depth = 1
  Gmisc:::prPushMarginViewport(bottom = convertY(mar[1], unitTo = "npc"), 
                               left = convertX(mar[2], unitTo = "npc"), top = convertY(mar[3], 
                                                                                       unitTo = "npc"), right = convertX(mar[4], unitTo = "npc"), 
                               "main_margins")
  if (!is.null(main) && nchar(main) > 0) {
    prGridPlotTitle(main, cex[1])
    vp_depth %<>% +2
  }
  if (!is.null(box_label) && length(box_label) == 2) {
    left <- prTpGetBoxPositions(side = "left", no = 1, transitions = transition_flow[1, 
    ], prop_start_sizes = prop_start_sizes, prop_end_sizes = prop_end_sizes, 
    tot_spacing = tot_spacing, box_width = box_width)
    right <- prTpGetBoxPositions(side = "right", no = 1, 
                                 transitions = transition_flow[, 1], prop_start_sizes = prop_start_sizes, 
                                 prop_end_sizes = prop_end_sizes, tot_spacing = tot_spacing, 
                                 box_width = box_width)
    left_label <- textGrob(box_label[1], gp = gpar(cex = box_label_cex))
    right_label <- textGrob(box_label[2], gp = gpar(cex = box_label_cex))
    label_height <- convertY(max(grobHeight(left_label), 
                                 grobHeight(right_label)), unitTo = "npc", valueOnly = TRUE)
    label_height <- unit(label_height * 2 + label_height * 
                           0.1, "npc")
    width <- list(left = unit(left$right - left$left, "npc"), 
                  right = unit(right$right - right$left, "npc"))
    if (box_label_pos == "top") {
      gl <- grid.layout(nrow = 2, ncol = 3, heights = unit.c(label_height, 
                                                             unit(1, "npc") - label_height), widths = unit.c(width$left, 
                                                                                                             unit(1, "npc") - width$left - width$right, width$right))
      label_row_no <- 1
      main_row_no <- 2
    }
    else {
      gl <- grid.layout(nrow = 2, ncol = 3, heights = unit.c(unit(1, 
                                                                  "npc") - label_height, label_height), widths = unit.c(width$left, 
                                                                                                                        unit(1, "npc") - width$left - width$right, width$right))
      label_row_no <- 2
      main_row_no <- 1
    }
    pushViewport(viewport(layout = gl, name = "Label_layout"))
    pushViewport(viewport(layout.pos.row = label_row_no, 
                          layout.pos.col = 1, name = "Left_label"))
    grid.draw(left_label)
    popViewport()
    pushViewport(viewport(layout.pos.row = label_row_no, 
                          layout.pos.col = 3, name = "Right_label"))
    grid.draw(right_label)
    popViewport()
    pushViewport(viewport(layout.pos.row = main_row_no, layout.pos.col = 1:3, 
                          name = "Main_exc_label"))
    vp_depth %<>% +2
  }
  if (color_bar != "none" && type_of_arrow == "gradient") {
    if (color_bar == "bottom") {
      bar_height <- unit(0.05, "npc")
      colorAxis <- xaxisGrob(at = c(0, 0.25, 0.5, 0.75, 
                                    1), label = sprintf("%d %%", c(0, 0.25, 0.5, 
                                                                   0.75, 1) * 100), main = FALSE, gp = gpar(cex = color_bar_cex))
      axis_height <- grobHeight(colorAxis) + unit(0.01, 
                                                  "npc")
      bar_layout <- grid.layout(nrow = 3, ncol = 3, heights = unit.c(unit(1, 
                                                                          "npc") - axis_height - bar_height, axis_height, 
                                                                     bar_height), widths = unit.c(unit(box_width, 
                                                                                                       "npc"), unit(1, "npc") - unit(box_width * 2, 
                                                                                                                                     "npc"), unit(box_width, "npc")))
      pushViewport(viewport(layout = bar_layout, name = "Bar_layout"))
      pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 2, 
                            name = "Color_bar"))
      bar_clrs <- prTpGetColors(fill_start_box[1, ], space = color_bar_subspace)
      grid.raster(t(as.raster(bar_clrs)), width = 1, height = 1, 
                  interpolate = FALSE)
      grid.draw(colorAxis)
      if (!missing(color_bar_labels)) {
        lab_height <- convertY(grobHeight(textGrob("Ij")), 
                               "npc", valueOnly = TRUE)
        lab_cex_adjusted <- 1/(lab_height * 2)
        if (missing(txt_start_clr)) {
          color_bar_txt_clr <- c("black", "black")
        }
        else if (ncol(txt_start_clr) == 1) {
          color_bar_txt_clr <- rep(txt_start_clr[1], 
                                   2)
        }
        else {
          color_bar_txt_clr <- txt_start_clr[1, ]
        }
        lab_margin <- 0.05
        left <- textGrob(color_bar_labels[1], x = 0 + 
                           lab_margin, just = "left", y = 0.5, gp = gpar(cex = lab_cex_adjusted, 
                                                                         col = color_bar_txt_clr[1]))
        right <- textGrob(color_bar_labels[2], x = 1 - 
                            lab_margin, just = "right", y = 0.5, gp = gpar(cex = lab_cex_adjusted, 
                                                                           col = color_bar_txt_clr[2]))
        grid.draw(left)
        grid.draw(right)
      }
      popViewport()
      pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:3, 
                            name = "Main_exc_bar"))
      vp_depth %<>% +2
    }
    else {
      stop("The color bar position you want, '", color_bar, 
           "', is not yet supported")
    }
  }
  shift <- box_width * 0.02
  vp1 <- viewport(x = 0.5 + shift, y = 0.5 - shift, height = 1 - 
                    shift * 2, width = 1 - shift * 2, name = "shadow_boxes")
  pushViewport(vp1)
  shadow_clr <- rep(grey(0.8), length.out = no_boxes)
  Gmisc:::prTpPlotBoxes(overlap_order = overlap_order, transition_flow = transition_flow, 
                        no_boxes = no_boxes, box_width = box_width, tot_spacing = tot_spacing, 
                        txt = matrix("", nrow = no_boxes, ncol = 2), cex = cex, 
                        prop_start_sizes = prop_start_sizes, prop_end_sizes = prop_end_sizes, 
                        box_prop = box_prop, lwd_prop_total = lwd_prop_total, 
                        fill_start_clr = shadow_clr, fill_end_clr = shadow_clr, 
                        txt_start_clr = txt_start_clr, txt_end_clr = txt_end_clr, 
                        line_col = shadow_clr[1],
                        plot_arrows = FALSE, proportion = FALSE)
  popViewport()
  vp1 <- viewport(x = 0.5 - shift, y = 0.5 + shift, height = 1 - 
                    shift * 2, width = 1 - shift * 2, name = "actual_boxes")
  pushViewport(vp1)
  Gmisc:::prTpPlotBoxes(overlap_order = overlap_order, transition_flow = transition_flow, 
                        no_boxes = no_boxes, box_width = box_width, tot_spacing = tot_spacing, 
                        txt = box_txt, cex = cex, prop_start_sizes = prop_start_sizes, 
                        prop_end_sizes = prop_end_sizes, box_prop = box_prop, 
                        lwd_prop_total = lwd_prop_total, fill_start_clr = fill_start_box, 
                        line_col = "gray80", #luc
                        fill_end_clr = fill_end_box, txt_start_clr = txt_start_clr, 
                        txt_end_clr = txt_end_clr, min_lwd = min_lwd, max_lwd = max_lwd, 
                        overlap_add_width = overlap_add_width, overlap_bg_clr = overlap_bg_clr, 
                        type_of_arrow = type_of_arrow, abs_arrow_width = abs_arrow_width, 
                        arrow_clr = arrow_clr, transition_arrow_props = transition_arrow_props, 
                        color_bar_subspace = color_bar_subspace, plot_arrows = TRUE, 
                        proportion = TRUE)
  popViewport()
  popViewport(vp_depth)
}


#starburstPlot
starburstPlot <- function(DEGS, cut.off.p = 0.01, cut.off.fc.ge = 1, cut.off.fc.meth = 1, porq = "p"){
  DEGs <- DEGS
  DEGS[, grep("Value", colnames(DEGS))] <- -log10(DEGS[, grep("Value", colnames(DEGS))])
  DEGS[DEGS$ge_fc < 0, grep("ge_[p,q]Value", colnames(DEGS))] <- - DEGS[DEGS$ge_fc < 0, grep("ge_[p,q]Value", colnames(DEGS))]
  DEGS[DEGS$md_fc < 0, grep("md_[p,q]Value", colnames(DEGS))] <- - DEGS[DEGS$md_fc < 0, grep("md_[p,q]Value", colnames(DEGS))]
  
  cut.off.p <- -log10(cut.off.p)
  
  ccol <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "red", "#0072B2", "#D55E00", "#CC79A7", "purple")
  names(ccol) <- c("Not Significant", "Up regulated & Hypo methylated", "Down regulated & Hypo methylated",
                   "Hypo methylated", "Hyper methylated", "Up regulated", "Down regulated",
                   "Up regulated & Hyper methylated", "Down regulated & Hyper methylated")
  
  if(porq == "p"){
    x <- DEGS$md_pValue
    y <- DEGS$ge_pValue
  }
  if(porq == "q"){
    x <- DEGS$md_qValue
    y <- DEGS$ge_qValue
  }
  xx <- DEGS$md_fc
  yy <- DEGS$ge_fc
  
  ttype <- rep("Not Significant", length(x))
  ttype[x < -cut.off.p] <- "Hypo methylated"
  ttype[x > cut.off.p] <- "Hyper methylated"
  ttype[y < -cut.off.p] <- "Down regulated"
  ttype[y > cut.off.p] <- "Up regulated"
  ttype[x < -cut.off.p & y > cut.off.p] <- "Up regulated & Hypo methylated"
  ttype[x < -cut.off.p & y < -cut.off.p] <- "Down regulated & Hypo methylated"
  ttype[x > cut.off.p & y > cut.off.p] <- "Up regulated & Hyper methylated"
  ttype[x > cut.off.p & y < -cut.off.p] <- "Down regulated & Hyper methylated"
  ccircle <- rep("Not Significant", length(x))
  ccircle[xx <  -cut.off.fc.meth & yy > cut.off.fc.ge & ttype == "Up regulated & Hypo methylated"] <- "Significant"
  ccircle[xx <  -cut.off.fc.meth & yy < -cut.off.fc.ge & ttype == "Down regulated & Hypo methylated"] <- "Significant"
  ccircle[xx >  cut.off.fc.meth & yy > cut.off.fc.ge & ttype == "Up regulated & Hyper methylated"] <- "Significant"
  ccircle[xx >  cut.off.fc.meth & yy < -cut.off.fc.ge & ttype == "Down regulated & Hyper methylated"] <- "Significant"
  
  mmain <- "Starburst plot"
  if(porq == "p"){
    xLab <- "DNA Methylation\nlog10 pValue"
    yLab <- "Gene Expression\nlog10 pValue"
  } 
  if(porq == "q"){
    xLab <- "DNA Methylation\nlog10 qValue"
    yLab <- "Gene Expression\nlog10 qValue"
  }
  
  ccol <- ccol[ttype]
  ccol[ccircle != "Significant"] <- "gray80"
  ccol[ccol %in% c("#56B4E9", "#CC79A7")] <- "gray80"
  
  par(mar = c(5, 5, 4, 2) + 0.1)
  # plot(x = x, y = y, main = mmain, xlab = xLab, ylab = yLab, type = "n")
  plot(x = x, y = y, main = "", xlab = "", ylab = "", type = "n", xlim = c(-7, 7), ylim = c(-30, 30), cex = 3, cex.axis = 2)
  points(x = x[ccol == "gray80"], y = y[ccol == "gray80"], col = ccol[ccol == "gray80"], cex = 3)
  points(x = x[ccol != "gray80"], y = y[ccol != "gray80"], col = ccol[ccol != "gray80"], cex = 3)
  abline(h = cut.off.p, lty = 2)
  abline(h = -cut.off.p, lty = 2)
  abline(v = cut.off.p, lty = 2)
  abline(v = -cut.off.p, lty = 2)
  # points(x = x[ccircle == "Significant"], y = y[ccircle == "Significant"], pch = 1, cex = 2)
  
  DEGs <- cbind(DEGS, Type = ttype, Significant = ccircle, stringsAsFactors = F)
  
  return(DEGs)
}


#heatmap.3
heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      NumColSideColors = 1,
                      NumRowSideColors = 1,
                      KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*NumColSideColors, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*NumRowSideColors, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(colnames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), colnames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      #par(mar = c(0.5, 0, 0, 40))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      #par(mar = c(0.5, 0, 0, margins[2]))
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), cex.axis=0.7, 
             colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}
