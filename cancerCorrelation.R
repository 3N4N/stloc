source("./functions.R")

library(SingleCellExperiment)
library(scater)
library(scran)

require(ggforce)
require(patchwork)
require(ggpubr)

if (!file.exists("output")) {
  system("mkdir -p output")
}

## Process Cancer dataset

counts_raw <- read.delim("./datasets/skin_cancer/GSE144239_ST_P2_S1_counts.tsv", header = TRUE, row.names = 1)

coords_raw <- do.call(rbind, strsplit(rownames(counts_raw), "x"))
coords <- apply(coords_raw, 1:2, as.numeric)
colnames(coords) <- c("x","y")
rownames(coords) <- rownames(counts_raw)

for (i in 1: ncol(counts_raw)){
  counts_raw[,i]= counts_raw[,i]/max(counts_raw[,i])
}

counts <- t(counts_raw)

### Determine highly variable genes

# sce = SingleCellExperiment(assays = list(counts = counts), colData = coords)
# sce <- logNormCounts(sce)
# dec <- modelGeneVar(sce)

# hvg <- getTopHVGs(dec,fdr.threshold = 0.05)
# hvg <- sort(hvg)
# length(hvg)

# seqvals = seq(min(dec$mean), max(dec$mean), length.out = 1000)
# peakExp = seqvals[which.max(metadata(dec)$trend(seqvals))]

# pdf(file = "./output/HVG_selection.pdf", height = 8, width = 8)
# plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
# curve(metadata(dec)$trend(x), col = "blue", add = TRUE)
# points(dec$mean[ which(rownames(dec) %in% hvg)],
#        dec$total[which(rownames(dec) %in% hvg)],
#        col = "red", pch = 16)
# abline(v = peakExp, lty = 2, col = "black")
# dev.off()


### Select for downstream analysis those marker genes which are also highly variable

clusterData <- read.delim("./datasets/skin_cancer/reference_markers_for_NMF.tsv", header = TRUE)
clusterGenes <- clusterData[,8]
clusterGenes <- unique(clusterGenes)
# commonGenes <- intersect(hvg, clusterGenes)
commonGenes <- intersect(rownames(counts), clusterGenes)
clusterData <- as.data.frame(clusterData)
write.table(clusterData[clusterData$gene %in% commonGenes,],
            file = "./datasets/common.tsv", row.names = FALSE, sep = "\t")


clusterData <- read.delim("./datasets/common.tsv", header = TRUE)
clusterNames <- unique(clusterData[,7])
clusterNames <- sapply(clusterNames, function(i) i <- toString(i))
clusterGenes <- clusterData[,8]
clusterGenePair <- list()
if (length(clusterNames) == 1) {
  clusterGenePair[[clusterNames[1]]] <- clusterData[clusterData[,7] == clusterNames[1], 8]
} else {
  clusterGenePair <- sapply(clusterNames, function(i) {
        clusterGenePair[[i]] <- clusterData[clusterData[,7] == i, 8]
  })
}


if (!file.exists("output/pval_plots_cancer")) {
  system("mkdir output/pval_plots_cancer")
}

plotdf = function(df_res, vals, pvals, valLabel, pvalLabel) {

  # for (i in 1:length(vals)) {
  #   if (pvals[i] > 0.5) {
  #     vals[i] <- NA
  #   }
  # }

  plot_vals <- ggplot(df_res, aes(x = x, y = -y)) +
    geom_point(aes(colour = vals), size = 5) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    theme(axis.text = element_blank()) +
    xlab("") +
    ylab("") +
    labs(colour = "") +
    theme(legend.position = "bottom") +
    theme(plot.title = element_text(hjust = 0.5, face = "italic")) +
    scale_color_viridis_c(na.value = "black", option = "plasma") +
    coord_fixed() +
    guides(colour = guide_colourbar(title.position = "top",
                                    title.hjust = 0.5)) +
    theme(legend.key.width = unit(0.5, "inches")) +
    theme(plot.title = element_text(size = 20)) +
    theme(axis.title = element_text(size = 15)) +
    theme(legend.title = element_text(size = 15)) +
    labs(colour = valLabel) +
    NULL

  plot_pvals <- ggplot(df_res, aes(x = x, y = -y)) +
    geom_point(aes(colour = pvals), size = 5) +
    # geom_point(aes(colour = -log10(pvals)), size = 5) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    theme(axis.text = element_blank()) +
    xlab("") +
    ylab("") +
    labs(colour = "") +
    theme(legend.position = "bottom") +
    theme(plot.title = element_text(hjust = 0.5, face = "italic")) +
    scale_color_viridis_c(option = "plasma") +
    coord_fixed() +
    guides(colour = guide_colourbar(title.position = "top",
                                    title.hjust = 0.5)) +
    theme(legend.key.width = unit(0.5, "inches")) +
    theme(plot.title = element_text(size = 20)) +
    theme(axis.title = element_text(size = 15)) +
    theme(legend.title = element_text(size = 15)) +
    labs(colour = pvalLabel) +
    NULL

  vals_leg = as_ggplot(get_legend(plot_vals))
  pvals_leg = as_ggplot(get_legend(plot_pvals))

  scater::multiplot(plot_vals + theme(legend.position = "none")
                    + theme(plot.margin = margin(10,0,-10,0)),
                    plot_pvals + theme(legend.position = "none")
                    + theme(plot.margin = margin(10,0,-10,0)),
                    vals_leg, pvals_leg,
                    layout = matrix(
                      c(1,1,1,2,2,2,
                        1,1,1,2,2,2,
                        1,1,1,2,2,2,
                        3,3,3,4,4,4), ncol = 6, byrow = TRUE))
}


# W <- weightMatrix_nD(coords, span = 0.3)
W <- weightMatrix_gaussian(coords, l = 0.5)

for (x in clusterNames) {
  if(!(x == "Epithelial" | x == "Fibroblast" | x == "Myeloid")) next
  genes <- unlist(c(clusterGenePair[x]))
  genes <- sapply(genes, function(i) i <- toString(i))
  # genes = sort(genes)
  # genes = genes[1:2]
  if (length(genes) == 1) next

  print(x)
  print(genes)

  pairCount <- as.matrix(rbind(counts[genes,]))
  rownames(pairCount) <- genes

  meanOfZenes = c()
  sdOfZenes = c()
  for(itr in 1:nrow(pairCount)) {
    meanOfZenes = append(meanOfZenes, mean(pairCount[itr,]), length(meanOfZenes))
    sdOfZenes = append(sdOfZenes, sd(pairCount[itr,]), length(sdOfZenes))
  }

  zscr <- as.matrix(sapply(1:nrow(coords),
                           function(i) zScore(pairCount, i, meanOfZenes, sdOfZenes)))
  meig <- as.matrix(sapply(1:nrow(coords),
                           function(i) maxEigenVal(pairCount, W[i,])))

  message(paste0("Conducting permutation tests for ", x))

  # set.seed(500)
  # nitr = 1000

  # pwcor <- matrix(nrow = nitr, ncol = nrow(coords))
  # pwcor <- sapply(1:nitr, function(i) {
  #   x <- pairCount
  #   o = sample(1:nrow(coords))
  #   x <- t(sapply(1:nrow(pairCount), function(j) {
  #       x[j,] = pairCount[j,o]
  #   }))
  #   pwcor[i,] = sapply(1:nrow(W), function(j) corTaylor(x, W[j, ]))
  # })

  # pmeig <- matrix(nrow = nitr, ncol = nrow(coords))
  # pmeig <- sapply(1:nitr, function(i) {
  #   x <- pairCount
  #   o = sample(1:nrow(coords))
  #   x <- t(sapply(1:nrow(pairCount), function(j) {
  #       x[j,] = pairCount[j,o]
  #   }))
  #   pmeig[i,] = sapply(1:nrow(W), function(j) maxEigenVal(x, W[j, ]))
  # })

  # pvals_cor <- matrix(nrow = nrow(coords), ncol = 1)
  # pvals_cor <- as.matrix(sapply(1:nrow(wcor), function(i) {
  #   pvals_cor[i,] = (sum(pwcor[i,] > wcor[i]) + 1) / (nitr + 1)
  # }))
  # pvals_eig <- matrix(nrow = nrow(coords), ncol = 1)
  # pvals_eig <- as.matrix(sapply(1:nrow(meig), function(i) {
  #   pvals_eig[i,] = (sum(pmeig[i,] > meig[i]) + 1) / (nitr + 1)
  # }))


  df_res <- data.frame(x = coords[,"x"],
                       y = coords[,"y"],
                       zscr = zscr,
                       meig = meig)

  pdf(paste0("output/pval_plots_cancer/", x, ".pdf"),
      height = 6, width = 10, onefile = F)
  # plotdf(df_res,df_res$wcor,df_res$pvals_cor,"Weighted Correlation", "-log10(pval)")
  plotdf(df_res,df_res$zscr,df_res$meig,"Z-score", "Largest Eigenvalue")
  dev.off()

  # break
}
