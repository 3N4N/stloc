source("./functions.R")

library(SingleCellExperiment)
library(scater)
library(scran)

if (!file.exists("output")) {
  system("mkdir -p output")
}

## Process MOB dataset

counts_raw <- read.delim("./datasets/Rep11_MOB_count_matrix-1.tsv", header = TRUE, row.names = 1)
coords_raw <- do.call(rbind, strsplit(rownames(counts_raw), "x"))
coords <- apply(coords_raw, 1:2, as.numeric)
colnames(coords) <- c("x","y")
rownames(coords) <- rownames(counts_raw)
counts <- t(counts_raw)

## Determine highly variable genes

# This process is taken from https://github.com/MarioniLab/scHOT2019
sce = SingleCellExperiment(assays = list(counts = counts), colData = coords)
sce <- scater::normalize(sce)
var.fit <- trendVar(sce, parametric=TRUE, loess.args=list(span=0.3), use.spikes = FALSE)
var.out <- decomposeVar(sce, var.fit)
seqvals = seq(min(var.out$mean), max(var.out$mean), length.out = 1000)
peakExp = seqvals[which.max(var.fit$trend(seqvals))]
hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$mean > peakExp),]
hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),]
pdf(file = "./output/HVG_selection.pdf", height = 8, width = 8)
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression",
     ylab="Variance of log-expression")
curve(var.fit$trend(x), col="dodgerblue", lwd=2, add=TRUE)
points(var.out$mean[which(var.out$FDR <= 0.05 & var.out$mean > peakExp)],
       var.out$total[which(var.out$FDR <= 0.05 & var.out$mean > peakExp)], col="red", pch=16)
abline(v = peakExp, lty = 2, col = "black")
dev.off()
HVG = sort(rownames(hvg.out))
length(HVG)

# The process below should work with R v4.0.
# Or at least with R v3.6 with BiocManager v3.10.
# But I'm confused about how it works.

# sce = SingleCellExperiment(assays = list(counts = counts), colData = coords)
# sce <- logNormCounts(sce)
# dec <- modelGeneVar(sce)
# hvg <- getTopHVGs(dec,fdr.threshold=0.05)
# top.hvgs <- getTopHVGs(dec, prop=0.1)
# top.hvgs2 <- getTopHVGs(dec, n=304)
# top.hvgs3 <- getTopHVGs(dec, var.threshold=0)
# top.hvgs4 <- getTopHVGs(dec, fdr.threshold=0.05)
# HVG = sort(hvg)
# HVG1 = sort(top.hvgs)
# HVG2 = sort(top.hvgs2)
# HVG3 = sort(top.hvgs3)
# HVG4 = sort(top.hvgs4)
# length(HVG)
# seqvals = seq(min(dec$mean), max(dec$mean), length.out = 1000)
# peakExp = seqvals[which.max(metadata(dec)$trend(seqvals))]
# pdf(file = "./output/HVG_selection4.pdf", height = 8, width = 8)
# plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
# curve(metadata(dec)$trend(x), col="blue", add=TRUE)
# points(dec$mean[ which(rownames(dec) %in% HVG4 & dec$mean > peakExp)],
#        dec$total[which(rownames(dec) %in% HVG4 & dec$mean > peakExp)], col="red", pch=16)
# abline(v = peakExp, lty = 2, col = "black")
# dev.off()


## Get common genes from MOB's HVG and dataset of marker genes

clusterData <- read.delim("./datasets/mmc2.tsv", header=TRUE)
clusterGenes <- clusterData[,7]
clusterGenes <- unique(clusterGenes)
commonGenes <- intersect(HVG, clusterGenes)
clusterData <- as.data.frame(clusterData)
write.table(clusterData[clusterData$gene %in% commonGenes,], file = "./datasets/mmc2-2.tsv", row.names=FALSE, sep="\t")


## Calculate P-values of multiway weighted correlation of a set of genes

# a02.MyOligo: Cldn11 Cnp Plp1 Gatm Ptgds Mbp
# a03.Neuron01: Gng13 S100a5 Kif5b Calm1
# a05.MicroG1: Cst3 B2m Ier5 Cst3 Cst3 Tmsb4x B2m Actb
# a08.Mf: Tmsb4x Apoe Sepp1 B2m
# a25.Neuron14: Dcx Nrep Tubb2b Crmp1 Dlx2 Bcl11b Meis2 Sp9 Cpne4 Synpr Tuba1a Nsg2 Map1b Nrxn3 Tubb5 Eif4a1


# genes <- c("Cst3", "B2m", "Ier5", "Cst3", "Cst3", "Tmsb4x", "B2m", "Actb")
genes <- c("Dcx", "Nrep", "Tubb2b", "Crmp1", "Dlx2", "Bcl11b", "Meis2", "Sp9", "Cpne4", "Synpr", "Tuba1a", "Nsg2", "Map1b", "Nrxn3", "Tubb5", "Eif4a1")

pairCount <- as.matrix(rbind(counts[genes,]))
rownames(pairCount) <- genes


W <- weightMatrix_nD(coords, span = 0.05)

wcor <- as.matrix(sapply(1:nrow(coords), function(i) corTaylor(pairCount, W[i, ])))

set.seed(500)
nitr = 1000
pwcor <- matrix(nrow = nitr, ncol = nrow(coords))

pwcor <- sapply(1:nitr, function(i) {
  x <- pairCount
  o = sample(1:nrow(coords))
  x <- t(sapply(1:nrow(pairCount), function(j) {
    x[j,] = pairCount[j,o]
  }))

  pwcor[i,] = sapply(1:nrow(W), function(j) corTaylor(x, W[j, ]))
})

pvals <- matrix(nrow = nrow(coords), ncol = 1)
pvals <- as.matrix(sapply(1:nrow(wcor), function(i) {
  pvals[i,] = sum(pwcor[i,] > wcor[i])/nitr
}))



## Plot the P-values

require(ggforce)
require(patchwork)
require(ggpubr)

df_res <- data.frame(x = coords[,"x"],
                     y = coords[,"y"],
                     wcor = wcor,
                     pvals = pvals)

t = theme(legend.key.width = unit(0.5, "inches")) +
  theme(plot.title = element_text(size = 20)) +
  theme(axis.title = element_text(size = 15))

# plot_this <- ggplot(df_res, aes(x = -x, y = y, fill = pvals)) +
#   geom_voronoi_tile(max.radius = 5) +
#   theme_minimal() +
#   theme(panel.grid = element_blank()) +
#   theme(axis.text = element_blank()) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   theme(legend.position = "bottom") +
#   labs(colour = "",fill = "") +
#   labs(colour = "local correlation",fill = "local correlation") +
#   ylab("") +
#   xlab("") +
#   ggtitle("correlation of both genes") +
#   scale_alpha_continuous(range = c(0,0.5)) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(0,1)) +
#   t +
#   coord_fixed() +
#   guides(fill = guide_colourbar(title.position = "top",
#                                 title.hjust = 0.5)) +
#   theme(legend.title=element_text(size=15)) +
#   NULL

plot_wcor <- ggplot(df_res, aes(x = -x, y = y)) +
  geom_point(aes(colour = wcor), size = 5) +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  theme(axis.text = element_blank()) +
  xlab("") +
  ylab("") +
  labs(colour = "") +
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(hjust = 0.5, face = "italic")) +
  scale_color_viridis_c(breaks = c(0,max(df_res$wcor)),
                        limits = c(0,max(df_res$wcor)),
                        labels = c("Low","High")) +
  t +
  coord_fixed() +
  guides(colour = guide_colourbar(title.position = "top",
                                  title.hjust = 0.5)) +
  theme(legend.title=element_text(size=15)) +
  labs(colour = "Weighted Correlation") +
  NULL

plot_pvals <- ggplot(df_res, aes(x = -x, y = y)) +
  # geom_point(aes(colour = pvals), size = 5) +
  geom_point(aes(colour = -log10(pvals)), size = 5) +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  theme(axis.text = element_blank()) +
  xlab("") +
  ylab("") +
  labs(colour = "") +
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(hjust = 0.5, face = "italic")) +
  scale_color_viridis_c(breaks = c(0,max(df_res$pvals)),
                        limits = c(0,max(df_res$pvals)),
                        labels = c("Low","High")) +
  t +
  coord_fixed() +
  guides(colour = guide_colourbar(title.position = "top",
                                  title.hjust = 0.5)) +
  theme(legend.title=element_text(size=15)) +
  labs(colour = "-log(pval)") +
  NULL
