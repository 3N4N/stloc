#  ----------------------------------------------------------------------
#                                  Load packages
#  ----------------------------------------------------------------------

library(SingleCellExperiment)
library(scran)

library(ggplot2)
library(gridExtra)
library(grid)


#  ----------------------------------------------------------------------
#                              Load helper functions
#  ----------------------------------------------------------------------

source("./functions.R")


#  ----------------------------------------------------------------------
#                            Create output directories
#  ----------------------------------------------------------------------

if (!file.exists("output")) {
    system("mkdir -p output")
}

if (!file.exists("output/skin_cancer/plots")) {
    system("mkdir -p output/skin_cancer/plots")
}

system("rm -rf output/skin_cancer/dump")
system("mkdir -p output/skin_cancer/dump")


#  ----------------------------------------------------------------------
#                              Process Cancer dataset
#  ----------------------------------------------------------------------

counts.raw <- read.delim("./data/skin_cancer/GSE144239_ST_P2_S1_counts.tsv", header = TRUE, row.names = 1)

coords.raw <- do.call(rbind, strsplit(rownames(counts.raw), "x"))
coords <- apply(coords.raw, 1:2, as.numeric)
colnames(coords) <- c("x","y")
rownames(coords) <- rownames(counts.raw)

counts <- t(counts.raw)

sce <- SingleCellExperiment(assays = list(counts = counts), colData = coords)
sce <- logNormCounts(sce)
counts <- logcounts(sce)

#  ----------------------------------------------------------------------
#                         Determine highly variable genes
#  ----------------------------------------------------------------------

# sce <- SingleCellExperiment(assays = list(counts = counts), colData = coords)
# sce <- logNormCounts(sce)
# dec <- modelGeneVar(sce)

# hvg <- getTopHVGs(dec,fdr.threshold = 0.05)
# hvg <- sort(hvg)
# length(hvg)

# seqvals <- seq(min(dec$mean), max(dec$mean), length.out = 1000)
# peakExp <- seqvals[which.max(metadata(dec)$trend(seqvals))]

# pdf(file = "./output/skin_cancer/HVG_selection.pdf", height = 8, width = 8)
# plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
# curve(metadata(dec)$trend(x), col = "blue", add = TRUE)
# points(dec$mean[ which(rownames(dec) %in% hvg)],
#        dec$total[which(rownames(dec) %in% hvg)],
#        col = "red", pch = 16)
# abline(v = peakExp, lty = 2, col = "black")
# dev.off()

#  ----------------------------------------------------------------------
#                            Process reference dataset
#  ----------------------------------------------------------------------

## read reference data
clusters.data <- read.delim("./data/skin_cancer/reference_markers_for_NMF.tsv", header = TRUE)
clusters.data <- as.data.frame(clusters.data)
clusters.genes <- unique(clusters.data$gene)

## get genes common in reference and st data
# clusters.genes <- intersect(hvg, clusters.genes)
clusters.genes <- intersect(rownames(counts), clusters.genes)

## remove data of genes not common in reference and st data
clusters.data <- clusters.data[clusters.data$gene %in% clusters.genes,]
clusters.name <- unique(clusters.data$cluster)
clusters.name <- sapply(clusters.name, function(i) i <- toString(i))

## get a list of cluster-gene pairs
clusters.pair <- list()
if (length(clusters.name) == 1) {
    clusters.pair[[clusters.name[1]]] <- clusters.data[clusters.data[,7] == clusters.name[1], 8]
} else {
    clusters.pair <- sapply(clusters.name, function(i) {
        clusters.pair[[i]] <- clusters.data[clusters.data[,7] == i, 8]
  })
}

clusters.name <- append(clusters.name,"EpithelialFibroblast" ,length(clusters.name))
clusters.name <- append(clusters.name,"MyeloidFibroblast" ,length(clusters.name))
clusters.name <- append(clusters.name,"MyeloidEpithelial" ,length(clusters.name))
clusters.pair[["MyeloidFibroblast"]] <- unlist((list(clusters.pair[["Myeloid"]], clusters.pair[["Fibroblast"]])) )
clusters.pair[["EpithelialFibroblast"]] <- unlist((list(clusters.pair[["Epithelial"]], clusters.pair[["Fibroblast"]])) )
clusters.pair[["MyeloidEpithelial"]] <- unlist((list(clusters.pair[["Myeloid"]], clusters.pair[["Epithelial"]])) )



#  ----------------------------------------------------------------------
#                                     Analysis
#  ----------------------------------------------------------------------


d <- sort (as.numeric (dist (coords )))[1]
W <- weightMatrix.gaussian(coords, l = 0.5)

set.seed(500)
for (nitr in c(1e3)) {
# for (nitr in c(1e3, 1e5)) {
    for (cluster in clusters.name) {
        # if (!(cluster=="Epithelial" | cluster=="Fibroblast" | cluster=="Myeloid")) next
        # if (!(cluster=="Fibroblast" | cluster=="Myeloid" | cluster=="MyeloidFibroblast")) next
        if (!(cluster=="Epithelial")) next

        # if (!(cluster=="MyeloidFibroblast" | cluster=="EpithelialFibroblast" | cluster=="MyeloidEpithelial")) next
        # if (!(cluster=="MyeloidEpithelial" )) next

        genes <- unlist(c(clusters.pair[cluster]))
        genes <- sapply(genes, function(i) i <- toString(i))
        if (length(genes) == 1) next

        print(cluster)
        print(genes)

        pairCount <- as.matrix(rbind(counts[genes,]))
        rownames(pairCount) <- genes

        # st <- Sys.time()
        # zscr <- as.matrix(sapply(1:nrow(coords),
        #                          function(i) {
        #                              zscores <- apply(pairCount, 1, scale)
        #                              aggzscores <- sum(zscores[i,])/nrow(pairCount)
        # }))
        # et <- Sys.time()
        # message("Runtime of Z-score: ", difftime(et,st,units="mins"), " mins")

        st <- Sys.time()
        meig.real <- as.matrix(sapply(1:nrow(coords),
                                      function(i) maxEigenVal(pairCount, W[i,])))
        et <- Sys.time()
        # print(summary(meig.real))
        message("Runtime of eigenvalues: ", difftime(et,st,units="mins"), " mins")

        message(paste0("Conducting permutation tests for ", cluster))

        meig.perm <- c(1:nitr)
        st <- Sys.time()
        for (i in 1:nitr) {
            cat("\r", "Iteration step", i)
            o <- sample(1:nrow(coords))
            x <- pairCount
            x <- t(sapply(1:nrow(pairCount), function(j) {
                    x[j,] <- pairCount[j,o]
            }))
            randloc <- sample(1:nrow(coords), 1)
            meig.perm[i] <- maxEigenVal(x, W[randloc,])
        }
        cat("\n")
        message(paste0("Permutation tests for ", cluster, " completed"))
        et <- Sys.time()
        message("Runtime of ", cluster, " for ", nitr, " iterations: ", difftime(et,st,units="mins"), " mins")

        meig.pval <- matrix(nrow = nrow(coords), ncol = 1)
        meig.pval <- as.matrix(sapply(1:nrow(meig.real), function(i) {
                meig.pval[i,] <- (sum(meig.perm > meig.real[i]) + 1) / (nitr + 1)
        }))

        meig.fdr <- p.adjust(meig.pval, method="BH")

        df <- data.frame(x = coords[,"x"], y = coords[,"y"])

        pdf(paste0("output/skin_cancer/plots/", cluster, "x", log(nitr, 10), ".pdf"),
            height = 6, width = 10, onefile = F)
        plotvals(3, df, cluster,
                 vals=list(meig.real, -log10(meig.pval), -log10(meig.fdr)),
                 c("Largest Eigenvalue","-log10(pval)", "-log10(fdr)"), 3, 1, 3)
        dev.off()
    }
}
