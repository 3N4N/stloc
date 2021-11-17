#  ----------------------------------------------------------------------
#                                  Load packages
#  ----------------------------------------------------------------------

library(SingleCellExperiment)
library(scran)
library(ks)


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
clusters.pair[["MyeloidFibroblast"]] <- unique(unlist((list(clusters.pair[["Myeloid"]], clusters.pair[["Fibroblast"]])) ))
clusters.pair[["EpithelialFibroblast"]] <- unique(unlist((list(clusters.pair[["Epithelial"]], clusters.pair[["Fibroblast"]])) ))
clusters.pair[["MyeloidEpithelial"]] <- unique(unlist((list(clusters.pair[["Myeloid"]], clusters.pair[["Epithelial"]])) ))


#  ----------------------------------------------------------------------
#                                       kde
#  ----------------------------------------------------------------------

counts.flat <- colSums(counts)

data <- rep(names(counts.flat), counts.flat)
data <- do.call(rbind, strsplit(data, "x"))
data <- apply(data, c(1,2), as.numeric)
colnames(data) <- c("x", "y")

hpi <- Hscv(x=data)
kde <- kde(x=data, H=hpi, eval.points=coords)



#  ----------------------------------------------------------------------
#                                     Analysis
#  ----------------------------------------------------------------------


d <- sort (as.numeric (dist (coords )))[1]
W <- weightMatrix.gaussian(coords, l = 1)

set.seed(500)
nitr <- 1e3
for (cluster in clusters.name) {
    # if (!(cluster=="Epithelial" | cluster=="Fibroblast" | cluster=="Myeloid")) next
    if (!(cluster=="Epithelial" | cluster=="Fibroblast" | cluster=="EpithelialFibroblast")) next
    # if (!(cluster=="EpithelialFibroblast")) next

    # if (!(cluster=="MyeloidFibroblast" | cluster=="EpithelialFibroblast" | cluster=="MyeloidEpithelial")) next
    # if (!(cluster=="MyeloidEpithelial" )) next

    genes <- unlist(c(clusters.pair[cluster]))
    genes <- sapply(genes, function(i) i <- toString(i))
    if (length(genes) == 1) next

    message("Analyzing cluster ", cluster)
    message("Number of marker genes: ", length(genes))

    pairCount <- as.matrix(rbind(counts[genes,]))
    rownames(pairCount) <- genes

    # st <- Sys.time()
    # zscores <- apply(pairCount, 1, scale)
    # zscr <- as.matrix(sapply(1:nrow(coords),
    #                          function(i) {
    #                              aggzscores <- sum(zscores[i,])/nrow(pairCount)
    # }))
    # et <- Sys.time()
    # message("Runtime of Z-score: ", difftime(et,st,units="mins"), " mins")

    st <- Sys.time()
    meig.real <- as.matrix(sapply(1:nrow(coords),
            function(i) (maxEigenVal(pairCount, W[i,]) / kde$estimate[i])))
    et <- Sys.time()

    # write.table(meig.real, file=paste0("output/skin_cancer/", cluster, ".txt"), row.names=FALSE, col.names=FALSE)

    message("Runtime of eigenvalues: ", difftime(et,st,units="mins"), " mins")

    message(paste0("Conducting permutation tests for ", cluster))

    meig.perm <- c(1:nitr)
    st <- Sys.time()
    for (i in 1:nitr) {
        cat("\r", "Iteration step", i)
        o <- sample(1:nrow(coords))
        x <- pairCount[,o]
        randloc <- sample(1:nrow(coords), 1)
        meig.perm[i] <- maxEigenVal(x, W[randloc,]) / kde$estimate[o[randloc]]
    }
    cat("\n")
    message(paste0("Permutation tests for ", cluster, " completed"))
    et <- Sys.time()
    message("Runtime of ", cluster, " for ", nitr, " iterations: ", difftime(et,st,units="mins"), " mins")

    meig.pval <- matrix(nrow = nrow(coords), ncol = 1)
    meig.pval <- as.matrix(sapply(1:nrow(meig.real), function(i) {
            meig.pval[i,] <- (sum(meig.perm > meig.real[i]) + 1) / (nitr + 1)
    }))

    # meig.fdr <- p.adjust(meig.pval, method="BH")

    df <- data.frame(x = coords[,"x"], y = -coords[,"y"])

    pdf(paste0("output/skin_cancer/plots/", cluster, "x", log(nitr, 10), ".pdf"),
        height = 6, width = 10, onefile = F)
    plotvals(2, df, "", vals=list(meig.real, -log10(meig.pval)),
                c("Largest Eigenvalue","-log10(pval)"), 3, 1, 2)
    dev.off()
}
