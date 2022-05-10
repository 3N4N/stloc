#  ----------------------------------------------------------------------
#  This file preprocesses the Sking Cancer dataset
#  and uses `Celltype_Analysis.R` for eigenvalue analysis.
#  ----------------------------------------------------------------------



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
counts <- t(apply(counts, 1, minmax))


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

clusters.name <- list("Epithelial", "Fibroblast", "Myeloid")
names(clusters.name) <- clusters.name

## get a list of cluster-gene pairs
clusters.pair <- list()
if (length(clusters.name) == 1) {
    clusters.pair[[clusters.name[1]]] <- clusters.data[clusters.data[,7] == clusters.name[1], 8]
} else {
    clusters.pair <- sapply(clusters.name, function(i) {
        clusters.pair[[i]] <- clusters.data[clusters.data[,7] == i, 8]
  })
}

# clusters.name <- append(clusters.name,"EpithelialFibroblast" ,length(clusters.name))
# clusters.name <- append(clusters.name,"MyeloidFibroblast" ,length(clusters.name))
# clusters.name <- append(clusters.name,"MyeloidEpithelial" ,length(clusters.name))
# clusters.pair[["EpithelialFibroblast"]] <- unique(unlist((list(clusters.pair[["Epithelial"]], clusters.pair[["Fibroblast"]])) ))
# clusters.pair[["MyeloidFibroblast"]] <- unique(unlist((list(clusters.pair[["Myeloid"]], clusters.pair[["Fibroblast"]])) ))
# clusters.pair[["MyeloidEpithelial"]] <- unique(unlist((list(clusters.pair[["Myeloid"]], clusters.pair[["Epithelial"]])) ))



source("Celltype_Analysis.R")
rf <- analyze(clusters.name, clusters.pair, counts, "output/skin_cancer/plots/")
