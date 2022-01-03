library(SingleCellExperiment)
library(scran)
library(ks)


library(ggplot2)
library(gridExtra)


source("./functions.R")


if (!file.exists("output")) {
    system("mkdir -p output")
}

if (!file.exists("output/merfish/plots")) {
    system("mkdir -p output/merfish/plots")
}


counts.raw <- read.table("./data/merfish/merfishSpatial.csv", header = TRUE, row.names = 1)
cellCount <- counts.raw[(length(counts.raw) - 15):length(counts.raw)]
correlationMatrix <- cor(cellCount)
counts.raw <- counts.raw[1:(length(counts.raw) - 16)]

coords.raw <- do.call(rbind, strsplit(rownames(counts.raw), "x"))
coords <- apply(coords.raw, 1:2, as.numeric)
colnames(coords) <- c("x", "y")
rownames(coords) <- rownames(counts.raw)

counts <- t(counts.raw)

sce <- SingleCellExperiment(assays = list(counts = counts), colData = coords)
sce <- logNormCounts(sce)
counts <- logcounts(sce)
counts <- t(apply(counts, 1, minmax))


clusters.data <- read.table("./data/merfish/markerGene_for_merfish_data.csv", sep = ",", header = TRUE, check.names = FALSE)
clusters.data <- as.data.frame(clusters.data, check.names = FALSE)
clusters.name <- unique(clusters.data$cell_type)
clusters.name <- sapply(clusters.name, function(i) i <- toString(i))

clusters.pair <- list()
if (length(clusters.name) == 1) {
    clusters.pair[[clusters.name[1]]] <- clusters.data[clusters.data[, 1] == clusters.name[1], 2]
} else {
    clusters.pair <- sapply(clusters.name, function(i) {
        clusters.pair[[i]] <- clusters.data[clusters.data[, 1] == i, 2]
    })
}


source("Celltype_Analysis.R")

analyze(clusters.name, clusters.pair, counts, "output/merfish/plots/")