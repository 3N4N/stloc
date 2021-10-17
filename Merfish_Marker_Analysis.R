#  ======================================================================
#                                  load libraries
#  ======================================================================

library(SingleCellExperiment)
library(dplyr)
library(Seurat)
library(patchwork)


#  ======================================================================
#                                  Preprocessing
#  ======================================================================

dataset = read.csv(file="./data/moffitt/s7.csv")

counts_raw = dataset[,-c(1:5, 8,9)]
counts = counts_raw[, colSums(counts_raw) != 0]

coords = -1 * counts[,c(1,2)]
colnames(coords) = c("x", "y")
rownames(coords) = apply(coords, 1, function(i) {
   paste0(i[1], "x", i[2])
})

rownames(counts) = rownames(coords)
counts = counts[,-c(1,2)]


sce = SingleCellExperiment(assays = list(counts = t(counts)), colData = coords)
sce = logNormCounts(sce)
counts = logcounts(sce)




plot(counts[,4])

genes = colnames(counts)[-c(1,2,3)]
message("Number of genes: ", length(genes))

clusters = unique(counts[,3])
message("Number of clusters: ", length(clusters))
