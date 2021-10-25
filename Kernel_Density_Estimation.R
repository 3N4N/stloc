#  ======================================================================
#                            Kernel Density Estimation
#  ======================================================================


#  ----------------------------------------------------------------------
#                                 loading packages
#  ----------------------------------------------------------------------


# We can use stats::density to determine kde
# but kdensity packages gives us more freedom

library(kdensity)
library(SingleCellExperiment)
library(scran)


#  ----------------------------------------------------------------------
#                             source helper functions
#  ----------------------------------------------------------------------

source("./functions.R")


#  ----------------------------------------------------------------------
#                                   merfish data
#  ----------------------------------------------------------------------

# dataset = read.csv("./data/moffitt/s7.csv")

# counts_raw = dataset[,-c(1:5, 8,9)]
# counts = counts_raw[, colSums(counts_raw) != 0]

# coords = -1 * counts[,c(1,2)]
# colnames(coords) = c("x", "y")
# rownames(coords) = apply(coords, 1, function(i) {
#    paste0(i[1], "x", i[2])
# })

# rownames(counts) = rownames(coords)
# counts = counts[,-c(1,2)]


# clusters = unique(dataset[,8])
# genes = colnames(counts)
# message(paste0("Number of clusters: ", length(clusters)))
# message(paste0("Number of genes: ", length(genes)))


#  ----------------------------------------------------------------------
#                                 skin cancer data
#  ----------------------------------------------------------------------

counts.raw <- read.delim("./data/skin_cancer/GSE144239_ST_P2_S1_counts.tsv", header = T, row.names = 1)

coords.raw <- do.call(rbind, strsplit(rownames(counts.raw), "x"))
coords <- apply(coords.raw, 1:2, as.numeric)
colnames(coords) <- c("x","y")
rownames(coords) <- rownames(counts.raw)

counts <- t(counts.raw)

sce <- SingleCellExperiment(assays = list(counts = counts), colData = coords)
sce <- logNormCounts(sce)
counts <- logcounts(sce)

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

message(paste0("Number of clusters: ", length(clusters.name)))
message(paste0("Number of genes: ", length(clusters.genes)))





#  ----------------------------------------------------------------------
#                                   kde plotting
#  ----------------------------------------------------------------------

cluster = "Epithelial"

counts.flat.all <- colSums(counts)
counts.flat.epi <- colSums(counts[clusters.pair[[cluster]],])

kde.all.norm = kdensity(minmax(counts.flat.all), kernel = "gaussian")
kde.epi.norm = kdensity(minmax(counts.flat.epi), kernel = "gaussian")

kde.all.std = kdensity(scale(counts.flat.all), kernel = "gaussian")
kde.epi.std = kdensity(scale(counts.flat.epi), kernel = "gaussian")


pdf(paste0("output/kde_cancer_minmax.pdf"),
    height = 6, width = 10, onefile = F)
plot(kde.all.norm, col="red", main="KDE with Min-Max Normalization")
lines(kde.epi.norm, col="blue")
legend("topleft", legend=c("All", cluster),
        col=c("red", "blue"), lty=1)
dev.off()

pdf(paste0("output/kde_cancer_scale.pdf"),
    height = 6, width = 10, onefile = F)
plot(kde.all.std, col="red", main="KDE with Z-score Standardization")
lines(kde.epi.std, col="blue")
legend("topleft", legend=c("All", cluster),
        col=c("red", "blue"), lty=1)
dev.off()
