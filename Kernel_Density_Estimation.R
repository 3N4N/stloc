#  ======================================================================
#                            Kernel Density Estimation
#  ======================================================================


#  ----------------------------------------------------------------------
#                                  load packages
#  ----------------------------------------------------------------------

message("Loading packages . . .")

# We can use stats::density to determine kde
# but kdensity packages gives us more freedom
library(kdensity)

# And ks::kde can handle multivariate data
library(ks)

# SCE is for log-normalization
library(SingleCellExperiment)
library(scran)
library(scater)

library(viridis)

message("Packages loaded.")


#  ----------------------------------------------------------------------
#                                   merfish data
#  ----------------------------------------------------------------------

# message("Processing MERFISH data . . .")

# dataset <- read.csv("./data/merfish/s7.csv")

# counts_raw <- dataset[,-c(1:5, 8,9)]
# counts <- counts_raw[, colSums(counts_raw) != 0]

# coords <- -1 * counts[,c(1,2)]
# colnames(coords) <- c("x", "y")
# rownames(coords) <- apply(coords, 1, function(i) {
#    paste0(i[1], "x", i[2])
# })

# rownames(counts) <- rownames(coords)
# counts <- counts[,-c(1,2)]


# clusters <- unique(dataset[,8])
# genes <- colnames(counts)
# message(paste0("Number of clusters: ", length(clusters)))
# message(paste0("Number of genes: ", length(genes)))

# message("MERFISH data processed.")


#  ----------------------------------------------------------------------
#                                 skin cancer data
#  ----------------------------------------------------------------------

message("Processing skin cancer data . . .")

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
clusters.data <- read.delim("./data/skin_cancer/reference_markers_for_NMF.tsv", header=T)
clusters.data <- as.data.frame(clusters.data)
clusters.genes <- unique(clusters.data$gene)

## get genes common in reference and st data
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

message("Skin cancer data processed.")

message(paste0("Number of clusters: ", length(clusters.name)))
message(paste0("Number of genes: ", length(clusters.genes)))





#  ----------------------------------------------------------------------
#                                       KDE
#  ----------------------------------------------------------------------

message("Calculating KDE . . .")

counts.flat <- colSums(counts)

data <- rep(names(counts.flat), counts.flat)
data <- do.call(rbind, strsplit(data, "x"))
data <- apply(data, c(1,2), as.numeric)
colnames(data) <- c("x", "y")

hpi <- Hscv(x=data)
kde <- kde(x=data, H=hpi, eval.points=coords)
kde <- kde(x=data, H=hpi)
# kde <- kde(x = data, gridsize = c(200, 200), xmin = c(-4, -3), xmax = c(4, 3))

message ("KDE calculated. Plotting . . .")

# # plot w/o kde.plot
# image(kde$eval.points[[1]], kde$eval.points[[2]], kde$estimate, col = viridis(20))

# # Perspective plot
# tiff(file="./output/skin_cancer/kde_persp.tif", height=8, width=8, units="in", res=400)
# plot(kde, display = "persp", col.fun = viridis)
# dev.off()

tiff(file="./output/skin_cancer/kde.tif", height=8, width=16, units="in", res=400)
par(mfrow = c(1,2))
plot(kde, display = "slice", cont = c(25, 50, 75), main="Contour Plot")
plot(kde, display = "image", col = viridis(20), main="Colour Plot")
dev.off()

message("Exiting . . .")
