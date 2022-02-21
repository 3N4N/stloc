library(SingleCellExperiment)
library(scran)
library(ks)
require(dplyr)


library(ggplot2)
library(gridExtra)


source("./functions.R")


if (!file.exists("output")) {
    system("mkdir -p output")
}

if (!file.exists("output/merfish/plots")) {
    system("mkdir -p output/merfish/plots")
}

if (!file.exists("output/merfish/Bregma/scatterPlot")) {
    system("mkdir -p output/merfish/Bregma/scatterPlot")
}


if (!file.exists("output/merfish/Fault/scatterPlot")) {
    system("mkdir -p output/merfish/Fault/scatterPlot")
}


counts.raw <- read.table("./data/FaultSpatial.csv", header = TRUE, sep = ",", row.names = 1)
# counts.raw <- select(counts.raw, -Fos) # Remove Fos column only for bregma


# counts.raw <- read.table("./data/merfish/merfishSpatial.csv", header = TRUE, row.names = 1)
cellCount <- counts.raw[(length(counts.raw) - 1):length(counts.raw)] # cellCount of each cell type


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
for (row in 1:nrow(clusters.data)) {
    clusters.data[row, "cell_type"] <- gsub(" ", "", clusters.data[row, "cell_type"], fixed = TRUE)
}
write.table(clusters.data, "./data/merfish/markerGene_for_merfish_data.csv", sep = ",", row.names = FALSE, append = FALSE)

clusters.name <- list("Inhibitory", "Excitatory")

# clusters.name <- unique(clusters.name)

# clusters.name <- unique(clusters.data$cell_type)
clusters.name <- sapply(clusters.name, function(i) i <- toString(i))


clusters.pair <- list()
if (length(clusters.name) == 1) {
    clusters.pair[[clusters.name[1]]] <- clusters.data[clusters.data[, 1] == clusters.name[1], 2]
} else {
    clusters.pair <- sapply(clusters.name, function(i) {
        print(clusters.data[clusters.data[, 1] == i, 2])
        clusters.pair[[i]] <- list(clusters.data[clusters.data[, 1] == i, 2])
    })
}


source("Celltype_Analysis.R")

# bregmaOutputDirectory <- "output/merfish/Bregma/scatterPlot/"
faultOutputDirectory <- "output/merfish/Fault/scatterPlot/"
resultFrame <- analyze(clusters.name, clusters.pair, counts, cellCount, faultOutputDirectory)
# corOutput <- cor(resultFrame)
# write.table(corOutput, "svdCor.txt")
CV <- sapply(resultFrame, function(x) sd(x) / mean(x) * 100)