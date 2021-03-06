#  ----------------------------------------------------------------------
#  Preprocesses the MERFISH dataset
#  and uses `Celltype_Analysis.R` for eigenvalue analysis
#  ----------------------------------------------------------------------


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
if (!file.exists("output/merfish/combination/scatterPlot")) {
    system("mkdir -p output/merfish/combination/scatterPlot")
}


# read the gene expressions
counts.raw <- read.table("./data/combinationSpatial.csv", header = TRUE, sep = ",", row.names = 1)
# counts.raw <- read.table("./data/FaultSpatial.csv", header = TRUE, sep = ",", row.names = 1)

# counts.raw <- select(counts.raw, -Fos) # Remove Fos column only for bregma
# counts.raw <- read.table("./data/merfish/merfishSpatial.csv", header = TRUE, row.names = 1)


# cellCount of each cell type
cellCount <- counts.raw[(length(counts.raw) - 1):length(counts.raw)]


# get the coordinates of the spots
coords.raw <- do.call(rbind, strsplit(rownames(counts.raw), "x"))
coords <- apply(coords.raw, 1:2, as.numeric)
colnames(coords) <- c("x", "y")
rownames(coords) <- rownames(counts.raw)

counts <- t(counts.raw)

# Log normalize the count matrix
sce <- SingleCellExperiment(assays = list(counts = counts), colData = coords)
sce <- logNormCounts(sce)
counts <- logcounts(sce)
counts <- t(apply(counts, 1, minmax))


# read marker gene info
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


# create pairs of culter name and its marker genes
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
combinationOutputDirectory <- "output/merfish/combination/scatterPlot/"

# call the eigenvalue analysis function
resultFrame <- analyze(clusters.name, clusters.pair, counts, cellCount, faultOutputDirectory)
# corOutput <- cor(resultFrame)
# write.table(corOutput, "svdCor.txt")
# CV <- sapply(resultFrame, function(x) sd(x) / mean(x) * 100)
