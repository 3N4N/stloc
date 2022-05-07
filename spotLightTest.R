library(ggplot2)
library(SPOTlight)
library(SingleCellExperiment)
# library(SpatialExperiment)
library(scater)
library(scran)
library(NMF)
library(nnls)


for (f in list.files("R/")) {
    source(paste0("R/", f))
}

sce <- mockSC(ng = 200, nc = 10, nt = 3)
spe <- mockSP(sce)
mgs <- getMGS(sce)

singleCellFile <- read.table("./data/twoType/allSingleCell/InhibitoryExcitatorySingle.csv", sep = ",", header = TRUE, check.names = FALSE)
coldata <- DataFrame(
    type = singleCellFile[, "type"]
)
singleCellFile <- subset(singleCellFile, select = -c(type))
singleCellFile <- t(singleCellFile)
singleCellData <- SingleCellExperiment(list(counts = singleCellFile))
colData(singleCellData) <- coldata

markers <- getMGS(singleCellData)
spatialCellFile <- read.table("./data/twoType/spatialData/InhibitoryExcitatorySpatial.csv", sep = ",", header = TRUE, check.names = FALSE)

coord <- spatialCellFile[, "coord"]

metaData <- spatialCellFile[(length(spatialCellFile) - 3):length(spatialCellFile)]
spatialCellFile <- subset(spatialCellFile, select = -c(coord))
spatialCellFile <- subset(spatialCellFile, select = -c((length(spatialCellFile) - 3):length(spatialCellFile)))
spatialCellFile <- t(spatialCellFile)
spatialData <- SingleCellExperiment((list(counts = spatialCellFile)))
colnames(spatialData) <- coord


singleCellData <- logNormCounts(singleCellData)
dec <- modelGeneVar(singleCellData)
hvg <- getTopHVGs(dec, n = 145)


spatialData <- logNormCounts(spatialData)

res <- SPOTlight(
    x = counts(singleCellData),
    y = counts(spatialData),
    groups = as.character(singleCellData$type),
    mgs = markers,
    hvg = hvg,
    weight_id = "weight",
    group_id = "type",
    gene_id = "gene"
)
result <- res$mat
faultVector <- c()
for (i in 1:nrow(result)) {
    if (abs(result[i, 1] - 0.5) > 0.1) {
        faultVector <- append(faultVector, i)
    }
}
rownames <- row.names(result)