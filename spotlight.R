library(Matrix)
library(data.table)
library(Seurat)
library(SeuratData)
library(dplyr)
library(gt)
library(SPOTlight)
library(igraph)
library(RColorBrewer)


counts.raw <- read.table("./data/FaultSingleCell.csv", header = TRUE, sep = ",")

# counts.raw <- read.table("./data/merfish/s7.csv", header = TRUE, sep = ",")
data <- subset(counts.raw, select = -c(1, 2, 3, 4, 5, 6, 7, 8, 9))
#remove NA
# data = subset(data, select = -c(Fos))

cortex_sc <- CreateSeuratObject(counts = t(data))

cortex_sc[["subclass"]] <- counts.raw[[8]]


for (i in 1:ncol(data)) {
  print(i)
  data[is.infinite(data[, i]), i] <- mean(data[, i], na.rm = TRUE)
}

for (i in 1:ncol(data)) {
  if (!is.finite(data[, i])) {
    print(data[, i])
  }
}


data <- do.call(data.frame, lapply(DT, function(x) replace(x, is.infinite(x), NA)))



#spatial data
rawSpatial <- read.table("./data/FaultSpatial.csv", header = TRUE, sep = ',')
# spatialData <- rawSpatial[1:(length(rawSpatial)-16)]
spatialData <- rawSpatial[1:(length(rawSpatial)-2)]
# spatialData = subset(spatialData, select = -c(Fos))
anterior = subset(spatialData, select = -c(1))
anterior <- CreateSeuratObject(counts = t(anterior))

#preprocessing

set.seed(123)
cortex_sc <- Seurat::SCTransform(cortex_sc, verbose = FALSE) %>%
  Seurat::RunPCA(., verbose = FALSE) %>%
  Seurat::RunUMAP(., dims = 1:30, verbose = FALSE)

Seurat::DimPlot(cortex_sc,
                group.by = "subclass",
                label = TRUE) + Seurat::NoLegend()

#marker gene
Seurat::Idents(object = cortex_sc) <- cortex_sc@meta.data$subclass
cluster_markers_all <- Seurat::FindAllMarkers(object = cortex_sc, 
                                              assay = "SCT",
                                              slot = "data",
                                              verbose = TRUE, 
                                              only.pos = TRUE)

saveRDS(object = cluster_markers_all,
        file = here::here("inst/markers_sc.RDS"))


set.seed(123)

spotlight_ls <- spotlight_deconvolution(
  se_sc = cortex_sc,
  counts_spatial = anterior@assays$RNA@counts,
  clust_vr = "subclass", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cluster_markers_all, # Dataframe with the marker genes
  cl_n = 100, # number of cells per cell type to use
  hvg = 3000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold
)

saveRDS(object = spotlight_ls, file = here::here("inst/spotlight_ls.rds"))


#deconvolution
spotlight_ls <- readRDS(file = here::here("inst/spotlight_ls.rds"))

nmf_mod <- spotlight_ls[[1]]
decon_mtrx <- spotlight_ls[[2]]
spotlight_correlation = cor(decon_mtrx)


spot_counts <- anterior@assays$RNA@counts


for (r in 1:nrow(decon_mtrx)) {
  for (c in 1:ncol(decon_mtrx)) {
    decon_mtrx1[r, c] <- (decon_mtrx[r, c] / sum(decon_mtrx[r, ])) * 100
  }
}

#cellCount
cellCount = rawSpatial[(length(rawSpatial)-15):length(rawSpatial)]
cellCountCorrelation = cor(cellCount)

cellTypes = colnames(cellCount)

#singular value correlation
singularValueCorrelation <- read.table("./data/merfish/Bregma/singularData/svdCor06.txt", header = TRUE, sep = " ")

fractionCorrelationList = c()
cellCountCorelationList = c()
singularValueCorrelationList = c()

for(r in 1:16){
  for(c in (r + 1):16){
    if(c == 17) break
    cellCountCorelationList = append(cellCountCorelationList, cellCountCorrelation[cellTypes[r], cellTypes[c]])
    fractionCorrelationList = append(fractionCorrelationList, spotlight_correlation[cellTypes[r], cellTypes[c]])
    # singularValueCorrelationList = append(singularValueCorrelationList, singularValueCorrelation[cellTypes[r], cellTypes[c]])
  }
}

cor(fractionCorrelationList, cellCountCorelationList)
cor(singularValueCorrelationList, cellCountCorelationList)


#writing files to 
write.table(cellCountCorrelation, './data/merfish/Bregma/correlations/cellcount/06.txt', append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
write.table(spotlight_correlation, './data/merfish/Bregma/correlations/franction/06.txt', append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
write.table(cellCountCorrelation, './data/merfish/Bregma/correlations/singularval/06.txt', append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)


#scatterplot
#excitatory
scatterDf <- data.frame(x= cellCount[,1], y = decon_mtrx[,1], check.names = FALSE)
scatterPlotSpotlight(scatterDf)
#ependymal
scatterDf <- data.frame(x= cellCount[,16], y = decon_mtrx[,6], check.names = FALSE)
scatterPlotSpotlight(scatterDf)

#inihibitory
scatterDf <- data.frame(x= cellCount[,2], y = decon_mtrx[,8], check.names = FALSE)
scatterPlotSpotlight(scatterDf)

#astrocyte
scatterDf <- data.frame(x= cellCount[,5], y = decon_mtrx[,2], check.names = FALSE)


#to see if 0 cel count has more thena 0 %
positive = 0
negative = 0
for(r in 1: nrow(decon_mtrx)){
  for(celltype in colnames(decon_mtrx)){
    if(cellCount[r, celltype] == 0){
      if(decon_mtrx[r, celltype] < 15)
        positive = positive + 1
      else negative = negative + 1
    } 
  }
}

#cell to percentage
for (r in 1:nrow(cellCount)) {
  for (c in 1:ncol(cellCount)) {
    cellCount1[r, c] <- (cellCount[r, c] / sum(cellCount[r, ])) * 100
  }
}

for(celltype in colnames(cellCount)){
  templist = c()
  for(spot in 1:nrow(cellCount)){
    if(abs(cellCount1[spot, celltype] - decon_mtrx[spot, celltype]) > 5){
      templist = c(templist, cellCount[spot, celltype])
    }
  }

  # templist = cellCount[, celltype]
  paths = "images/spotlight/mispredictions/"
  ttemp = paste(paths, celltype, sep = "")
  pngg = ".png"
  ttemp = paste(ttemp, pngg, sep = "")
  png(ttemp)
  hist(templist, breaks = seq(from=0, to=15, by=1), xlabel=celltype)
  dev.off()  
}

