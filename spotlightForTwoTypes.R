library(Matrix)
library(data.table)
library(Seurat)
library(SeuratData)
library(dplyr)
library(gt)
library(SPOTlight)
library(igraph)
library(RColorBrewer)
source('funcions.R')

#read all file names from data/twoType/singleCell/
singleCellNames <- list.files(path = "data/twoType/singleCell/", pattern = "*.csv")
#read all file names from data/twoType/spatialData/
spatialNames <- list.files(path = "data/twoType/spatialData/", pattern = "*.csv")

#for loop iterate singleCellNames
for(i in 1:length(singleCellNames)){
  print(i)
  #read in data
  counts.raw <- read.table(paste0("data/twoType/singleCell/", singleCellNames[i]), header = TRUE, sep = ",")
  #remove NA
  data <- subset(counts.raw, select = -c(1, 2, 3, 4, 5, 6, 7, 8, 9))
  #remove NA
  # data = subset(data, select = -c(Fos))
  # data = subset(data, select = -c(Blank_1, Blank_2, Blank_3, Blank_4, Blank_5))
  #create Seurat object
  cortex_sc <- CreateSeuratObject(counts = t(data))
  cortex_sc[["subclass"]] <- counts.raw[[8]]
  #normalize data
  cortex_sc <- NormalizeData(cortex_sc)
  #variable gene finding
  cortex_sc <- FindVariableFeatures(cortex_sc, selection.method = "vst", nfeatures = 2000)
  #save object
  # save(cortex_sc, file = paste0("data/twoType/singleCell/", singleCellNames[i], ".RData"))


  #read in spatial data
  rawSpatial <- read.table(paste0("data/twoType/spatialData/", spatialNames[i]), header = TRUE, sep = ',')

  spatialData <- rawSpatial[1:(length(rawSpatial)-2)]

  anterior = subset(spatialData, select = -c(1))
  anterior <- CreateSeuratObject(counts = t(anterior))
  anterior <- NormalizeData(anterior)
  #variable gene finding
  anterior <- FindVariableFeatures(anterior, selection.method = "vst", nfeatures = 2000)

  # write a try catch
  possibleError = tryCatch(
    {
      set.seed(123)
      cortex_sc <- Seurat::SCTransform(cortex_sc, verbose = FALSE) %>%
        Seurat::RunPCA(., verbose = FALSE) %>%
        Seurat::RunUMAP(., dims = 1:30, verbose = FALSE)
      # save(anterior, file = paste0("data/twoType/spatialData/", spatialNames[i], ".RData"))
    },
    error = function(e)
    {
      # print(e)
    }

  )

  #check if error
  if(is.null(possibleError)){
    #save object
    print('-----------------------------------------------------')
    next
  }



  # set.seed(123)
  # cortex_sc <- Seurat::SCTransform(cortex_sc, verbose = FALSE) %>%
  #   Seurat::RunPCA(., verbose = FALSE) %>%
  #   Seurat::RunUMAP(., dims = 1:30, verbose = FALSE)

  # Seurat::DimPlot(cortex_sc,
  #                 group.by = "subclass",
  #                 label = TRUE) + Seurat::NoLegend()

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
  # cellCount = rawSpatial[(length(rawSpatial)-15):length(rawSpatial)]
  cellCount = rawSpatial[(length(rawSpatial)-1):length(rawSpatial)]
  cellTypes = colnames(cellCount)

  #most bizarre cells in a cell type
  #find the index of the row
  max1 = -10
  index1 = 1
  index2 = 1
  max2 = -10
  for(row in 1: nrow(decon_mtrx)){
    if(cellCount[row, 1] == 1){
      if(abs(decon_mtrx[row, cellTypes[1]] - .5) > max1){
        max1 = abs(decon_mtrx[row, cellTypes[1]] - .5)
        index1 = row
      }
      if(abs(decon_mtrx[row, cellTypes[2]] - .5) > max1){
        max2 = abs(decon_mtrx[row, cellTypes[2]] - .5)
        index2 = row
      }
    }
  }


  #create a folder named singleCellNames[i] in "images/spotlight/twoTypes"
  # mkdir(paste0("images/spotlight/twoTypes/", singleCellNames[i]))
  dir.create(paste0("images/spotlight/twoTypes/", singleCellNames[i]))

  #save the directory name in paths variable
  paths <- paste0("images/spotlight/twoTypes/", singleCellNames[i], "/")

  #loop through cellTypes
  for(cellType in cellTypes){
    ttemp = paste(paths, cellType, ".png", sep = "")
    png(ttemp)
    scatterDf <- data.frame(x= cellCount[,cellType], y = decon_mtrx[,cellType], check.names = FALSE)
    sp = scatterPlotSpotlight(scatterDf)
    #save this plot to  "images/spotlight/twoTypes/singleCellNames[i]/cellType.png"
    # savePlot(sp, file = paste0("images/spotlight/twoTypes/", singleCellNames[i], "/", cellType, ".png"))
    dev.off()
  }


}
