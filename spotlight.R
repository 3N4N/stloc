library(Matrix)
library(data.table)
library(Seurat)
library(SeuratData)
library(dplyr)
library(gt)
library(SPOTlight)
library(igraph)
library(RColorBrewer)
 

counts.raw <- read.table("./data/merfish/s7.csv", header = TRUE, sep = ',')
data = subset(counts.raw, select = -c(1, 2, 3, 4, 5, 6, 7, 8, 9))
cortex_sc <- CreateSeuratObject(counts = t(data))

cortex_sc[["subclass"]] = counts.raw[[8]]


for(i in 1:ncol(data)){
  print(i)
  data[is.infinite(data[,i]), i] <- mean(data[,i], na.rm = TRUE)
}

for(i in 1:ncol(data)) {  
  if(!is.finite(data[, i])){
    print(data[, i])
  }
}


data = do.call(data.frame,lapply(DT, function(x) replace(x, is.infinite(x),NA)))








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


start_time <- Sys.time()
nmf_mod_ls <- train_nmf(cluster_markers = cluster_markers_all, 
                        se_sc = se_sc_down, 
                        mtrx_spatial = anterior@assays$RNA@counts,
                        clust_vr = "subclass",
                        ntop = NULL,
                        hvg = 3000,
                        transf = "uv",
                        method = "nsNMF")

nmf_mod <- nmf_mod_ls[[1]]


spot_counts <- anterior@assays$RNA@counts


for (r in 1:row) {
  for(c in 1:col){
    decon_mtrx1[r, c] = (decon_mtrx[r, c] / sum(decon_mtrx[r,])) * 100
  }
}