library(dplyr)
library(Seurat)

#reading merfish data
dataset = read.csv(file="./data/merfish/s7.csv")
data = dataset[, colSums(dataset[,-c(1:9)]) != 0]
counts = data[,-c(1:9)]

#reading pbmc data
pbmc.data <- Read10X(data.dir = "./data/pbmc3k/filtered_gene_bc_matrices/hg19")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 2)
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

markerGenes = list()

for(x in 1:16){
  print(x)
  markerGenes[[x]] = intersect((pbmc.markers[pbmc.markers[6] == x-1,][[7]]), toupper(colnames(counts)))
}



