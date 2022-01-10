dataset = read.csv(file="./data/merfish/s7.csv")
data = dataset[, colSums(dataset[,-c(1:5, 8,9)]) != 0]

counts = data[,-c(1:5, 8,9)]

coords = -1 * counts[,c(1,2)]
colnames(coords) = c("x", "y")
rownames(coords) = apply(coords, 1, function(i) {
   paste0(i[1], "x", i[2])
})

rownames(counts) = rownames(coords)
counts = t(counts[,-c(1,2)])

genes = colnames(counts)
message("Number of genes: ", length(genes))
clusters = unique(data[,8])
message("Number of clusters: ", length(clusters))

# sce = SingleCellExperiment(assays = list(counts = counts), colData = coords)
# sce = logNormCounts(sce)
# counts = logcounts(sce)

so <- CreateSeuratObject(counts = counts, project = "merfish")
                         # min.cells = 3, min.features = 200)

# clusterList = c()
# for (x in 1:ncol(so)) {
#   randVal = sample(1:6, 1)
#   clusterList = c(clusterList, randVal)
# }
pbmc = so


pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-") # for nothing i donno
#normalize data
pbmc <- NormalizeData(pbmc)
#variable gene finding
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#scaling
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
#Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
#clustering(it is not needed actually, will be rewritten according to our data)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 1.2)
#UAMP
pbmc <- RunUMAP(pbmc, dims = 1:10)


#Assign each cluster a number

clusterList = c()

for (x in 1:ncol(so)) {
  if(dataset[x,'Cell_class'] == "Astrocyte") clusterList = c(clusterList, 0)
  else if(dataset[x,'Cell_class'] == "Inhibitory") clusterList = c(clusterList, 1)
  else if(dataset[x,'Cell_class'] == "Pericytes") clusterList = c(clusterList, 2)
  else if(dataset[x,'Cell_class'] == "Ambiguous") clusterList = c(clusterList, 3)
  else if(dataset[x,'Cell_class'] == "Endothelial 1") clusterList = c(clusterList, 4)
  else if(dataset[x,'Cell_class'] == "Excitatory") clusterList = c(clusterList, 5)
  else if(dataset[x,'Cell_class'] == "OD Immature 1") clusterList = c(clusterList, 6)
  else if(dataset[x,'Cell_class'] == "OD Immature 2") clusterList = c(clusterList, 7)
  else if(dataset[x,'Cell_class'] == "Microglia") clusterList = c(clusterList, 8)
  else if(dataset[x,'Cell_class'] == "OD Mature 2") clusterList = c(clusterList, 9)
  else if(dataset[x,'Cell_class'] == "OD Mature 1") clusterList = c(clusterList, 10)
  else if(dataset[x,'Cell_class'] == "Endothelial 3") clusterList = c(clusterList, 11)
  else if(dataset[x,'Cell_class'] == "OD Mature 3") clusterList = c(clusterList, 12)
  else if(dataset[x,'Cell_class'] == "OD Mature 4") clusterList = c(clusterList, 13)
  else if(dataset[x,'Cell_class'] == "Endothelial 2") clusterList = c(clusterList, 14)
  else if(dataset[x,'Cell_class'] == "Ependymal") clusterList = c(clusterList, 15)
}


#overwrite part according to our merfish dataset

pbmc[['seurat_clusters']] <- clusterList
pbmc[['RNA_snn_res.1.2']] <- clusterList

#finding marker genes
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers


