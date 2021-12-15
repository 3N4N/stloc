library("splatter")
library("scater")
library("ggplot2")
set.seed(1)
sce <- mockSCE()
params <- newSplatParams()
# counts <- counts(sce)
# counts
sim <- splatSimulate(params, batchCells = c(100, 100,400,25), verbose = FALSE,nGenes=1000)

# sim.groups <- splatSimulate(group.prob = c(0.3, 0.5,0.2), method = "groups",
#                             verbose = FALSE)
# sim.groups <- normalize(sim.groups)
# sim.groups
sim <- logNormCounts(sim)
simPCA <- runPCA(sim)
mat<- counts(sim)
mat = t(mat)
p<-colData(sim)
df_mtx <- as.data.frame(mat)
df_mtx["cellType"]<-as.list(p['Batch'])
# mat<-  cbind(mat,Class= as.list(p['Batch'])
r<-rownames(mat)
write.csv(df_mtx, file = "splatterCSV.csv", append = FALSE, quote = TRUE, sep = ",")

plotPCA(simPCA, colour_by = "Batch")
