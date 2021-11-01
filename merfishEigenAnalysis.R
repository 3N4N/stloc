#  ----------------------------------------------------------------------
#                                  Load packages
#  ----------------------------------------------------------------------

library(SingleCellExperiment)
library(scran)

library(ggplot2)
library(gridExtra)


#  ----------------------------------------------------------------------
#                              Load helper functions
#  ----------------------------------------------------------------------

source("./functions.R")


#  ----------------------------------------------------------------------
#                            Create output directories
#  ----------------------------------------------------------------------

if (!file.exists("output")) {
    system("mkdir -p output")
}

if (!file.exists("output/merfish/plots")) {
    system("mkdir -p output/merfish/plots")
}
if (!file.exists("output/merfish/scatterPlots")) {
    system("mkdir -p output/merfish/scatterPlots")
}



counts.raw <- read.table("./data/merfish/merfishSpatial.csv", header = TRUE, row.names = 1)
coords.raw <- do.call(rbind, strsplit(rownames(counts.raw), "x"))
coords <- apply(coords.raw, 1:2, as.numeric)
colnames(coords) <- c("x","y")
rownames(coords) <- rownames(counts.raw)

counts <- t(counts.raw)

sce <- SingleCellExperiment(assays = list(counts = counts), colData = coords)
sce <- logNormCounts(sce)
counts <- logcounts(sce)


clusters.data <- read.table("./data/merfish/markerGene_for_merfish_data.csv", sep=",",header = TRUE)
clusters.data <- as.data.frame(clusters.data)
clusters.name <- unique(clusters.data$cell_type)
clusters.name <- sapply(clusters.name, function(i) i <- toString(i))

clusters.pair <- list()
if (length(clusters.name) == 1) {
    clusters.pair[[clusters.name[1]]] <- clusters.data[clusters.data[,1] == clusters.name[1], 2]
} else {
    clusters.pair <- sapply(clusters.name, function(i) {
        clusters.pair[[i]] <- clusters.data[clusters.data[,1] == i, 2]
  })
}

d <- sort (as.numeric (dist (coords )))[1]
W <- weightMatrix.gaussian(coords, l = 0.1)


set.seed(500)
# for (nitr in c(1e3)) {
for (nitr in c(1e3)) {
    brk = 0

    for (cluster in clusters.name) {
        

        # if (!(cluster=="Epithelial" | cluster=="Fibroblast" | cluster=="Myeloid")) next
        if (!(cluster=="Excitatory" | cluster=="Inhibitory" | cluster== "Astrocyte" | cluster==  "Inhibitory" | cluster== "Pericytes" | cluster== "Ambiguous" | cluster=="Endothelial 1"| cluster==  "Excitatory"| cluster=="OD Immature 1" | cluster=="OD Immature 2" | cluster== "Microglia" | cluster=="OD Mature 2" | cluster== "OD Mature 1" | cluster== "Endothelial 3" | cluster=="OD Mature 3" | cluster== "OD Mature 4" | cluster== "Endothelial 2" | cluster== "Ependymal")) next
        # if (!(cluster=="Astrocyte")) next
        genes <- unlist(c(clusters.pair[cluster]))
        genes <- sapply(genes, function(i) i <- toString(i))
        if (length(genes) == 1) next

        print(cluster)
        print(genes)

        pairCount <- as.matrix(rbind(counts[genes,]))
        rownames(pairCount) <- genes

        # st = Sys.time()
        # zscr <- as.matrix(sapply(1:nrow(coords),
        #                          function(i) {
        #                              zscores <- apply(pairCount, 1, scale)
        #                              aggzscores <- sum(zscores[i,])/nrow(pairCount)
        # }))
        # et = Sys.time()
        # message("Runtime of Z-score: ", difftime(et,st,units="mins"), " mins")

        st = Sys.time()
        meig.real <- as.matrix(sapply(1:nrow(coords),
                                      function(i) maxEigenVal(pairCount, W[i,])))
        et = Sys.time()
        # print(summary(meig.real))
        message("Runtime of eigenvalues: ", difftime(et,st,units="mins"), " mins")

        message(paste0("Conducting permutation tests for ", cluster))

        cnt = 1
        meig.perm = c(1:nitr)
        brk = 0
        st = Sys.time()
        for (i in 1:nitr) {
            if (brk) break
            cat("\r", "Iteration step", i)
            o = sample(1:nrow(coords))
            x = pairCount
            x = t(sapply(1:nrow(pairCount), function(j) {
                    x[j,] = pairCount[j,o]
            }))
            randloc = sample(1:nrow(coords), 1)
            meig.perm[i] = maxEigenVal(x, W[randloc,])
        }
        cat("\n")
        message(paste0("Permutation tests for ", cluster, " completed"))
        et = Sys.time()
        message("Runtime of ", cluster, " for ", nitr, " iterations: ", difftime(et,st,units="mins"), " mins")

        # if (brk == 1) next

        meig.pval <- matrix(nrow = nrow(coords), ncol = 1)
        meig.pval <- as.matrix(sapply(1:nrow(meig.real), function(i) {
                meig.pval[i,] = (sum(meig.perm > meig.real[i]) + 1) / (nitr + 1)
        }))

        meig.fdr = p.adjust(meig.pval, method="BH")

        df <- data.frame(x = coords[,"x"], y = coords[,"y"])
        
        pdf(paste0("output/merfish/plots/", cluster, ".pdf"),
            height = 6, width = 10, onefile = F)
        plotvals(3, df, vals=list(meig.real, -log10(meig.pval), -log10(meig.fdr)),
                 c("Largest Eigenvalue","-log10(pval)", "-log10(fdr)"), 3)
        # plotvals(1, df, vals=list(meig.real), c("Largest Eigenvalue"), 3)
        dev.off()

         scatterDf <- data.frame(x= counts.raw$cellCount, y = meig.real)


        pdf(paste0("output/merfish/scatterPlots/", cluster, ".pdf"),
            height = 8, width = 8, onefile = F)

        scatterPlot(scatterDf)
        dev.off()


    }
}




