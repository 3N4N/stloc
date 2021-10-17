source("./functions.R")

library(SingleCellExperiment)
library(scater)
library(scran)

require(ggforce)
require(patchwork)
require(ggpubr)

if (!file.exists("output")) {
  system("mkdir -p output")
}

## Process Cancer dataset

counts_raw <- read.delim("./data/skin_cancer/GSE144239_ST_P2_S1_counts.tsv", header = TRUE, row.names = 1)

coords_raw <- do.call(rbind, strsplit(rownames(counts_raw), "x"))
coords <- apply(coords_raw, 1:2, as.numeric)
colnames(coords) <- c("x","y")
rownames(coords) <- rownames(counts_raw)

counts <- t(counts_raw)

sce <- SingleCellExperiment(assays = list(counts = counts), colData = coords)
sce <- logNormCounts(sce)
counts <- logcounts(sce)

### Determine highly variable genes

# sce = SingleCellExperiment(assays = list(counts = counts), colData = coords)
# sce <- logNormCounts(sce)
# dec <- modelGeneVar(sce)

# hvg <- getTopHVGs(dec,fdr.threshold = 0.05)
# hvg <- sort(hvg)
# length(hvg)

# seqvals = seq(min(dec$mean), max(dec$mean), length.out = 1000)
# peakExp = seqvals[which.max(metadata(dec)$trend(seqvals))]

# pdf(file = "./output/HVG_selection.pdf", height = 8, width = 8)
# plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
# curve(metadata(dec)$trend(x), col = "blue", add = TRUE)
# points(dec$mean[ which(rownames(dec) %in% hvg)],
#        dec$total[which(rownames(dec) %in% hvg)],
#        col = "red", pch = 16)
# abline(v = peakExp, lty = 2, col = "black")
# dev.off()


### Select for downstream analysis those marker genes which are also highly variable

clusterData <- read.delim("./data/skin_cancer/reference_markers_for_NMF.tsv", header = TRUE)
clusterGenes <- clusterData[,8]
clusterGenes <- unique(clusterGenes)
# commonGenes <- intersect(hvg, clusterGenes)
commonGenes <- intersect(rownames(counts), clusterGenes)
clusterData <- as.data.frame(clusterData)
write.table(clusterData[clusterData$gene %in% commonGenes,],
            file = "./data/common.tsv", row.names = FALSE, sep = "\t")

clusterData <- read.delim("./data/common.tsv", header = TRUE)
clusterNames <- unique(clusterData[,7])
clusterNames <- sapply(clusterNames, function(i) i <- toString(i))
clusterGenes <- clusterData[,8]
clusterGenePair <- list()
if (length(clusterNames) == 1) {
    clusterGenePair[[clusterNames[1]]] <- clusterData[clusterData[,7] == clusterNames[1], 8]
} else {
    clusterGenePair <- sapply(clusterNames, function(i) {
        clusterGenePair[[i]] <- clusterData[clusterData[,7] == i, 8]
  })
}


clusterGenePair[["MyeloidFibroblast"]] = unlist((list(clusterGenePair[["Myeloid"]], clusterGenePair[["Fibroblast"]])) )
clusterGenePair[["EpithelialFibroblast"]] = unlist((list(clusterGenePair[["Epithelial"]], clusterGenePair[["Fibroblast"]])) )
clusterGenePair[["MyeloidEpithelial"]] = unlist((list(clusterGenePair[["Myeloid"]], clusterGenePair[["Epithelial"]])) )



if (!file.exists("output/pval_plots_cancer")) {
  system("mkdir output/pval_plots_cancer")
}

if (!file.exists("output/pval_distribution")) {
  system("mkdir output/pval_distribution")
}

if (!file.exists("output/dump")) {
  system("mkdir output/dump")
}


plotdf = function(df_res, vals1, vals2, vals3, label1, label2, label3) {
# plotdf = function(df_res, vals1, vals2, label1, label2) {

    # for (i in 1:length(vals1)) {
    #   if (vals2[i] > 0.5) {
    #     vals1[i] <- NA
    #   }
    # }

    plot_vals1 <- ggplot(df_res, aes(x = x, y = -y)) +
        geom_point(aes(colour = vals1), size = 5) +
        theme_minimal() +
        theme(panel.grid = element_blank()) +
        theme(axis.text = element_blank()) +
        xlab("") +
        ylab("") +
        labs(colour = "") +
        theme(legend.position = "bottom") +
        theme(plot.title = element_text(hjust = 0.5, face = "italic")) +
        scale_color_viridis_c(option = "plasma") +
        coord_fixed() +
        guides(colour = guide_colourbar(title.position = "top",
                                        title.hjust = 0.5)) +
        theme(legend.key.width = unit(0.5, "inches")) +
        theme(plot.title = element_text(size = 20)) +
        theme(axis.title = element_text(size = 15)) +
        theme(legend.title = element_text(size = 15)) +
        labs(colour = label1) +
        NULL

    plot_vals2 <- ggplot(df_res, aes(x = x, y = -y)) +
        # geom_point(aes(colour = vals2), size = 5) +
        geom_point(aes(colour = -log10(vals2)), size = 5) +
        theme_minimal() +
        theme(panel.grid = element_blank()) +
        theme(axis.text = element_blank()) +
        xlab("") +
        ylab("") +
        labs(colour = "") +
        theme(legend.position = "bottom") +
        theme(plot.title = element_text(hjust = 0.5, face = "italic")) +
        scale_color_viridis_c(option = "plasma") +
        coord_fixed() +
        guides(colour = guide_colourbar(title.position = "top",
                                        title.hjust = 0.5)) +
        theme(legend.key.width = unit(0.5, "inches")) +
        theme(plot.title = element_text(size = 20)) +
        theme(axis.title = element_text(size = 15)) +
        theme(legend.title = element_text(size = 15)) +
        labs(colour = label2) +
        NULL

    plot_vals3 <- ggplot(df_res, aes(x = x, y = -y)) +
        geom_point(aes(colour = -log10(vals3)), size = 5) +
        theme_minimal() +
        theme(panel.grid = element_blank()) +
        theme(axis.text = element_blank()) +
        xlab("") +
        ylab("") +
        labs(colour = "") +
        theme(legend.position = "bottom") +
        theme(plot.title = element_text(hjust = 0.5, face = "italic")) +
        scale_color_viridis_c(option = "plasma") +
        coord_fixed() +
        guides(colour = guide_colourbar(title.position = "top",
                                        title.hjust = 0.5)) +
        theme(legend.key.width = unit(0.5, "inches")) +
        theme(plot.title = element_text(size = 20)) +
        theme(axis.title = element_text(size = 15)) +
        theme(legend.title = element_text(size = 15)) +
        labs(colour = label3) +
        NULL


    vals1_leg = as_ggplot(get_legend(plot_vals1))
    vals2_leg = as_ggplot(get_legend(plot_vals2))
    vals3_leg = as_ggplot(get_legend(plot_vals3))

    gridExtra::grid.arrange(plot_vals1 + theme(legend.position = "none")
                            + theme(plot.margin = margin(10,0,-10,0)),
                            plot_vals2 + theme(legend.position = "none")
                            + theme(plot.margin = margin(10,0,-10,0)),
                            plot_vals3 + theme(legend.position = "none")
                            + theme(plot.margin = margin(10,0,-10,0)),
                            vals1_leg, vals2_leg, vals3_leg,
                            layout_matrix = matrix(c(1,1,2,2,3,3,
                                                     1,1,2,2,3,3,
                                                     4,4,5,5,6,6),
                                                   ncol = 6,
                                                   byrow = TRUE))
}

ploteig = function(df.eig, vals, loc, label) {

    df.eig$vals[loc] = NA

    plot.vals = ggplot(df.eig, aes(x = x, y = -y)) +
        geom_point(aes(colour = vals), size = 4) +
        theme_minimal() +
        theme(panel.grid = element_blank()) +
        theme(axis.text = element_blank()) +
        xlab("") +
        ylab("") +
        labs(colour = "") +
        theme(legend.position = "bottom") +
        theme(plot.title = element_text(hjust = 0.5, face = "italic")) +
        scale_color_viridis_c(option = "plasma", na.value="red") +
        coord_fixed() +
        guides(colour = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
        theme(legend.key.width = unit(0.5, "inches")) +
        theme(plot.title = element_text(size = 20)) +
        theme(axis.title = element_text(size = 15)) +
        theme(legend.title = element_text(size = 15)) +
        labs(colour = label) +
        NULL

    vals.leg = as_ggplot(get_legend(plot.vals))

    gridExtra::grid.arrange(plot.vals + theme(legend.position = "none")
                            + theme(plot.margin = margin(10,0,-10,0)),
                            vals.leg,
                            layout_matrix = rbind(c(1,1,1),
                                                  c(1,1,1),
                                                  c(2,2,2)))
}


# W <- weightMatrix_nD(coords, span = 0.3)
W <- weightMatrix_gaussian(coords, l = 0.5)

# clusterNames = append(clusterNames,"EpithelialFibroblast" ,length(clusterNames))
# clusterNames = append(clusterNames,"MyeloidFibroblast" ,length(clusterNames))
# clusterNames = append(clusterNames,"MyeloidEpithelial" ,length(clusterNames))

set.seed(500)

for (nitr in c(1e3, 1e5)) {
    for (cluster in clusterNames) {
        # if (!(cluster=="Epithelial" | cluster=="Fibroblast" | cluster=="Myeloid")) next
        if (!(cluster=="Epithelial" | cluster=="Fibroblast")) next
        # if (!(cluster=="Epithelial")) next

        # if (!(cluster=="MyeloidFibroblast" | cluster=="EpithelialFibroblast" | cluster=="MyeloidEpithelial")) next
        # if (!(cluster=="MyeloidEpithelial" )) next

        genes <- unlist(c(clusterGenePair[cluster]))
        genes <- sapply(genes, function(i) i <- toString(i))
        if (length(genes) == 1) next

        print(cluster)
        print(genes)

        pairCount <- as.matrix(rbind(counts[genes,]))
        rownames(pairCount) <- genes
        # print(dim(pairCount))

        # st = Sys.time()
        # zscr <- as.matrix(sapply(1:nrow(coords),
        #                          function(i) {
        #                              zscores <- apply(pairCount, 1, scale)
        #                              aggzscores <- sum(zscores[i,])/nrow(pairCount)
        # }))
        # et = Sys.time()
        # message("Runtime of Z-score: ", et-st)

        st = Sys.time()
        meig <- as.matrix(sapply(1:nrow(coords), function(i) maxEigenVal(pairCount, W[i,])))
        et = Sys.time()
        message("Runtime of eigenvalues: ", et-st)

        # st = Sys.time()
        # msvd <- as.matrix(sapply(1:nrow(coords), function(i) maxSingVal(pairCount, W[i,])))
        # et = Sys.time()
        # message("Runtime of singular values: ", et-st)

        message(paste0("Conducting permutation tests for ", cluster))

        # pmeig <- matrix(nrow = nitr, ncol = nrow(coords))
        # pmeig <- sapply(1:nitr, function(i) {
        #   cat("\r", "Iteration step", i)
        #   x <- pairCount
        #   o = sample(1:nrow(coords))
        #   x <- t(sapply(1:nrow(pairCount), function(j) {
        #       x[j,] = pairCount[j,o]
        #   }))
        #   pmeig[i,] = sapply(1:nrow(W), function(j) maxEigenVal(x, W[j, ]))
        # })

        cnt = 1
        pmeig = c(1:nitr)
        for (i in 1:nitr) {
            cat("\r", "Iteration step", i)
            o = sample(1:nrow(coords))
            x = pairCount
            x = t(sapply(1:nrow(pairCount), function(j) {
                             x[j,] = pairCount[j,o]
                                                  }))
            c = coords
            c = sapply(1:ncol(coords), function(j) {
                           c[,j] = coords[o,j]
                            })
            randloc = sample(1:nrow(W), 1)
            pmeig[i] = maxEigenVal(x, W[randloc, ])
            cutoff = if (cluster == "Epithelial") 200 else 150
            if (pmeig[i] > cutoff & cnt <= 10) {
                tmeig = as.matrix(sapply(1:ncol(x), function(i) maxEigenVal(x, W[i,])))
                df.eig = data.frame(x = c[,1],
                                    y = c[,2],
                                    vals = tmeig)
                if (!file.exists(paste0("output/dump/", cluster, "x", log(nitr,10)))) {
                    system((paste0("mkdir output/dump/", cluster, "x", log(nitr,10))))
                }
                pdf(paste0("output/dump/", cluster, "x", log(nitr, 10), "/", cluster, "_", i, ".pdf"),
                    height = 6, width = 10, onefile = F)
                ploteig(df.eig, df.eig$meig, randloc, "Largest Eigenvalue")
                dev.off()
                cnt = cnt + 1
            }
        }
        cat("\n")
        message(paste0("Permutation tests for ", cluster, " completed"))

        # pvals_cor <- matrix(nrow = nrow(coords), ncol = 1)
        # pvals_cor <- as.matrix(sapply(1:nrow(wcor), function(i) {
        #   pvals_cor[i,] = (sum(pwcor[i,] > wcor[i]) + 1) / (nitr + 1)
        # }))
        pvals_eig <- matrix(nrow = nrow(coords), ncol = 1)
        pvals_eig <- as.matrix(sapply(1:nrow(meig), function(i) {
                                          # pvals_eig[i,] = (sum(pmeig[i,] > meig[i]) + 1) / (nitr + 1)
                                          pvals_eig[i,] = (sum(pmeig > meig[i]) + 1) / (nitr + 1)
                                                  }))

        fdrbh_eig = p.adjust(pvals_eig, method="BH")


        df_res <- data.frame(x = coords[,"x"],
                             y = coords[,"y"],
                             meig = meig,
                             pval = pvals_eig,
                             fdr = fdrbh_eig
        )

        pdf(paste0("output/pval_plots_cancer/", cluster, ".pdf"),
            height = 6, width = 10, onefile = F)
        # plotdf(df_res,df_res$zscr,df_res$meig,"Z-score", "Largest Eigenvalue")
        # plotdf(df_res,df_res$zscr,df_res$meig,"Z-score", "Max Singular Value")
        # plotdf(df_res,df_res$wcor,df_res$pvals_cor,"Weighted Correlation", "-log10(pval)")
        # plotdf(df_res,df_res$meig,df_res$pval,"Largest Eigenvalue","-log10(pval)")
        plotdf(df_res,df_res$meig,df_res$pval,df_res$fdr,"Largest Eigenvalue","-log10(pval)", "-log10(fdr)")
        dev.off()

        # break
    }
}
