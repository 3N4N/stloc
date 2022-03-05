#  ----------------------------------------------------------------------
#  This file contains the code for eigenvalue analysis. Other files can
#  preprocess the dataset and call the `analyze` function here.
#  ----------------------------------------------------------------------



source("./functions.R")


analyze <- function(clusters.name, clusters.pair, counts, cellCount, outdir) {
    counts.flat <- colSums(counts)

    data <- rep(names(counts.flat), counts.flat)
    data <- do.call(rbind, strsplit(data, "x"))
    data <- apply(data, c(1, 2), as.numeric)
    colnames(data) <- c("x", "y")


    # hpi <- Hscv(x = data)
    # kde <- kde(x = data, H = hpi, eval.points = coords)


    d <- sort(as.numeric(dist(coords)))[1]
    W <- weightMatrix.gaussian(coords, l = d * 0.1)

    set.seed(500)
    nitr <- 1e3

    # _________________________Cell type names for two datasets_______________________
    # if (!(cluster=="Excitatory" | cluster=="Inhibitory" | cluster== "Astrocyte" | cluster==  "Inhibitory" | cluster== "Pericytes" | cluster== "Ambiguous" | cluster=="Endothelial1"| cluster==  "Excitatory"| cluster=="ODImmature1" | cluster=="ODImmature2" | cluster== "Microglia" | cluster=="ODMature2" | cluster== "ODMature1" | cluster== "Endothelial3" | cluster=="ODMature3" | cluster== "ODMature4" | cluster== "Endothelial2" | cluster== "Ependymal")) next
    # if (!(cluster=="Epithelial" | cluster=="Fibroblast" | cluster=="Myeloid")) next
    mat <- matrix(ncol = 0, nrow = ncol(counts))

    resultFrame <- data.frame(mat)
    print(clusters.name)

    for (cluster in clusters.name) {

        # if (!(cluster == "Excitatory")) next
        genes <- unlist(c(clusters.pair[cluster]))
        print(genes)
        genes <- sapply(genes, function(i) i <- toString(i))
        # if (length(genes) == 1) next

        message("Analyzing cluster ", cluster)
        message("Number of marker genes: ", length(genes))

        pairCount <- as.matrix(rbind(counts[genes, ]))
        rownames(pairCount) <- genes

        # st <- Sys.time()
        # zscores <- apply(pairCount, 1, scale)
        # zscr <- as.matrix(sapply(1:nrow(coords),
        #                          function(i) {
        #                              aggzscores <- sum(zscores[i,])/nrow(pairCount)
        # }))
        # et <- Sys.time()
        # message("Runtime of Z-score: ", difftime(et,st,units="mins"), " mins")

        st <- Sys.time()
        meig.real <- as.matrix(sapply(
            1:nrow(coords),
            # function(i) (maxByMinSingVal(pairCount, W[i, ]))
            function(i) (maxSingVal(pairCount, W[i, ]))
        ))
        resultFrame[cluster] <- meig.real
        # function(i) (maxEigenVal(pairCount, W[i,]) * kde$estimate[i])))
        # tmp <- colSums(pairCount) #/ kde$estimate
        # for (i in 1:length(meig.real)) meig.real[i] <- sum(tmp * W[i,])
        et <- Sys.time()

        message("Runtime of eigenvalues: ", difftime(et, st, units = "mins"), " mins")

        # message(paste0("Conducting permutation tests for ", cluster))

        # meig.perm <- c(1:nitr)
        # st <- Sys.time()
        # for (i in 1:nitr) {
        #     cat("\r", "Iteration step", i)
        #     o <- sample(1:nrow(coords))
        #     x <- pairCount[, o]
        #     randloc <- sample(1:nrow(coords), 1)
        #     meig.perm[i] <- maxSingVal(x, W[randloc, ])
        #     # meig.perm[i] <- maxEigenVal(x, W[randloc,]) * kde$estimate[o[randloc]]
        #     # meig.perm[i] <- colSums(x)[randloc] / kde$estimate[o[randloc]]
        #     # meig.perm[i] <- sum(colSums(x) #/ kde$estimate[o[randloc]] * W[randloc,])
        # }
        # cat("\n")
        # message(paste0("Permutation tests for ", cluster, " completed"))
        # et <- Sys.time()
        # message("Runtime of ", cluster, " for ", nitr, " iterations: ", difftime(et, st, units = "mins"), " mins")

        # meig.pval <- meig.real
        # for (i in 1:nrow(coords)) {
        #     meig.pval[i] <- (sum(meig.perm > meig.real[i]) + 1) / (nitr + 1)
        # }

        # meig.fdr <- p.adjust(meig.pval, method = "BH")

        # df <- data.frame(x = coords[, "x"], y = -coords[, "y"])

        # pdf(paste0(outdir, cluster, "x", log(nitr, 10), ".pdf"),
        #     height = 6, width = 10, onefile = F
        # )
        # # plotvals(2, df, "", vals=list(meig.real, -log10(meig.pval)), c("Largest Eigenvalue",-log10(meig.pval)), 3, 1, 2)

        # plotvals(3, df, cluster,
        #     vals = list(meig.real, -log10(meig.pval), -log10(meig.fdr)),
        #     c("Largest Eigenvalue", "-log10(pval)", "-log10(fdr)"), 3, 1, 3
        # )
        # dev.off()
        pdf(paste0(outdir, cluster, ".pdf"),
            height = 6, width = 10
        )
        scatterDf <- data.frame(x = cellCount[, cluster], y = meig.real)
        scatterPlot(scatterDf)
        dev.off()
    }
    SVDcor <- resultFrame
    return(SVDcor)
}
