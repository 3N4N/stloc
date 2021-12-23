if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")


packages.cran = c("ggplot2", "gridExtra", "ks")
packages.bioc = c("SingleCellExperiment", "scran")

install.packages(packages.cran, method="wget")
BiocManager::install(packages.bioc, method="wget")
