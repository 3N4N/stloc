if !requireNamespace("BiocManager", quietly = TRUE)
    install.packages("BiocManager")


packages.cran = c("ggforce", "patchwork", "ggpubr", "Seurat")
packages.bioc = c("SingleCellExperiment", "scater", "scran")

install.packages(packages.cran, method="wget")
BiocManager::install(packages.bioc, method="wget")
