if !requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleCellExperiment")
BiocManager::install("scater")
BiocManager::install("scran")

install.packages("ggforce")
install.packages("patchwork")
install.packages("ggpubr")

