library(ggplot2)
library(SPOTlight)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)

source("R/data.R")

sce <- mockSC(ng = 200, nc = 10, nt = 3)
spe <- mockSP(sce)
mgs <- getMGS(sce)