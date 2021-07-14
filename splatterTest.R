library("splatter")
library("scater")
library("ggplot2")
set.seed(1)
sce <- mockSCE()
params <- newSplatParams()
params
counts <- counts(sce)
counts
sim <- splatSimulate(params, batchCells = c(100, 100), verbose = FALSE)

sim.groups <- splatSimulate(group.prob = c(0.3, 0.5,0.2), method = "groups",
                            verbose = FALSE)
sim.groups <- normalize(sim.groups)
plotPCA(sim.groups, colour_by = "Group")
