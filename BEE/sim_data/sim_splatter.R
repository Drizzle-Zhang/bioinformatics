library(splatter)
library(scater)

params1 <- newSplatParams(batchCells = c(50, 50))
params2 <- newSplatParams(batchCells = c(50, 50), batch.facLoc = c(0.2, 0.2))
params3 <- newSplatParams(batchCells = c(50, 50, 50))
params4 <- newSplatParams(batchCells = c(50, 50), batch.facLoc = c(1, 1),
                          batch.facScale = c(1, 1))


sim1 <- splatSimulate(params1, verbose = F)
sim1 <- normalize(sim1)
plotPCA(sim1, colour_by = "Batch")
sim2 <- splatSimulate(params2, verbose = F)
sim2 <- normalize(sim2)
plotPCA(sim2, colour_by = "Batch")
sim3 <- splatSimulate(params3, verbose = F)
sim3 <- normalize(sim3)
plotPCA(sim3, colour_by = "Batch")
sim4 <- splatSimulate(params4, verbose = F)
sim4 <- normalize(sim4)
plotPCA(sim4, colour_by = "Batch")
