library(BiocParallel)
library(doParallel)
source('functions.R')

# evaluate multiple methods for imbalance data
num.batch <- 2
num.cell <- 200
facLoc <- 0.06
ratios <- c(1, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05, 0.02)
list.cell <- list()
for (i in 1:length(ratios)) {
    list.cell[[i]] <- c(num.cell, num.cell*ratios[i])
}
var.path <- '/home/zy/single_cell/BEE/sim_data/variance'
var.path <- '/lustre/tianlab/zhangyu/BEE/sim_data/variance'
# var.path <- '/home/drizzle_zhang/Desktop/single_cell/BEE/sim_data/variance'

# parallel computing
BPParam <- BatchtoolsParam(workers = 8)
var.list <- do.call(
    bpmapply,
    c(list(FUN = generate_data_evaluate, SIMPLIFY = FALSE, USE.NAMES = FALSE,
           BPPARAM = BPParam,
           MoreArgs = list(num.batch = num.batch, facLoc = facLoc)),
      list(num.cell = list.cell)))





