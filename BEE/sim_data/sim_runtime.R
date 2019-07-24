library(BiocParallel)
source('functions_test_runtime.R')

# calculate runtime of multiple methods
num.batch <- 2
# vec.num.cell <- c(100, 200, 300, 400, 500, 600, 800, 1000, 1200)
vec.num.cell <- c(100, 200, 400, 700, 1000, 1500, 2000, 3000, 4000, 5000)
# vec.num.cell <- 
#     c(100, 200, 400, 700, 1000, 1500, 2000, 3000, 5000, 10000, 20000)
num.gene <- 15000
facLoc <- 0.06
facScale <- 0
num.pc <- 30
num.cell = 1500
seed.splatter = 1234
var.path <- '/lustre/tianlab/zhangyu/BEE/sim_data/runtime'
#var.path <- '/home/zy/single_cell/BEE/sim_data/runtime'

# parallel computing
BPParam <- BatchtoolsParam(workers = 10)
var.list <- do.call(
    bpmapply,
    c(list(FUN = calculate_evaluate_runtime, 
           SIMPLIFY = FALSE, USE.NAMES = FALSE, BPPARAM = BPParam,
           MoreArgs = list()),
      list(num.cell = vec.num.cell)))





