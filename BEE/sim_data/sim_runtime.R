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
runtime.path <- '/lustre/tianlab/zhangyu/BEE/sim_data/runtime'
#runtime.path <- '/home/zy/single_cell/BEE/sim_data/runtime'

# parallel computing
BPParam <- BatchtoolsParam(workers = 10)
runtime.list <- do.call(
    bpmapply,
    c(list(FUN = calculate_evaluate_runtime, 
           SIMPLIFY = FALSE, USE.NAMES = FALSE, BPPARAM = BPParam,
           MoreArgs = list()),
      list(num.cell = vec.num.cell)))

# transfer to data.frame
for (i in 1:length(runtime.list)) {
    if (i == 1) {
        df.runtime <- runtime.list[[i]]
    } else{
        df.runtime <- rbind(df.runtime, runtime.list[[i]])
    }
}

# save results
path.write <- paste(runtime.path, 'runtime_out.txt', sep = '/')
write.table(df.runtime, file = path.write, quote = F, sep = '\t')

df.runtime <- read.table(path.write, sep = '\t')

# plot
vec.num.cell <- vec.num.cell * 2
runtime.names <- names(df.runtime)
for (j in 1:dim(df.runtime)[2]) {
    if (j == 1) {
        df.runtime.plot <- data.frame(
            CellNum = vec.num.cell, RunTime = df.runtime[,j],
            Label = rep(runtime.names[j], length(vec.num.cell)))
    } else {
        df.runtime.plot <- rbind(
            df.runtime.plot, 
            data.frame(
                CellNum = vec.num.cell, 
                RunTime = df.runtime[,j],
                Label = rep(runtime.names[j], length(vec.num.cell))))
    }
}
runtime.plot <- 
    ggplot(
        df.runtime.plot, 
        aes(x = CellNum, y = RunTime, 
            color = Label, shape = Label)) + 
    geom_line() + 
    geom_point(size = 2) + 
    labs(title = "Evaluate Runtime of multiple methods for different cell numbers") +
    scale_colour_discrete(name = 'Evaluation methods', 
                          breaks = c('pcReg', 'sil', 'kBET', 'mixent', 
                                     'ARI', 'NMI', 'ldaReg'), 
                          labels = c('PC Regression', 'Silhouettes', 'kBET',
                                     'Entropy of batch mixing', 
                                     'Adjusted rand index', 
                                     'Mormalized mutual information',
                                     'LDA Regression')) + 
    scale_shape_manual(values = c(0, 1, 2, 4, 20, 5, 7),
                       name = 'Evaluation methods', 
                       breaks = c('pcReg', 'sil', 'kBET', 'mixent', 
                                  'ARI', 'NMI', 'ldaReg'), 
                       labels = c('PC Regression', 'Silhouettes', 'kBET',
                                  'Entropy of batch mixing', 
                                  'Adjusted rand index', 
                                  'Mormalized mutual information',
                                  'LDA Regression'))
ggsave(
    plot = runtime.plot, path = runtime.path, filename = "runtime_out.png",
    units = 'cm', width = 30, height = 15)







