library(BiocParallel)
library(doParallel)
source('functions.R')

# evaluate multiple methods for different batch numbers
num.batch <- 2:10
num.cell <- 200
facLoc <- 0.06
batchnum.path <- '/home/zy/single_cell/BEE/sim_data/batchnum'
batchnum.path <- '/lustre/tianlab/zhangyu/BEE/sim_data/batchnum'
batchnum.path <- '/home/drizzle_zhang/Desktop/single_cell/BEE/sim_data/batchnum'

# parallel computing
BPParam <- BatchtoolsParam(workers = 9)
batchnum.list <- do.call(
    bpmapply,
    c(list(FUN = generate_data_evaluate, SIMPLIFY = FALSE, USE.NAMES = FALSE,
           BPPARAM = BPParam,
           MoreArgs = list(num.cell = num.cell, facLoc = facLoc)),
      list(num.batch = num.batch)))

# registerDoParallel(cores = 10)
# batchnum.list <-
#     foreach(num.batch = num.batch, .combine = rbind) %dopar%
#     generate_data_evaluate(num.batch = num.batch)


# PCA and umap plot
for (i in 1:length(batchnum.list)) {
  if (i == 1) {
    df.batchnum <- batchnum.list[[i]][['access.res']]
  } else{
    df.batchnum <- rbind(df.batchnum, batchnum.list[[i]][['access.res']])
  }
  # title.pca <- paste0('BatchNum parameter: ', vector.facLoc[i])
  ggplot.pca <- DimPlot(
    batchnum.list[[i]][['object.pca']],
    reductions = 'pca',
    group.by = "batch",
    pt.size = 1
  )
  ggsave(
    plot = ggplot.pca, path = batchnum.path,
    filename = paste0('PCA_plot_BatchNum_', num.batch[i], '.png'),
    units = 'cm', width = 15, height = 8
  )
  ggplot.umap <- DimPlot(
    batchnum.list[[i]][['object.pca']], 
    reductions = 'umap', group.by = "batch", pt.size = 1
  )
  ggsave(
    plot = ggplot.pca, path = batchnum.path,
    filename = paste0('UMAP_plot_BatchNum_', num.batch[i], '.png'),
    units = 'cm', width = 15, height = 8
  )
  
}

# save results
path.write <- paste0(batchnum.path, 'batchnum_out.txt')
write.table(df.batchnum, file = path.write, quote = F, sep = '\t')

# plot
batchnum.names <- names(df.batchnum)
for (j in 1:dim(df.batchnum)[2]) {
    if (j == 1) {
        df.batchnum.plot <- data.frame(
            BatchNumParameter = num.batch, BatchEffectIndex = df.batchnum[,j],
            Label = rep(batchnum.names[j], length(num.batch)))
    } else {
        df.batchnum.plot <- rbind(
            df.batchnum.plot, 
            data.frame(
                BatchNumParameter = num.batch, 
                BatchEffectIndex = df.batchnum[,j],
                Label = rep(batchnum.names[j], length(num.batch))))
    }
}
batchnum.plot <- ggplot(
    df.batchnum.plot, 
    aes(x = BatchNumParameter, y = BatchEffectIndex, 
        color = Label, shape = Label)) + 
    geom_line() + geom_point(size = 2) + 
    labs(title = "evaluate multiple methods for different batch numbers")
ggsave(
    plot = batchnum.plot, path = batchnum.path, filename = "batchnum_out.png",
    units = 'cm', width = 25, height = 15)

