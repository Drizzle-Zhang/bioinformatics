# test two-dims method
library(BiocParallel)
library(ggplot2)
source('functions.R')
source('functions_evaluator.R')

# parameters
num.batch <- 2
num.cell <- 200
num.gene <- 15000
vec.facLoc <- c(0.03, 0.06, 0.1, 0.15)
facScale <- 0
num.pc <- 30
num.group <- 3
group.prob <- c(1/3, 1/3, 1/3)
de.prob <- 0.1
vec.de.facLoc <- c(0.1, 0.15, 0.25, 0.4)
de.facScale <- 0
seed.splatter <- 1234
path <- '/lustre/tianlab/zhangyu/BEE/sim_data/twodims'

# combine variables
input <- list()
idx <- 0
for (facLoc in vec.facLoc) {
    for (de.facLoc in vec.de.facLoc) {
        idx <- idx + 1
        input.list <- list(facLoc = facLoc, de.facLoc = de.facLoc)
        input[[idx]] <- input.list
    }
}

# main function
main <- function(
    input, facScale, num.batch, num.group, group.prob, de.prob, 
    de.facScale, num.cell, num.gene, seed.splatter, num.pc) {
    facLoc <- input$facLoc
    de.facLoc <- input$de.facLoc
    # simulation
    list.sim <- simulate.count.data(
        facLoc = facLoc, facScale = facScale, num.batch = num.batch, 
        num.group = num.group, group.prob = group.prob, 
        de.prob = de.prob, de.facLoc = de.facLoc, de.facScale = de.facScale,
        num.cell = num.cell, num.gene = num.gene, seed.splatter = seed.splatter)
    # prepare data
    object <- preprocess.data(
        mtx.count = list.sim$mtx.count, batches = list.sim$batches, 
        groups = list.sim$groups, num.pc = num.pc)
    # evaluation
    output <- evaluate.two.dims(object)
    # PCA and umap plot
    ggplot.pca <- DimPlot(
        object,
        reductions = 'pca',
        group.by = "group",
        shape.by = 'batch',
        pt.size = 1.5
    )
    ggsave(
        plot = ggplot.pca, path = path,
        filename = paste0('PCA_plot_', facLoc, '_', de.facLoc, '.png'),
        units = 'cm', width = 15, height = 8
    )
    ggplot.umap <- DimPlot(
        object, reductions = 'umap', 
        group.by = "group", shape.by = 'batch', pt.size = 1.5
    )
    ggsave(
        plot = ggplot.umap, path = path,
        filename = paste0('UMAP_plot_', facLoc, '_', de.facLoc, '.png'),
        units = 'cm', width = 15, height = 8
    )
    
    output$facLoc <- facLoc
    output$de.facLoc <- de.facLoc
    
    return(output)
}

BPParam <- BatchtoolsParam(workers = 16)
twodims.list <- do.call(
    bpmapply,
    c(list(FUN = main, SIMPLIFY = FALSE, USE.NAMES = FALSE,
           BPPARAM = BPParam,
           MoreArgs = list(facScale = facScale, 
                           num.batch = num.batch, 
                           num.group = num.group, 
                           group.prob = group.prob, 
                           de.prob = de.prob, 
                           de.facScale = de.facScale, num.cell = num.cell, 
                           num.gene = num.gene, seed.splatter = seed.splatter, 
                           num.pc = num.pc)),
      list(input = input)))

# results
df.twodims <- data.frame()
for (i in 1:length(twodims.list)) {
    df.twodims <- rbind(df.twodims, twodims.list[[i]])
}

# save results
path.write <- paste(path, 'diffexp_out.txt', sep = '/')
write.table(df.twodims, file = path.write, quote = F, sep = '\t')

df.twodims <- read.table(path.write, sep = '\t')

# plot
df.twodims$Label <- 
    c(paste0("(", df.twodims$facLoc, ",", df.twodims$de.facLoc, ")"))
twodims.plot <- 
    ggplot(
        df.twodims, 
        aes(x = batch.effect.factor, y = cell.distance, color = Label)
        ) + geom_point(size = 2) +
    labs(title = "Evaluate Batch effect in two dimensions",
         x = 'Degree of uniform mixing', y = 'Cell distance')
ggsave(
    plot = twodims.plot, path = path, filename = "twodims_result.png",
    units = 'cm', width = 30, height = 15)


