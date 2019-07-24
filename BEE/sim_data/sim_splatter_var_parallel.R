library(BiocParallel)
library(doParallel)
source('functions.R')

# evaluate multiple methods for different variances
num.batch <- 2
num.cell <- 200
num.gene <- 15000
vector.facLoc <- seq(0, 0.2, 0.002)
#vector.facLoc <- seq(0.02, 0.08, 0.002)
facScale <- 0
num.pc <- 30
num.group <- 1
group.prob <- NULL
de.prob <- 0.1
de.facLoc <- 0.05
de.facScale <- 0
seed.splatter <- 1234
#var.path <- '/home/zy/single_cell/BEE/sim_data/variance'
var.path <- '/lustre/tianlab/zhangyu/BEE/sim_data/variance_2'
#var.path <- '/home/drizzle_zhang/Desktop/single_cell/BEE/sim_data/variance'

# main function
main <- function(
    facLoc, facScale, num.batch, num.group, group.prob, de.prob, de.facLoc, 
    de.facScale, num.cell, num.gene, seed.splatter, num.pc, method) {
    # simulation
    list.sim <- simulate.count.data(
        facLoc = facLoc, facScale = facScale, num.batch = num.batch, 
        num.group = num.group, group.prob = group.prob, 
        de.prob = de.prob, de.facLoc = de.facLoc, de.facScale = de.facScale,
        num.cell = num.cell, num.gene = num.gene, seed.splatter = seed.splatter)
    # prepare data
    object <- prepare_data(
        mtx.count = list.sim$mtx.count, batches = list.sim$batches, 
        groups = list.sim$groups, num.pc = num.pc)
    # evaluation
    output <- evaluate.batch.effect(object, method = method)
    
    return(output)
}

# parallel computing
BPParam <- BatchtoolsParam(workers = 34)
var.list <- do.call(
    bpmapply,
    c(list(FUN = main, SIMPLIFY = FALSE, USE.NAMES = FALSE,
           BPPARAM = BPParam,
           MoreArgs = list(facScale = facScale, num.batch = num.batch, 
                           num.group = num.group, group.prob = group.prob, 
                           de.prob = de.prob, de.facLoc = de.facLoc, 
                           de.facScale = de.facScale, num.cell = num.cell, 
                           num.gene = num.gene, seed.splatter = seed.splatter, 
                           num.pc = num.pc, method = 'batch')),
      list(facLoc = vector.facLoc)))

# registerDoParallel(cores = 10)
# var.list <-
#     foreach(facLoc = vector.facLoc, .combine = rbind) %dopar%
#     generate_data_evaluate(facLoc)

# calculate results
# var.list <- list()
# for (i in 1:length(vector.facLoc)) {
#     var.list[[i]] <-
#         generate_data_evaluate(vector.facLoc[i], vector.facScale, var.path)
#     if (i == 1) {
#         df.var <- var.list[[i]]
#     } else{
#         df.var <- rbind(df.var, var.list[[i]])
#     }
# }

# PCA and umap plot
for (i in 1:length(var.list)) {
  if (i == 1) {
    df.var <- var.list[[i]][['access.res']]
  } else{
    df.var <- rbind(df.var, var.list[[i]][['access.res']])
  }
  # title.pca <- paste0('Variance parameter: ', vector.facLoc[i])
  ggplot.pca <- DimPlot(
      var.list[[i]][['object.pca']],
      reductions = 'pca',
      group.by = "batch",
      pt.size = 1
  )
  ggsave(
      plot = ggplot.pca, path = var.path,
      filename = paste0('PCA_plot_Variance_', vector.facLoc[i], '.png'),
      units = 'cm', width = 15, height = 8
  )
  ggplot.umap <- DimPlot(
      var.list[[i]][['object.pca']], 
      reductions = 'umap', group.by = "batch", pt.size = 1
      )
  ggsave(
      plot = ggplot.pca, path = var.path,
      filename = paste0('UMAP_plot_Variance_', vector.facLoc[i], '.png'),
      units = 'cm', width = 15, height = 8
  )
  
}

# save results
path.write <- paste(var.path, 'var_out.txt', sep = '/')
write.table(df.var, file = path.write, quote = F, sep = '\t')

df.var <- read.table(path.write, sep = '\t')

# plot
var.names <- names(df.var)
for (j in 1:dim(df.var)[2]) {
    if (j == 1) {
        df.var.plot <- data.frame(
            VarianceParameter = vector.facLoc, BatchEffectIndex = df.var[,j],
            Label = rep(var.names[j], length(vector.facLoc)))
    } else {
        df.var.plot <- rbind(
            df.var.plot, 
            data.frame(
                VarianceParameter = vector.facLoc, 
                BatchEffectIndex = df.var[,j],
                Label = rep(var.names[j], length(vector.facLoc))))
    }
}
var.plot <- 
    ggplot(
    df.var.plot, 
    aes(x = VarianceParameter, y = BatchEffectIndex, 
        color = Label, shape = Label)) + 
    geom_line() + 
    geom_point(size = 2) + 
    labs(title = "Evaluate multiple methods for different variances") +
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
    plot = var.plot, path = var.path, filename = "var_out.png",
    units = 'cm', width = 30, height = 15)

# effective intervel and correlation coefficient
intervel <- function(df.method) {
    x <- df.method[, 1]
    y <- df.method[, 2]
    len <- length(x)
    diff <- c()
    for (i in 1:len) {
        if (i = 1) {
            diff <- c(diff, (y[i + 1] - y[i])/(x[i + 1] - x[i]))
        }
        if (i = len - 1) {
            diff <- c(diff, (y[i] - y[i - 1])/(x[i] - x[i - 1]))
        }
        diff <- c(diff, (y[i + 1] - y[i - 1])/(x[i + 1] - x[i - 1]))
    }
    
}
df.eval <- data.frame()
df.eval.plot <- data.frame()
for (method in unique(df.var.plot$Label)) {
    df.method <- df.var.plot[df.var.plot$Label == method,]
    start.interval <- 
        max(df.method$VarianceParameter[df.method$BatchEffectIndex < 0.05])
    end.interval <- 
        max(df.method$VarianceParameter[df.method$BatchEffectIndex < 0.95])
    var.interval <- 
        intersect(df.method$VarianceParameter[
            df.method$VarianceParameter > start.interval],
            df.method$VarianceParameter[
                df.method$VarianceParameter < end.interval])
    index.interval <- 
        df.method$BatchEffectIndex[
            df.method$VarianceParameter %in% var.interval]
    pcc = cor(var.interval, index.interval, method = 'pearson')
    scc = cor(var.interval, index.interval, method = 'spearman')
    kcc = cor(var.interval, index.interval, method = 'kendall')
    df.eval <- rbind(df.eval, 
                     data.frame(method = method, 
                                start = start.interval, end = end.interval,
                                pcc = pcc, scc = scc, kcc = kcc))
    df.eval.plot <- 
        rbind(
            df.eval.plot, 
            data.frame(method = method, point = start.interval, pcc = pcc))
    df.eval.plot <- 
        rbind(
            df.eval.plot, 
            data.frame(method = method, point = end.interval, pcc = pcc))
    }


eval.plot <- 
    ggplot(
        df.eval.plot, 
        aes(x = point, y = pcc, color = method, shape = method)) + 
    geom_line() + 
    geom_point(size = 2) + 
    labs(title = "Evaluate multiple methods for different variances",
         x = 'Effective intervel', y = 'PCC') +
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
    plot = eval.plot, path = var.path, filename = "eval_out.png",
    units = 'cm', width = 25, height = 15)

# interval
interval <- c()
for (i in 1:dim(df.var.plot)[1]) {
    if (df.var.plot$VarianceParameter[i] >= 0 && 
        df.var.plot$VarianceParameter[i] < 0.05) {
        interval <- c(interval, '0-0.05')
    }
    if (df.var.plot$VarianceParameter[i] >= 0.05 && 
        df.var.plot$VarianceParameter[i] < 0.1) {
        interval <- c(interval, '0.05-0.1')
    }
    if (df.var.plot$VarianceParameter[i] >= 0.1 && 
        df.var.plot$VarianceParameter[i] < 0.15) {
        interval <- c(interval, '0.1-0.15')
    }
    if (df.var.plot$VarianceParameter[i] >= 0.15 && 
        df.var.plot$VarianceParameter[i] <= 0.2) {
        interval <- c(interval, '0.15-0.2')
    }
    
}
df.var.plot$Interval <- interval
df.barplot <- data.frame()
for (method in unique(df.var.plot$Label)) {
    subdata_m <- df.var.plot[df.var.plot$Label == method,]
    for (range in unique(df.var.plot$Interval)) {
        subdata <- subdata_m[df.var.plot$Interval == range,]
        pcc <- cor(subdata$VarianceParameter, subdata$BatchEffectIndex, 
                   method = 'pearson', use = 'complete.obs')
        scc <- cor(subdata$VarianceParameter, subdata$BatchEffectIndex, 
                   method = 'spearman', use = 'complete.obs')
        kcc <- cor(subdata$VarianceParameter, subdata$BatchEffectIndex, 
                   method = 'kendall', use = 'complete.obs')
        df.barplot <- 
            rbind(df.barplot,
                  data.frame(Method = method, Interval = range,
                             PCC = pcc, SCC = scc, KCC = kcc))
    }
    pcc <- cor(subdata_m$VarianceParameter, subdata_m$BatchEffectIndex, 
               method = 'pearson', use = 'complete.obs')
    scc <- cor(subdata_m$VarianceParameter, subdata_m$BatchEffectIndex, 
               method = 'spearman', use = 'complete.obs')
    kcc <- cor(subdata_m$VarianceParameter, subdata_m$BatchEffectIndex, 
               method = 'kendall', use = 'complete.obs')
    df.barplot <- 
        rbind(df.barplot,
              data.frame(Method = method, Interval = '0-0.2',
                         PCC = pcc, SCC = scc, KCC = kcc))
}
df.barplot[is.na(df.barplot)] <- 0

bar.plot <- 
    ggplot(
        df.barplot, 
        aes(x = Interval, y = KCC, fill = Method)) + 
    geom_bar(position = 'dodge', stat = 'identity') + 
    labs(title = "Evaluate multiple methods for different variances",
         y = 'Kendall Correlation Coefficient') +
    scale_fill_discrete(name = 'Evaluation methods', 
                          breaks = c('pcReg', 'sil', 'kBET', 'mixent', 
                                     'ARI', 'NMI', 'ldaReg'), 
                          labels = c('PC Regression', 'Silhouettes', 'kBET',
                                     'Entropy of batch mixing', 
                                     'Adjusted rand index', 
                                     'Mormalized mutual information',
                                     'LDA Regression'))

ggsave(
    plot = bar.plot, path = var.path, filename = "barplot_KCC.png",
    units = 'cm', width = 25, height = 15)

bar.plot <- 
    ggplot(
        df.barplot, 
        aes(x = Interval, y = SCC, fill = Method)) + 
    geom_bar(position = 'dodge', stat = 'identity') + 
    labs(title = "Evaluate multiple methods for different variances",
         y = 'Spearman Correlation Coefficient') +
    scale_fill_discrete(name = 'Evaluation methods', 
                        breaks = c('pcReg', 'sil', 'kBET', 'mixent', 
                                   'ARI', 'NMI', 'ldaReg'), 
                        labels = c('PC Regression', 'Silhouettes', 'kBET',
                                   'Entropy of batch mixing', 
                                   'Adjusted rand index', 
                                   'Mormalized mutual information',
                                   'LDA Regression'))

ggsave(
    plot = bar.plot, path = var.path, filename = "barplot_SCC.png",
    units = 'cm', width = 25, height = 15)




