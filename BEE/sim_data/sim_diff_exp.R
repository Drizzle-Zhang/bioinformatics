setwd('/lustre/tianlab/zhangyu/my_git/bioinformatics/BEE/sim_data')
library(BiocParallel)
source('functions.R')

# evaluate multiple methods for different degree of different expression
num.batch <- 2
num.cell <- 200
num.gene <- 15000
facLoc <- 0.06
facScale <- 0
num.pc <- 30
num.group <- 3
group.prob <- c(1/3, 1/3, 1/3)
de.prob <- 0.1
de.facLoc <- seq(0, 0.5, 0.005)
de.facScale <- 0
seed.splatter <- 1234
method = 'group'
#diffexp.path <- '/home/zy/single_cell/BEE/sim_data/variance'
diffexp.path <- '/lustre/tianlab/zhangyu/BEE/sim_data/diffexp_3'
#diffexp.path <- '/home/drizzle_zhang/Desktop/single_cell/BEE/sim_data/variance'

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
diffexp.list <- do.call(
    bpmapply,
    c(list(FUN = main, SIMPLIFY = FALSE, USE.NAMES = FALSE,
           BPPARAM = BPParam,
           MoreArgs = list(facLoc = facLoc, facScale = facScale, 
                           num.batch = num.batch, 
                           num.group = num.group, 
                           group.prob = group.prob, 
                           de.prob = de.prob, 
                           de.facScale = de.facScale, num.cell = num.cell, 
                           num.gene = num.gene, seed.splatter = seed.splatter, 
                           num.pc = num.pc, method = method)),
      list(de.facLoc = de.facLoc)))

# registerDoParallel(cores = 10)
# diffexp.list <-
#     foreach(facLoc = vector.facLoc, .combine = rbind) %dopar%
#     generate_data_evaluate(facLoc)

# calculate results
# diffexp.list <- list()
# for (i in 1:length(vector.facLoc)) {
#     diffexp.list[[i]] <-
#         generate_data_evaluate(vector.facLoc[i], vector.facScale, diffexp.path)
#     if (i == 1) {
#         df.diffexp <- diffexp.list[[i]]
#     } else{
#         df.diffexp <- rbind(df.diffexp, diffexp.list[[i]])
#     }
# }

# PCA and umap plot
for (i in 1:length(diffexp.list)) {
    if (i == 1) {
        df.diffexp <- diffexp.list[[i]][['access.res']]
    } else {
        df.diffexp <- rbind(df.diffexp, diffexp.list[[i]][['access.res']])
    }
    ggplot.pca <- DimPlot(
        diffexp.list[[i]][['object.pca']],
        reductions = 'pca',
        group.by = "group",
        shape.by = 'batch',
        pt.size = 1.5
    )
    ggsave(
        plot = ggplot.pca, path = diffexp.path,
        filename = paste0('PCA_plot_diffexp_', de.facLoc[i], '.png'),
        units = 'cm', width = 15, height = 8
    )
    ggplot.umap <- DimPlot(
        diffexp.list[[i]][['object.pca']], reductions = 'umap', 
        group.by = "group", shape.by = 'batch', pt.size = 1.5
    )
    ggsave(
        plot = ggplot.umap, path = diffexp.path,
        filename = paste0('UMAP_plot_diffexp_', de.facLoc[i], '.png'),
        units = 'cm', width = 15, height = 8
    )
}

# save results
path.write <- paste(diffexp.path, 'diffexp_out.txt', sep = '/')
write.table(df.diffexp, file = path.write, quote = F, sep = '\t')

df.diffexp <- read.table(path.write, sep = '\t')

# plot
diffexp.names <- names(df.diffexp)
for (j in 1:dim(df.diffexp)[2]) {
    if (j == 1) {
        df.diffexp.plot <- data.frame(
            DiffexpParameter = de.facLoc, BatchEffectIndex = df.diffexp[,j],
            Label = rep(diffexp.names[j], length(de.facLoc)))
    } else {
        df.diffexp.plot <- rbind(
            df.diffexp.plot, 
            data.frame(
                DiffexpParameter = de.facLoc, 
                BatchEffectIndex = df.diffexp[,j],
                Label = rep(diffexp.names[j], length(de.facLoc))))
    }
}
diffexp.plot <- 
    ggplot(
    df.diffexp.plot, 
    aes(x = DiffexpParameter, y = BatchEffectIndex, 
        color = Label, shape = Label)) + 
    geom_line() + 
    geom_point(size = 2) + 
    labs(title = "Evaluate multiple methods for different degree of different expression") +
    scale_colour_discrete(name = 'Evaluation methods', 
                        breaks = c('pcReg', 'sil', 'kBET', 'mixent', 
                                   'ARI', 'NMI', 'ldaReg', 'sd'), 
                        labels = c('PC Regression', 'Silhouettes', 'kBET',
                                   'Entropy of batch mixing', 
                                   'Adjusted rand index', 
                                   'Mormalized mutual information',
                                   'LDA Regression', 'Standard Deviation')) + 
    scale_shape_manual(values = c(0, 1, 2, 4, 20, 5, 7, 24),
                       name = 'Evaluation methods', 
                       breaks = c('pcReg', 'sil', 'kBET', 'mixent', 
                                  'ARI', 'NMI', 'ldaReg', 'sd'), 
                       labels = c('PC Regression', 'Silhouettes', 'kBET',
                                  'Entropy of batch mixing', 
                                  'Adjusted rand index', 
                                  'Mormalized mutual information',
                                  'LDA Regression', 'Standard Deviation'))
ggsave(
    plot = diffexp.plot, path = diffexp.path, filename = "diffexp_out.png",
    units = 'cm', width = 30, height = 15)

# effective intervel and correlation coefficient
get.intervel <- function(df.method) {
    x <- df.method[, 1]
    y <- df.method[, 2]
    len <- length(x)
    step <- x[2] - x[1]
    diff <- c()
    for (i in 1:len) {
        if (i == 1) {
            diff <- c(diff, (y[i + 1] - y[i])/(x[i + 1] - x[i]))
        }
        if (i == len - 1) {
            diff <- c(diff, (y[i] - y[i - 1])/(x[i] - x[i - 1]))
        }
        diff <- c(diff, (y[i + 1] - y[i - 1])/(x[i + 1] - x[i - 1]))
    }
    effective.point <- c()
    for (i in 1:(len - 4)) {
        num.point <- sum(diff[seq(i, i + 4)] > 1)
        if (num.point >= 4) {
            effective.point <- c(effective.point, i)
        }
    }
    if (is.null(effective.point)) {
        return(list(start.interval = NA, 
                    end.interval = NA))
    }
    
    start.interval <- max(effective.point[1], 1) * step
    end.interval <- min(tail(effective.point, 1) + 4, len) * step
    
    return(list(start.interval = start.interval, end.interval = end.interval))
}

df.eval <- data.frame()
df.eval.plot <- data.frame()
for (method in unique(df.diffexp.plot$Label)) {
    df.method <- df.diffexp.plot[df.diffexp.plot$Label == method,]
    intervel <- get.intervel(df.method)
    start.interval <- intervel$start.interval
    end.interval <- intervel$end.interval
    if (is.na(start.interval) || is.na(end.interval)) {
        df.eval <- rbind(df.eval, 
                         data.frame(method = method, 
                                    start = start.interval, end = end.interval,
                                    pcc = NA, scc = NA, kcc = NA))
        df.eval.plot <- 
            rbind(
                df.eval.plot, 
                data.frame(method = method, point = start.interval, 
                           PCC = NA, SCC = NA, KCC = NA))
        df.eval.plot <- 
            rbind(
                df.eval.plot, 
                data.frame(method = method, point = end.interval, 
                           PCC = NA, SCC = NA, KCC = NA))
        next()
    }
    diffexp.interval <- 
        intersect(df.method$DiffexpParameter[
            df.method$DiffexpParameter > start.interval],
            df.method$DiffexpParameter[
                df.method$DiffexpParameter < end.interval])
    index.interval <- 
        df.method$BatchEffectIndex[
            df.method$DiffexpParameter %in% diffexp.interval]
    pcc = cor(diffexp.interval, index.interval, method = 'pearson')
    scc = cor(diffexp.interval, index.interval, method = 'spearman')
    kcc = cor(diffexp.interval, index.interval, method = 'kendall')
    df.eval <- rbind(df.eval, 
                     data.frame(method = method, 
                                start = start.interval, end = end.interval,
                                pcc = pcc, scc = scc, kcc = kcc))
    df.eval.plot <- 
        rbind(
            df.eval.plot, 
            data.frame(method = method, point = start.interval, 
                       PCC = pcc, SCC = scc, KCC = kcc))
    df.eval.plot <- 
        rbind(
            df.eval.plot, 
            data.frame(method = method, point = end.interval, 
                       PCC = pcc, SCC = scc, KCC = kcc))
}

eval.plot <-
    ggplot(df.eval.plot,
           aes(
               x = point,
               y = KCC,
               color = method,
               shape = method
           )) +
    geom_line() +
    geom_point(size = 2) +
    labs(title = "Evaluate multiple methods for different degree of different expression",
         x = 'Effective intervel', y = 'KCC') +
    scale_colour_discrete(
        name = 'Evaluation methods',
        breaks = c('pcReg', 'sil', 'kBET', 'mixent',
                   'ARI', 'NMI', 'ldaReg', 'sd'),
        labels = c(
            'PC Regression',
            'Silhouettes',
            'kBET',
            'Entropy of batch mixing',
            'Adjusted rand index',
            'Mormalized mutual information',
            'LDA Regression',
            'Standard Deviation'
        )
    ) +
    scale_shape_manual(
        values = c(0, 1, 2, 4, 20, 5, 7, 24),
        name = 'Evaluation methods',
        breaks = c('pcReg', 'sil', 'kBET', 'mixent',
                   'ARI', 'NMI', 'ldaReg', 'sd'),
        labels = c(
            'PC Regression',
            'Silhouettes',
            'kBET',
            'Entropy of batch mixing',
            'Adjusted rand index',
            'Mormalized mutual information',
            'LDA Regression',
            'Standard Deviation'
        )
    )
ggsave(
    plot = eval.plot, path = diffexp.path, filename = "interval.png",
    units = 'cm', width = 30, height = 15)

# # interval
# interval <- c()
# for (i in 1:dim(df.diffexp.plot)[1]) {
#     if (df.diffexp.plot$DiffexpParameter[i] >= 0 && 
#         df.diffexp.plot$DiffexpParameter[i] < 0.05) {
#         interval <- c(interval, '0-0.05')
#     }
#     if (df.diffexp.plot$DiffexpParameter[i] >= 0.05 && 
#         df.diffexp.plot$DiffexpParameter[i] < 0.1) {
#         interval <- c(interval, '0.05-0.1')
#     }
#     if (df.diffexp.plot$DiffexpParameter[i] >= 0.1 && 
#         df.diffexp.plot$DiffexpParameter[i] < 0.15) {
#         interval <- c(interval, '0.1-0.15')
#     }
#     if (df.diffexp.plot$DiffexpParameter[i] >= 0.15 && 
#         df.diffexp.plot$DiffexpParameter[i] <= 0.2) {
#         interval <- c(interval, '0.15-0.2')
#     }
#     
# }
# df.diffexp.plot$Interval <- interval
# df.barplot <- data.frame()
# for (method in unique(df.diffexp.plot$Label)) {
#     subdata_m <- df.diffexp.plot[df.diffexp.plot$Label == method,]
#     for (range in unique(df.diffexp.plot$Interval)) {
#         subdata <- subdata_m[df.diffexp.plot$Interval == range,]
#         pcc <- cor(subdata$DiffexpParameter, subdata$BatchEffectIndex, 
#                    method = 'pearson', use = 'complete.obs')
#         scc <- cor(subdata$DiffexpParameter, subdata$BatchEffectIndex, 
#                    method = 'spearman', use = 'complete.obs')
#         kcc <- cor(subdata$DiffexpParameter, subdata$BatchEffectIndex, 
#                    method = 'kendall', use = 'complete.obs')
#         df.barplot <- 
#             rbind(df.barplot,
#                   data.frame(Method = method, Interval = range,
#                              PCC = pcc, SCC = scc, KCC = kcc))
#     }
#     pcc <- cor(subdata_m$DiffexpParameter, subdata_m$BatchEffectIndex, 
#                method = 'pearson', use = 'complete.obs')
#     scc <- cor(subdata_m$DiffexpParameter, subdata_m$BatchEffectIndex, 
#                method = 'spearman', use = 'complete.obs')
#     kcc <- cor(subdata_m$DiffexpParameter, subdata_m$BatchEffectIndex, 
#                method = 'kendall', use = 'complete.obs')
#     df.barplot <- 
#         rbind(df.barplot,
#               data.frame(Method = method, Interval = '0-0.2',
#                          PCC = pcc, SCC = scc, KCC = kcc))
# }
# df.barplot[is.na(df.barplot)] <- 0
# 
# bar.plot <- 
#     ggplot(
#         df.barplot, 
#         aes(x = Interval, y = KCC, fill = Method)) + 
#     geom_bar(position = 'dodge', stat = 'identity') + 
#     labs(title = "Evaluate multiple methods for different degree of different expression",
#          y = 'Kendall Correlation Coefficient') +
#     scale_fill_discrete(name = 'Evaluation methods', 
#                           breaks = c('pcReg', 'sil', 'kBET', 'mixent', 
#                                      'ARI', 'NMI', 'ldaReg', 'sd'), 
#                           labels = c('PC Regression', 'Silhouettes', 'kBET',
#                                      'Entropy of batch mixing', 
#                                      'Adjusted rand index', 
#                                      'Mormalized mutual information',
#                                      'LDA Regression', 'Standard Deviation'))
# 
# ggsave(
#     plot = bar.plot, path = diffexp.path, filename = "barplot_KCC.png",
#     units = 'cm', width = 25, height = 15)
# 
# bar.plot <- 
#     ggplot(
#         df.barplot, 
#         aes(x = Interval, y = SCC, fill = Method)) + 
#     geom_bar(position = 'dodge', stat = 'identity') + 
#     labs(title = "Evaluate multiple methods for different degree of different expression",
#          y = 'Spearman Correlation Coefficient') +
#     scale_fill_discrete(name = 'Evaluation methods', 
#                         breaks = c('pcReg', 'sil', 'kBET', 'mixent', 
#                                    'ARI', 'NMI', 'ldaReg', 'sd'), 
#                         labels = c('PC Regression', 'Silhouettes', 'kBET',
#                                    'Entropy of batch mixing', 
#                                    'Adjusted rand index', 
#                                    'Mormalized mutual information',
#                                    'LDA Regression', 'Standard Deviation'))
# 
# ggsave(
#     plot = bar.plot, path = diffexp.path, filename = "barplot_SCC.png",
#     units = 'cm', width = 25, height = 15)
# 



