library(ggplot2)

datasets <- c('Campbell', 'panc8_indrop', 'pbmcsca_10Xv2', 'Tasic')
# GSE_ids <- c('GSE93374', 'GSE84133', 'GSE132044', 'GSE71585')
GSE_ids <- c('Campbell, Mouse hypothalamic Arc-ME', 
             'Baron, Human pancreas', 
             'Ding, PBMC', 'Tasic, Primary visual cortex')
path.output <- '/home/zy/scRef/cross_validation/train4_test1/'
path.fig <- '/home/zy/scRef/figure/'

# barplot
df.plot <- data.frame(stringsAsFactors = F)
for (i in 1:length(datasets)) {
    dataset <- datasets[i]
    sub.res <- read.delim(paste0(path.output, dataset, '/results_', dataset, '.txt'), stringsAsFactors = F)
    sub.res <- sub.res[sub.res$term == 'Accuracy',]
    sub.res$dataset <- rep(GSE_ids[i], dim(sub.res)[1])
    df.plot <- rbind(df.plot, sub.res)
}

accuracy.mean <- aggregate(df.plot$value, by = list(df.plot$method), FUN = mean)
df.plot$method <- factor(df.plot$method, 
                            levels = accuracy.mean$Group.1[order(accuracy.mean$x, decreasing = T)])

plot.bar <- ggplot(df.plot,
                   aes(x = method, y = value, fill = dataset)) +
    geom_bar(position = 'dodge', stat = 'identity') +
    theme_bw() +
    # facet_wrap(~ Evaluator, scales = 'free_x', ncol = 2) +
    labs(title = "", y = 'Accuracy', x = 'DataSet', fill = 'Methods') + 
    theme(axis.text = element_text(size = 9),
          panel.grid = element_blank(),
          panel.grid.major.y = element_line(color = 'grey', size = 0.2),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(size = 12))
ggsave(filename = 'barplot_CV.png', 
       path = '/home/zy/scRef/figure/', plot = plot.bar,
       units = 'cm', height = 12, width = 18)


# scatter plot
# accuracy
df.plot <- data.frame(stringsAsFactors = F)
for (i in 1:length(datasets)) {
    dataset <- datasets[i]
    GSE_id <- GSE_ids[i]
    sub.res <- read.delim(paste0(path.output, dataset, '/results_', dataset, '.txt'), stringsAsFactors = F)
    sub.res <- sub.res[sub.res$term == 'Accuracy',]
    acc.scRef <- sub.res[sub.res$method == 'scRef', 'value']
    sub.plot <- data.frame(other = sub.res[sub.res$method != 'scRef', 'value'], 
                           scRef = rep(acc.scRef, nrow(sub.res)-1), 
                           method = sub.res[sub.res$method != 'scRef', 'method'],
                           dataset = rep(GSE_id, nrow(sub.res)-1))
    df.plot <- rbind(df.plot, sub.plot)
}
plot.scatter <- 
    ggplot(data = df.plot, aes(x = other, y = scRef, color = method, shape = dataset)) + 
    geom_point(size = 3) + 
    geom_abline(intercept = 0, slope = 1, linetype = 2) + 
    ylim(0, 1) + xlim(0, 1) + 
    labs(x = 'Accuracy(Other method)', y = 'Accuracy(scMAGIC)', 
         shape = 'Query dataset', color = 'Method') + 
    theme(panel.background = element_rect(fill = 'transparent', color = 'gray'))
ggsave(filename = paste0('scatter_CV.png'), 
       path = path.fig, plot = plot.scatter,
       units = 'cm', height = 12, width = 20)

# Macro F1
df.plot <- data.frame(stringsAsFactors = F)
for (i in 1:length(datasets)) {
    dataset <- datasets[i]
    GSE_id <- GSE_ids[i]
    sub.res <- read.delim(paste0(path.output, dataset, '/results_', dataset, '.txt'), stringsAsFactors = F)
    sub.res <- sub.res[sub.res$term == 'macro F1',]
    acc.scRef <- sub.res[sub.res$method == 'scRef', 'value']
    sub.plot <- data.frame(other = sub.res[sub.res$method != 'scRef', 'value'], 
                           scRef = rep(acc.scRef, nrow(sub.res)-1), 
                           method = sub.res[sub.res$method != 'scRef', 'method'],
                           dataset = rep(GSE_id, nrow(sub.res)-1))
    df.plot <- rbind(df.plot, sub.plot)
}
plot.scatter <- 
    ggplot(data = df.plot, aes(x = other, y = scRef, color = method, shape = dataset)) + 
    geom_point(size = 3) + 
    geom_abline(intercept = 0, slope = 1, linetype = 2) + 
    ylim(0, 1) + xlim(0, 1) + 
    labs(x = 'Macro F1(Other method)', y = 'Macro F1(scMAGIC)', 
         shape = 'Query dataset', color = 'Method') + 
    theme(panel.background = element_rect(fill = 'transparent', color = 'gray'))
ggsave(filename = paste0('scatter_MacroF1_CV.png'), 
       path = path.fig, plot = plot.scatter,
       units = 'cm', height = 12, width = 20)







