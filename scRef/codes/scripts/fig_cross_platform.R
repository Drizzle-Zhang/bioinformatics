library(ggplot2)

ref.dataset <- c('Tasic', 'celseq2', 'celseq2', '10Xv2')
# ref.GSE_id <- c('GSE71585', 'GSE85241', 'GSE85241', 'GSE132044')
ref.GSE_id <- c('Tasic', 'Muraro', 'Muraro', 'Ding')
ref.platform <- c('SMART-seq', 'CEL-seq2', 'CEL-seq2', '10x Chromium (v2)')

sc.dataset <- c('Campbell_equal', 'indrop', 'smartseq2', 'indrop')
# sc.GSE_id <- c('GSE93374', 'GSE84133', 'E-MTAB-5061', 'GSE132044')
sc.GSE_id <- c('Campbell', 'Baron', 'Segerstolpe', 'Ding')
sc.platform <- c('Drop-seq', 'inDrop', 'SMART-seq2', 'inDrop')

folders <- c('mouse_brain', 'human_panc', 'human_panc', 'human_pbmsca')
tissues <- c('Mouse brain', 'Human pancreas', 'Human pancreas', 'Human PBMC')

path.res <- '/home/zy/scRef/Benchmark/'
path.fig <- '/home/zy/scRef/figure/cross_platform/'

# scatter plot
# accuracy
df.plot <- data.frame(stringsAsFactors = F)
for (i in 1:length(ref.dataset)) {
    ref.name <- ref.dataset[i]
    sc.name <- sc.dataset[i]
    ref.plat <- ref.platform[i]
    sc.plat <- sc.platform[i]
    ref.GSE <- ref.GSE_id[i]
    sc.GSE <- sc.GSE_id[i]
    folder <- folders[i]
    tissue <- tissues[i]
    sub.res <- read.delim(paste0(path.res, folder, '/results_', ref.name, '_', sc.name, '.txt'), 
                          stringsAsFactors = F)
    sub.res <- sub.res[sub.res$term == 'Accuracy',]
    acc.scRef <- sub.res[sub.res$method == 'scRef', 'value']
    legend <- paste0(ref.GSE, '(', ref.plat, ')', ' -> ', sc.GSE, '(', sc.plat, ')', ', ', tissue)
    sub.plot <- data.frame(other = sub.res[sub.res$method != 'scRef', 'value'], 
                           scRef = rep(acc.scRef, nrow(sub.res)-1), 
                           method = sub.res[sub.res$method != 'scRef', 'method'],
                           dataset = rep(legend, nrow(sub.res)-1))
    df.plot <- rbind(df.plot, sub.plot)
}
plot.scatter <- 
    ggplot(data = df.plot, aes(x = other, y = scRef, color = method, shape = dataset)) + 
    geom_point(size = 3) + 
    geom_abline(intercept = 0, slope = 1, linetype = 2) + 
    ylim(0, 1) + xlim(0, 1) + 
    labs(x = 'Accuracy(Other method)', y = 'Accuracy(scMAGIC)', 
         shape = 'Query dataset', color = 'Method') + 
    theme(panel.background = element_rect(fill = 'transparent', color = 'gray'),
          legend.text = element_text(size = 7))
ggsave(filename = paste0('scatter_crossplat_acc.png'), 
       path = path.fig, plot = plot.scatter,
       units = 'cm', height = 12, width = 22)

# Macro F1
df.plot <- data.frame(stringsAsFactors = F)
for (i in 1:length(ref.dataset)) {
    ref.name <- ref.dataset[i]
    sc.name <- sc.dataset[i]
    ref.plat <- ref.platform[i]
    sc.plat <- sc.platform[i]
    ref.GSE <- ref.GSE_id[i]
    sc.GSE <- sc.GSE_id[i]
    folder <- folders[i]
    tissue <- tissues[i]
    sub.res <- read.delim(paste0(path.res, folder, '/results_', ref.name, '_', sc.name, '.txt'), 
                          stringsAsFactors = F)
    sub.res <- sub.res[sub.res$term == 'Macro F1',]
    acc.scRef <- sub.res[sub.res$method == 'scRef', 'value']
    legend <- paste0(ref.GSE, '(', ref.plat, ')', ' -> ', sc.GSE, '(', sc.plat, ')', ', ', tissue)
    sub.plot <- data.frame(other = sub.res[sub.res$method != 'scRef', 'value'], 
                           scRef = rep(acc.scRef, nrow(sub.res)-1), 
                           method = sub.res[sub.res$method != 'scRef', 'method'],
                           dataset = rep(legend, nrow(sub.res)-1))
    df.plot <- rbind(df.plot, sub.plot)
}
plot.scatter <- 
    ggplot(data = df.plot, aes(x = other, y = scRef, color = method, shape = dataset)) + 
    geom_point(size = 3) + 
    geom_abline(intercept = 0, slope = 1, linetype = 2) + 
    ylim(0, 1) + xlim(0, 1) + 
    labs(x = 'Macro F1(Other method)', y = 'Macro F1(scMAGIC)', 
         shape = 'Query dataset', color = 'Method') + 
    theme(panel.background = element_rect(fill = 'transparent', color = 'gray'),
          legend.text = element_text(size = 7))
ggsave(filename = paste0('scatter_crossplat_MacroF1.png'), 
       path = path.fig, plot = plot.scatter,
       units = 'cm', height = 12, width = 22)


