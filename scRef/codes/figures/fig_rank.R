library(ggplot2)
library(RColorBrewer)

path.fig <- '/mdshare/node9/zy/scRef/figure/Rank/'
if (!file.exists(path.fig)) {
    dir.create(path.fig)
}


methods <- c("scMAGIC", "SingleR", "scmap-cell", "scmap-cluster", "CHETAH", "scPred",
             "sciBet", "singleCellNet", "scID", "scClassify")

df_plot <- data.frame(stringsAsFactors = F)
# cross validation
# Accuracy
datasets <- c('Campbell', 'panc8_indrop', 'pbmcsca_10Xv2', 'Tasic')
# GSE_ids <- c('GSE93374', 'GSE84133', 'GSE132044', 'GSE71585')
# GSE_ids <- c('Campbell, Mouse hypothalamic Arc-ME',
#              'Baron, Human pancreas',
#              'Ding, PBMC', 'Tasic, Primary visual cortex')
GSE_ids <- c('Campbell, Mouse HArc-ME',
             'Baron, Human pancreas',
             'Ding, Human PBMC',
             'Tasic, Mouse PVC')
path.output <- '/mdshare/node9/zy/scRef/cross_validation/train4_test1/'

df.plot.CV <- data.frame(stringsAsFactors = F)
for (i in 1:length(datasets)) {
    dataset <- datasets[i]
    sub.res <- read.delim(paste0(path.output, dataset, '/results_', dataset, '.txt'), stringsAsFactors = F)
    sub.acc <- sub.res[sub.res$term == 'Accuracy',]
    sub.acc$type <- rep('cross validation', nrow(sub.acc))
    sub.acc$rank <- nrow(sub.acc) - ceiling(rank(sub.acc$value)) + 1
    sub.bacc <- sub.res[sub.res$term == 'Balanced accuracy',]
    sub.bacc$type <- rep('cross validation', nrow(sub.bacc))
    sub.bacc$rank <- nrow(sub.acc) - ceiling(rank(sub.bacc$value)) + 1
    sub.plot <- rbind(sub.acc, sub.bacc)
    sub.plot$dataset <- rep(GSE_ids[i], nrow(sub.plot))
    sub.plot$data_idx <- rep(paste0('A', i), nrow(sub.plot))
    df.plot.CV <- rbind(df.plot.CV, sub.plot)
}
df_plot <- rbind(df_plot, df.plot.CV)


# cross platform
ref.dataset <- c('smartseq2', 'celseq2', 
                 'Dropseq', '10Xv2', 'Campbell', 'Campbell')
ref.GSE_id <- c('Segerstolpe', 'Muraro', 
                'Ding', 'Ding', 'Campbell', 'Campbell')
ref.platform <- c('SMART-seq2', 'CEL-seq2',
                  'Drop-seq', '10x Chromium (v2)', 'Drop-seq', 'Drop-seq')

sc.dataset <- c('celseq2', 'indrop', 
                'CELSeq2', 'Smartseq2', 'Zeisel', 'Tasic')
sc.GSE_id <- c('Muraro', 'Baron', 
               'Ding', 'Ding', 'Zeisel', 'Tasic')
sc.platform <- c('CEL-seq2', 'inDrop',
                 'CEL-seq2', 'SMART-seq2', 'Fluidigm C1', 'SMART-seq')

folders <- c('human_panc', 'human_panc',
             'human_pbmsca', 'human_pbmsca', 'mouse_brain', 'mouse_brain')
tissues <- c('Human pancreas', 'Human pancreas', 
             'Human PBMC', 'Human PBMC', 'Mouse brain', 'Mouse brain')

path.res <- '/mdshare/node9/zy/scRef/Benchmark/'

# accuracy
df.plot.CP <- data.frame(stringsAsFactors = F)
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
    sub.acc <- sub.res[sub.res$term == 'Accuracy',]
    sub.acc$type <- rep('cross platform', nrow(sub.acc))
    sub.acc$rank <- nrow(sub.acc) - ceiling(rank(sub.acc$value)) + 1
    sub.bacc <- sub.res[sub.res$term == 'Balanced accuracy',]
    sub.bacc$type <- rep('cross platform', nrow(sub.bacc))
    sub.bacc$rank <- nrow(sub.acc) - ceiling(rank(sub.bacc$value)) + 1
    sub.plot <- rbind(sub.acc, sub.bacc)
    sub.plot$dataset <- rep(paste0(ref.plat, ' -> ', sc.plat, ', ', tissue), nrow(sub.plot))
    sub.plot$data_idx <- rep(paste0('B', i), nrow(sub.plot))
    df.plot.CP <- rbind(df.plot.CP, sub.plot)
}
df_plot <- rbind(df_plot, df.plot.CP)


# Tasic
ref.data <- 'Tasic'
sc.dataset <- c('Tasic2018', 'Campbell', 'HochgernerA', 'Mizrak')
sc.legends <- c('Tasic, mouse Neocortex',
                'Campbell, mouse HArc-ME',
                'Hochgerner, mouse DGH',
                'Mizrak, mouse V-SVZ')

path.res <- '/mdshare/node9/zy/scRef/Benchmark/mouse_brain/'

df.plot.Tasic <- data.frame(stringsAsFactors = F)
for (i in 1:length(sc.dataset)) {
    sc.data <- sc.dataset[i]
    sc.legend <- sc.legends[i]
    file.res <- paste0(path.res, 'results_', ref.data, '_', sc.data, '.txt')
    sub.res <- read.delim(file = file.res, row.names = 1, stringsAsFactors = F)
    sub.acc <- sub.res[sub.res$term == 'Accuracy',]
    sub.acc$type <- rep('Tasic', nrow(sub.acc))
    sub.acc$rank <- nrow(sub.acc) - ceiling(rank(sub.acc$value)) + 1
    sub.bacc <- sub.res[sub.res$term == 'Balanced accuracy',]
    sub.bacc$type <- rep('Tasic', nrow(sub.bacc))
    sub.bacc$rank <- nrow(sub.acc) - ceiling(rank(sub.bacc$value)) + 1
    sub.plot <- rbind(sub.acc, sub.bacc)
    sub.plot$dataset <- rep(paste0('Tasic_', sc.legend), nrow(sub.plot))
    sub.plot$data_idx <- rep(paste0('C', i), nrow(sub.plot))
    df.plot.Tasic <- rbind(df.plot.Tasic, sub.plot)
}
df_plot <- rbind(df_plot, df.plot.Tasic)

# MCA
ref.data <- 'MCA'
sc.dataset <- c('Tasic2018', 'Campbell', 'HochgernerA', 'Mizrak')
sc.legends <- c('Tasic, mouse Neocortex',
                'Campbell, mouse HArc-ME',
                'Hochgerner, mouse DGH',
                'Mizrak, mouse V-SVZ')

df.plot.MCA <- data.frame(stringsAsFactors = F)
for (i in 1:length(sc.dataset)) {
    sc.data <- sc.dataset[i]
    sc.legend <- sc.legends[i]
    file.res <- paste0(path.res, 'results_', ref.data, '_', sc.data, '.txt')
    sub.res <- read.delim(file = file.res, row.names = 1, stringsAsFactors = F)
    sub.acc <- sub.res[sub.res$term == 'Accuracy',]
    sub.acc$type <- rep('MCA', nrow(sub.acc))
    sub.acc$rank <- nrow(sub.acc) - ceiling(rank(sub.acc$value)) + 1
    sub.bacc <- sub.res[sub.res$term == 'Balanced accuracy',]
    sub.bacc$type <- rep('MCA', nrow(sub.bacc))
    sub.bacc$rank <- nrow(sub.acc) - ceiling(rank(sub.bacc$value)) + 1
    sub.plot <- rbind(sub.acc, sub.bacc)
    sub.plot$dataset <- rep(paste0('MCA_', sc.legend), nrow(sub.plot))
    sub.plot$data_idx <- rep(paste0('D', i), nrow(sub.plot))
    df.plot.MCA <- rbind(df.plot.MCA, sub.plot)
}
df_plot <- rbind(df_plot, df.plot.MCA)


# cross species
ref.dataset <- c('BaronM', 'BaronM', 'MCA', 'MCA')
sc.dataset <- c('panc8_celseq2', 'panc8_smartseq2', 'pbmc3k', 'pbmcsca_10Xv2')
sc.legends <- c('Baron -> Muraro, Pancreas',
                'Baron -> Segerstolpe, Pancreas',
                'Han -> Butler, PBMC',
                'Han -> Ding, PBMC')

path.res <- '/mdshare/node9/zy/scRef/Benchmark/cross_species/'

df.plot.CS <- data.frame(stringsAsFactors = F)
for (i in 1:length(sc.dataset)) {
    ref.data <- ref.dataset[i]
    sc.data <- sc.dataset[i]
    sc.legend <- sc.legends[i]
    file.res <- paste0(path.res, 'results_', ref.data, '_', sc.data, '.txt')
    sub.res <- read.delim(file = file.res, row.names = 1, stringsAsFactors = F)
    sub.acc <- sub.res[sub.res$term == 'Accuracy',]
    sub.acc$type <- rep('cross species', nrow(sub.acc))
    sub.acc$rank <- nrow(sub.acc) - ceiling(rank(sub.acc$value)) + 1
    sub.bacc <- sub.res[sub.res$term == 'Balanced accuracy',]
    sub.bacc$type <- rep('cross species', nrow(sub.bacc))
    sub.bacc$rank <- nrow(sub.acc) - ceiling(rank(sub.bacc$value)) + 1
    sub.plot <- rbind(sub.acc, sub.bacc)
    sub.plot$dataset <- rep(paste0('CS_', sc.legend), nrow(sub.plot))
    sub.plot$data_idx <- rep(paste0('E', i), nrow(sub.plot))
    df.plot.CS <- rbind(df.plot.CS, sub.plot)
}
df_plot <- rbind(df_plot, df.plot.CS)

# plot
df_plot$method <- gsub('singleR', 'SingleR', df_plot$method)
vec.color <- brewer.pal(10,"Paired")
# vec.color <- brewer.pal(10,"Spectral")

# accuracy
df_plot_acc <- df_plot[df_plot$term == 'Accuracy',]
df_rank_acc <- aggregate(df_plot_acc$rank, by = list(df_plot_acc$method), mean)
colnames(df_rank_acc) <- c('method', 'mean')
df_rank_acc <- as.data.frame(df_rank_acc)
df_plot_acc$method <- factor(df_plot_acc$method, 
                             levels = df_rank_acc$method[order(df_rank_acc$mean, decreasing = T)])
plot.rank <- 
    ggplot(data = df_plot_acc, aes(x = data_idx, y = rank, 
                                   color = method, fill = method, 
                                   group = 1)) + 
    geom_area(alpha = 0.6) + 
    scale_color_manual(values = vec.color) + 
    scale_fill_manual(values = vec.color) + 
    facet_grid(method ~ ., scales = 'free', space = 'free') + 
    ylim(c(0, 10)) + 
    theme(panel.background = element_rect(color = 'transparent',
                                          fill = 'transparent'), 
          strip.text.y = element_text(
              size = 13, color = "black", family = 'Arial', angle = 0), 
          strip.background = element_rect(
              color = "transparent", fill = "transparent"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks.y = element_blank(), 
          legend.position = 'none')
ggsave(filename = 'Line_accuracy.png', 
       path = path.fig, plot = plot.rank,
       units = 'cm', height = 10, width = 23)

