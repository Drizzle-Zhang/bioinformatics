library(ggplot2)
library(RColorBrewer)

# data
path.fig <- '/mdshare/node9/zy/scRef/figure/equal_celltype/'

# CV
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

methods <- c("scMAGIC", "SingleR", "scmap-cell", "scmap-cluster", "CHETAH", "scPred",
             "sciBet", "singleCellNet", "scID", "scClassify")

# Accuracy
df.plot.CV <- data.frame(stringsAsFactors = F)
for (i in 1:length(datasets)) {
    dataset <- datasets[i]
    sub.res <- read.delim(paste0(path.output, dataset, '/results_', dataset, '.txt'), stringsAsFactors = F)
    sub.res <- sub.res[sub.res$term == 'Accuracy',]
    rownames(sub.res) <- sub.res$method
    sub.plot <- data.frame(sub.res[methods, 'value'], row.names = methods)
    names(sub.plot) <- GSE_ids[i]
    if (i == 1) {
        df.plot.CV <- sub.plot
    } else {
        df.plot.CV <- cbind(df.plot.CV, sub.plot)
    }
}

df.plot.CV$`Mean Accuracy (Cross validation)` <- rowMeans(df.plot.CV)
df.plot.acc <- df.plot.CV


# Balanced accuracy
df.plot.CV <- data.frame(stringsAsFactors = F)
for (i in 1:length(datasets)) {
    dataset <- datasets[i]
    sub.res <- read.delim(paste0(path.output, dataset, '/results_', dataset, '.txt'), stringsAsFactors = F)
    sub.res <- sub.res[sub.res$term == 'Balanced accuracy',]
    rownames(sub.res) <- sub.res$method
    sub.plot <- data.frame(sub.res[methods, 'value'], row.names = methods)
    names(sub.plot) <- GSE_ids[i]
    if (i == 1) {
        df.plot.CV <- sub.plot
    } else {
        df.plot.CV <- cbind(df.plot.CV, sub.plot)
    }
}

df.plot.CV$`Mean Balanced Accuracy (Cross validation)` <- rowMeans(df.plot.CV)
df.plot.ban_acc <- df.plot.CV

mat.plot.CV <- cbind(df.plot.acc, df.plot.ban_acc)


# cross platform
ref.dataset <- c('smartseq2', 'celseq2', 
                 'Dropseq', '10Xv2', 'Campbell', 'Campbell')
# ref.GSE_id <- c('GSE71585', 'GSE85241', 'GSE85241', 'GSE132044')
ref.GSE_id <- c('Segerstolpe', 'Muraro', 
                'Ding', 'Ding', 'Campbell', 'Campbell')
ref.platform <- c('SMART-seq2', 'CEL-seq2',
                  'Drop-seq', '10x Chromium (v2)', 'Drop-seq', 'Drop-seq')

sc.dataset <- c('celseq2', 'indrop', 
                'CELSeq2', 'Smartseq2', 'Zeisel', 'Tasic')
# sc.GSE_id <- c('GSE93374', 'GSE84133', 'E-MTAB-5061', 'GSE132044')
sc.GSE_id <- c('Muraro', 'Baron', 
               'Ding', 'Ding', 'Zeisel', 'Tasic')
# sc.platform <- c('Drop-seq', 'inDrop', 'CEL-seq2', '10x Chromium (v2)')
sc.platform <- c('CEL-seq2', 'inDrop',
                 'CEL-seq2', 'SMART-seq2', 'Fluidigm C1', 'SMART-seq')

folders <- c('human_panc', 'human_panc',
             'human_pbmsca', 'human_pbmsca', 'mouse_brain', 'mouse_brain')
tissues <- c('Human pancreas', 'Human pancreas', 
             'Human PBMC', 'Human PBMC', 'Mouse brain', 'Mouse brain')

path.res <- '/mdshare/node9/zy/scRef/Benchmark/'

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
    sub.res <- sub.res[sub.res$term == 'Accuracy',]
    sub.res$method <- gsub("singleR", "SingleR", sub.res$method)
    rownames(sub.res) <- sub.res$method
    sub.plot <- data.frame(sub.res[methods, 'value'], row.names = methods)
    # names(sub.plot) <- paste0(ref.GSE, '(', ref.plat, ')', ' -> ', sc.GSE, '(', sc.plat, ')', ', ', tissue)
    names(sub.plot) <- paste0(ref.plat, ' -> ', sc.plat, ', ', tissue)
    if (i == 1) {
        df.plot.CP <- sub.plot
    } else {
        df.plot.CP <- cbind(df.plot.CP, sub.plot)
    }
}

df.plot.CP$`Mean Accuracy (Cross platform)` <- rowMeans(df.plot.CP)
df.plot.acc <- df.plot.CP

# Balanced accuracy
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
    sub.res <- sub.res[sub.res$term == 'Balanced accuracy',]
    sub.res$method <- gsub("singleR", "SingleR", sub.res$method)
    rownames(sub.res) <- sub.res$method
    sub.plot <- data.frame(sub.res[methods, 'value'], row.names = methods)
    # names(sub.plot) <- paste0(ref.GSE, '(', ref.plat, ')', ' -> ', sc.GSE, '(', sc.plat, ')', ', ', tissue)
    names(sub.plot) <- paste0(ref.plat, ' -> ', sc.plat, ', ', tissue)
    if (i == 1) {
        df.plot.CP <- sub.plot
    } else {
        df.plot.CP <- cbind(df.plot.CP, sub.plot)
    }
}

df.plot.CP$`Mean Balanced Accuracy (Cross platform)` <- rowMeans(df.plot.CP)
df.plot.ban_acc <- df.plot.CP

mat.plot.CP <- cbind(df.plot.acc, df.plot.ban_acc)
mat.plot <- cbind(mat.plot.CV, mat.plot.CP)
mat.plot$`Mean Accuracy` <- rowMeans(mat.plot[, c(5, 17)])
mat.plot$`Mean Balanced Accuracy` <- rowMeans(mat.plot[, c(10, 24)])
# mat.plot <- cbind(mat.plot[, 1:5], mat.plot[, 6:10], 
#                   mat.plot[, c(11:15)], mat.plot[, c(16:20)], mat.plot[, c(21:22)])
mat.plot <- mat.plot[order(mat.plot$`Mean Accuracy`, decreasing = T),]

df.plot <- reshape2::melt(as.matrix(mat.plot))
names(df.plot) <- c('Method', 'Dataset', 'Accuracy')
vec.metrics <- c(rep('CV_Accuracy', 50),  rep('CV_Balanced Accuracy', 50),  
                 rep('CP_Accuracy', 70), rep('CP_Balanced Accuracy', 70),
                 rep('Mean_Accuracy', 20))
df.plot$Type <- factor(vec.metrics, 
                       levels = c('CV_Accuracy', 'CV_Balanced Accuracy', 
                                  'CP_Accuracy', 'CP_Balanced Accuracy', 
                                  'Mean_Accuracy'))
df.mean <- df.plot[df.plot$Dataset == 'Mean Accuracy', ]
df.plot$Method <- factor(df.plot$Method, 
                         levels = rev(df.mean$Method[order(df.mean$Accuracy, decreasing = F)]))
sort.methods <- df.mean$Method[order(df.mean$Accuracy, decreasing = F)]
color_type <- rep('0', nrow(df.plot))
color_type[(df.plot$Dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy',
                                   'Mean Accuracy (Cross validation)',
                                   'Mean Balanced Accuracy (Cross validation)',
                                   'Mean Accuracy (Cross platform)',
                                   'Mean Balanced Accuracy (Cross platform)')) & 
               (df.plot$Accuracy < 0.1)] <- '1'
color_type[(df.plot$Dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy',
                                   'Mean Accuracy (Cross validation)',
                                   'Mean Balanced Accuracy (Cross validation)',
                                   'Mean Accuracy (Cross platform)',
                                   'Mean Balanced Accuracy (Cross platform)')) & 
               (df.plot$Accuracy > 0.85)] <- '1'
color_type[(df.plot$Dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy',
                                   'Mean Accuracy (Cross validation)',
                                   'Mean Balanced Accuracy (Cross validation)',
                                   'Mean Accuracy (Cross platform)',
                                   'Mean Balanced Accuracy (Cross platform)')) & 
               (df.plot$Accuracy > 0.1) & (df.plot$Accuracy < 0.85)] <- '2'
df.plot$Color_Type <- color_type

df.plot$Dataset[df.plot$Dataset %in% c('Mean Accuracy (Cross validation)',
                                       'Mean Accuracy (Cross platform)')] <- 'Mean Accuracy'
df.plot$Dataset[df.plot$Dataset %in% c('Mean Balanced Accuracy (Cross validation)',
                                       'Mean Balanced Accuracy (Cross platform)')] <- 'Mean Balanced Accuracy'

plot.heatmap <-
    ggplot(data = df.plot, aes(x = Method, y = Dataset)) +
    geom_tile(aes(fill = Accuracy)) +
    facet_grid(Type ~ ., scales = 'free', space = 'free') +
    scale_fill_gradientn(
        colors = colorRampPalette(c(rev(brewer.pal(n = 9, name =  "RdYlBu"))[1:8], 'firebrick'))(100)) +
    scale_y_discrete(position = 'right') +
    labs(x = '', y = '', fill = '') +
    geom_text(aes(label = round(Accuracy, 2), color = Color_Type),
              family = "Arial", size = 3.8) +
    scale_color_manual(breaks = c('0', '1', '2'),
                       values = c('transparent', 'white', 'black')) +
    theme(panel.background = element_rect(color = 'transparent',
                                          fill = 'transparent'),
          axis.text.x = element_blank(),
          axis.text.y = element_text(
              size = 13, color = "black", family = 'Arial',
              vjust = 0.5, hjust = 0.5),
          strip.text.y = element_blank(),
          strip.background = element_rect(
              color = "transparent", fill = "transparent"),
          legend.title = element_text(
              size = 12, color = "black", family = 'Arial'),
          legend.text = element_text(size = 11, family = 'Arial'),
          legend.position = 'none') +
    guides(color = F)

ggsave(filename = 'heatmap_CV_CP.png', 
       path = path.fig, plot = plot.heatmap,
       units = 'cm', height = 20, width = 22)

# boxplot
df.box <- df.plot[df.plot$Color_Type == '0',]
df.box$Metrics <- unlist(lapply(as.character(df.box$Type), function(x) {strsplit(x, '_')[[1]][2]}))

plot.box <- 
    ggplot(df.box, aes(x = Method, y = Accuracy, color = Metrics)) +
    geom_boxplot(width = 0.7, outlier.size = 1) +
    scale_color_manual(breaks = c('Accuracy', 'Balanced Accuracy'),
                       values = c("#1F78B4", "#EE7700")) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    # coord_flip() + 
    labs(title = "", y = 'Accuracy \n Balanced Accuracy', x = '', color = 'Metrics', fill = '') + 
    theme_bw() + 
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = 'lightgray',
                                            size = 0.15, linetype = 'dashed'),
          legend.position = 'right',
          axis.text.y = element_text(size = 11, color = 'black', family = 'Arial'),
          axis.text.x = element_text(size = 16, color = 'black', family = 'Arial',
                                     angle = 315, vjust = 1, hjust = 0),
          axis.title = element_text(size = 13, family = 'Arial'),
          legend.title = element_text(
              size = 14, color = "black", family = 'Arial'),
          legend.text = element_text(size = 13, family = 'Arial'),
          # legend.key.size = unit(1, 'cm'),
          legend.key = element_blank()) +
    guides(color = guide_legend(ncol = 1))
ggsave(filename = 'boxplot_CV_CP.png', 
       path = path.fig, plot = plot.box,
       units = 'cm', height = 9, width = 20)
