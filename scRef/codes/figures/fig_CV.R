library(ggplot2)
require("RColorBrewer")
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

df.plot.CV$`Mean (Cross validation)` <- rowMeans(df.plot.CV)


# cross platform
ref.dataset <- c('Tasic', 'smartseq2', 'CELSeq2', 'celseq2', 
                 'Dropseq', '10Xv2', 'Campbell', 'Campbell')
# ref.GSE_id <- c('GSE71585', 'GSE85241', 'GSE85241', 'GSE132044')
ref.GSE_id <- c('Tasic', 'Segerstolpe', 'Ding', 'Muraro', 
                'Ding', 'Ding', 'Campbell', 'Campbell')
ref.platform <- c('SMART-seq', 'SMART-seq2', 'CEL-seq2', 'CEL-seq2',
                  'Drop-seq', '10x Chromium (v2)', 'Drop-seq', 'Drop-seq')

sc.dataset <- c('Campbell_equal', 'celseq2', '10Xv2', 'indrop', 
                'CELSeq2', 'Smartseq2', 'Zeisel', 'Tasic')
# sc.GSE_id <- c('GSE93374', 'GSE84133', 'E-MTAB-5061', 'GSE132044')
sc.GSE_id <- c('Campbell', 'Muraro', 'Ding', 'Baron', 
               'Ding', 'Ding', 'Zeisel', 'Tasic')
# sc.platform <- c('Drop-seq', 'inDrop', 'CEL-seq2', '10x Chromium (v2)')
sc.platform <- c('Drop-seq', 'CEL-seq2', '10x Chromium (v2)', 'inDrop',
                 'CEL-seq2', 'SMART-seq2', 'Fluidigm C1', 'SMART-seq')

folders <- c('mouse_brain', 'human_panc', 'human_pbmsca', 'human_panc',
             'human_pbmsca', 'human_pbmsca', 'mouse_brain', 'mouse_brain')
tissues <- c('Mouse brain', 'Human pancreas', 'Human PBMC', 'Human pancreas', 
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

df.plot.CP$`Mean (Cross platform)` <- rowMeans(df.plot.CP)

mat.plot <- cbind(df.plot.CV, df.plot.CP)
mat.plot$`Mean Accuracy` <- (mat.plot$`Mean (Cross validation)` +
                                mat.plot$`Mean (Cross platform)`)/2
df.plot <- reshape2::melt(as.matrix(mat.plot))
names(df.plot) <- c('Method', 'Dataset', 'Accuracy')
vec.metrics <- c(rep('Cross validation', 50), rep('Cross platform', 90),
                 rep('Mean', 10))
df.plot$Type <- factor(vec.metrics, levels = c('Cross validation', 'Cross platform',
                                               'Mean'))
df.mean <- df.plot[df.plot$Dataset == 'Mean Accuracy', ]
df.plot$Method <- factor(df.plot$Method,
                            levels = df.mean$Method[order(df.mean$Accuracy, decreasing = F)])
color_type <- rep('0', nrow(df.plot))
sel.type <- c('Mean (Cross validation)', 'Mean (Cross platform)', 'Mean Accuracy')
color_type[(df.plot$Dataset %in% sel.type) & (df.plot$Accuracy < 0.1)] <- '1'
color_type[(df.plot$Dataset %in% sel.type) & (df.plot$Accuracy > 0.85)] <- '1'
color_type[(df.plot$Dataset %in% sel.type) & (df.plot$Accuracy > 0.1) &
               (df.plot$Accuracy < 0.85)] <- '2'
df.plot$Color_Type <- color_type

plot.heatmap <-
    ggplot(data = df.plot, aes(x = Dataset, y = Method)) +
    geom_tile(aes(fill = Accuracy)) +
    facet_grid( ~ Type, scales = 'free', space = 'free') +
    scale_fill_gradientn(
        colors = colorRampPalette(c(rev(brewer.pal(n = 9, name =  "RdYlBu"))[1:8], 'firebrick'))(100)) +
    scale_y_discrete(position = 'right') +
    labs(x = '', y = '', fill = 'Accuracy') +
    geom_text(aes(label = round(Accuracy, 2), color = Color_Type),
              family = "Arial", size = 3.8) +
    scale_color_manual(breaks = c('0', '1', '2'),
                       values = c('transparent', 'white', 'black')) +
    theme(panel.background = element_rect(color = 'transparent',
                                          fill = 'transparent'),
          axis.text.x = element_text(size = 13, color = "black", family = 'Arial',
                                     angle = 315, vjust = 1, hjust = 0),
          axis.text.y = element_text(
              size = 13, color = "black", family = 'Arial',
              vjust = 0.5, hjust = 0.5),
          strip.text.x = element_text(
              size = 12.8, color = "black", family = 'Arial'),
          strip.background = element_rect(
              color = "transparent", fill = "transparent"),
          legend.title = element_text(
              size = 12, color = "black", family = 'Arial'),
          legend.position = 'left') +
    guides(color = F)

ggsave(filename = 'heatmap_accuracy.png',
       path = path.fig, plot = plot.heatmap,
       units = 'cm', height = 18, width = 17)



# boxplot
df.box <- df.plot[df.plot$Type %in% c('Cross validation', 'Cross platform'),]

plot.box <-
    ggplot(df.box, aes(x = Method, y = Accuracy, color = Type)) +
    geom_boxplot(width = 0.7, outlier.size = 1) +
    # scale_fill_manual(breaks = c('Cross validation', 'Cross platform'),
    #                   values = c('dimgray', 'gainsboro')) +
    scale_color_manual(breaks = c('Cross validation', 'Cross platform'),
                       values = c("#1F78B4", "#EE7700")) +
    coord_flip() +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    labs(title = "", y = 'Accuracy', x = '', color = '', fill = '') +
    theme_bw() +
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          axis.text.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = 'lightgray',
                                            size = 0.15, linetype = 'dashed'),
          legend.position = 'bottom',
          axis.text.x = element_text(size = 9, color = 'black', family = 'Arial'),
          axis.title = element_text(size = 10, family = 'Arial'),
          legend.text = element_text(size = 9, family = 'Arial'),
          # legend.key.size = unit(1, 'cm'),
          legend.key = element_blank()) +
    guides(color = guide_legend(ncol = 1))
ggsave(filename = 'boxplot_Accuracy.png',
       path = path.fig, plot = plot.box,
       units = 'cm', height = 10.5, width = 5)


# Balanced Accuracy
datasets <- c('Campbell', 'panc8_indrop', 'pbmcsca_10Xv2', 'Tasic')
GSE_ids <- c('Campbell, mHArc-ME',
             'Baron, Human pancreas',
             'Ding, hPBMC', 'Tasic, mPVC')
path.output <- '/mdshare/node9/zy/scRef/cross_validation/train4_test1/'

methods <- c("scMAGIC", "SingleR", "scmap-cell", "scmap-cluster", "CHETAH", "scPred",
             "sciBet", "singleCellNet", "scID", "scClassify")

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

df.plot.CV$`Mean (Cross validation)` <- rowMeans(df.plot.CV)


# cross platform
ref.dataset <- c('Tasic', 'smartseq2', 'CELSeq2', 'celseq2', 
                 'Dropseq', '10Xv2', 'Campbell', 'Campbell')
# ref.GSE_id <- c('GSE71585', 'GSE85241', 'GSE85241', 'GSE132044')
ref.GSE_id <- c('Tasic', 'Segerstolpe', 'Ding', 'Muraro', 
                'Ding', 'Ding', 'Campbell', 'Campbell')
ref.platform <- c('SMART-seq', 'SMART-seq2', 'CEL-seq2', 'CEL-seq2',
                  'Drop-seq', '10x Chromium (v2)', 'Drop-seq', 'Drop-seq')

sc.dataset <- c('Campbell_equal', 'celseq2', '10Xv2', 'indrop', 
                'CELSeq2', 'Smartseq2', 'Zeisel', 'Tasic')
# sc.GSE_id <- c('GSE93374', 'GSE84133', 'E-MTAB-5061', 'GSE132044')
sc.GSE_id <- c('Campbell', 'Muraro', 'Ding', 'Baron', 
               'Ding', 'Ding', 'Zeisel', 'Tasic')
# sc.platform <- c('Drop-seq', 'inDrop', 'CEL-seq2', '10x Chromium (v2)')
sc.platform <- c('Drop-seq', 'CEL-seq2', '10x Chromium (v2)', 'inDrop',
                 'CEL-seq2', 'SMART-seq2', 'Fluidigm C1', 'SMART-seq')

folders <- c('mouse_brain', 'human_panc', 'human_pbmsca', 'human_panc',
             'human_pbmsca', 'human_pbmsca', 'mouse_brain', 'mouse_brain')
tissues <- c('Mouse brain', 'Human pancreas', 'Human PBMC', 'Human pancreas', 
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

df.plot.CP$`Mean (Cross platform)` <- rowMeans(df.plot.CP)

mat.plot <- cbind(df.plot.CV, df.plot.CP)
mat.plot$`Mean Balanced Accuracy` <- (mat.plot$`Mean (Cross validation)` +
                                 mat.plot$`Mean (Cross platform)`)/2
df.plot <- reshape2::melt(as.matrix(mat.plot))
names(df.plot) <- c('Method', 'Dataset', 'Accuracy')
vec.metrics <- c(rep('Cross validation', 50), rep('Cross platform', 90),
                 rep('Mean', 10))
df.plot$Type <- factor(vec.metrics, levels = c('Cross validation', 'Cross platform',
                                               'Mean'))
df.mean <- df.plot[df.plot$Dataset == 'Mean Balanced Accuracy', ]
df.plot$Method <- factor(df.plot$Method,
                         levels = df.mean$Method[order(df.mean$Accuracy, decreasing = F)])
color_type <- rep('0', nrow(df.plot))
sel.type <- c('Mean (Cross validation)', 'Mean (Cross platform)', 'Mean Balanced Accuracy')
color_type[(df.plot$Dataset %in% sel.type) & (df.plot$Accuracy < 0.1)] <- '1'
color_type[(df.plot$Dataset %in% sel.type) & (df.plot$Accuracy > 0.85)] <- '1'
color_type[(df.plot$Dataset %in% sel.type) & (df.plot$Accuracy > 0.1) &
               (df.plot$Accuracy < 0.85)] <- '2'
df.plot$Color_Type <- color_type

plot.heatmap <-
    ggplot(data = df.plot, aes(x = Dataset, y = Method)) +
    geom_tile(aes(fill = Accuracy)) +
    facet_grid( ~ Type, scales = 'free', space = 'free') +
    scale_fill_gradientn(
        colors = colorRampPalette(c(rev(brewer.pal(n = 9, name =  "RdYlBu"))[1:8], 'firebrick'))(100)) +
    scale_y_discrete(position = 'right') +
    labs(x = '', y = '', fill = 'Balanced \nAccuracy') +
    geom_text(aes(label = round(Accuracy, 2), color = Color_Type),
              family = "Arial", size = 3.8) +
    scale_color_manual(breaks = c('0', '1', '2'),
                       values = c('transparent', 'white', 'black')) +
    theme(panel.background = element_rect(color = 'transparent',
                                          fill = 'transparent'),
          axis.text.x = element_text(size = 13, color = "black", family = 'Arial',
                                     angle = 315, vjust = 1, hjust = 0),
          axis.text.y = element_text(
              size = 13, color = "black", family = 'Arial',
              vjust = 0.5, hjust = 0.5),
          strip.text.x = element_text(
              size = 12.8, color = "black", family = 'Arial'),
          strip.background = element_rect(
              color = "transparent", fill = "transparent"),
          legend.title = element_text(
              size = 12, color = "black", family = 'Arial'),
          legend.position = 'left') +
    guides(color = F)

ggsave(filename = 'heatmap_Balanced_Accuracy.png',
       path = path.fig, plot = plot.heatmap,
       units = 'cm', height = 18, width = 17)



# boxplot
df.box <- df.plot[df.plot$Type %in% c('Cross validation', 'Cross platform'),]

plot.box <-
    ggplot(df.box, aes(x = Method, y = Accuracy, color = Type)) +
    geom_boxplot(width = 0.7, outlier.size = 1) +
    scale_color_manual(breaks = c('Cross validation', 'Cross platform'),
                       values = c("#1F78B4", "#EE7700")) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    coord_flip() +
    labs(title = "", y = 'Balanced Accuracy', x = '', color = '', fill = '') +
    theme_bw() +
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          axis.text.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = 'lightgray',
                                            size = 0.15, linetype = 'dashed'),
          legend.position = 'bottom',
          axis.text.x = element_text(size = 9, color = 'black', family = 'Arial'),
          axis.title = element_text(size = 10, family = 'Arial'),
          legend.text = element_text(size = 9, family = 'Arial'),
          # legend.key.size = unit(1, 'cm'),
          legend.key = element_blank()) +
    guides(color = guide_legend(ncol = 1))
ggsave(filename = 'boxplot_Balanced_Accuracy.png',
       path = path.fig, plot = plot.box,
       units = 'cm', height = 10.5, width = 5)







