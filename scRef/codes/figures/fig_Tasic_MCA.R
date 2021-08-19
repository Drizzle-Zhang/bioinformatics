library(ggplot2)
library(RColorBrewer)

# data
path.fig <- '/mdshare/node9/zy/scRef/figure/Tasic_MCA/'
if (!file.exists(path.fig)) {
    dir.create(path.fig)
}

ref.data <- 'Tasic'
sc.dataset <- c('Tasic2018', 'Campbell', 'HochgernerA', 'Mizrak')
# sc.legends <- c('Tasic, Neocortex, SMART-Seq v4',
#                 'Campbell, hypothalamic Arc-ME, Drop-seq',
#                 'Hochgerner, Denatate gyrus, 10X',
#                 'Mizrak, V-SVZ, Drop-seq')
sc.legends <- c('Tasic, mouse Neocortex',
                'Campbell, mouse HArc-ME',
                'Hochgerner, mouse DGH',
                'Mizrak, mouse V-SVZ')

path.res <- '/mdshare/node9/zy/scRef/Benchmark/mouse_brain/'
# path.fig <- '/home/zy/scRef/figure/mouse_brain/'

methods <- c("scMAGIC", "SingleR", "scmap-cell", "scmap-cluster", "CHETAH", "scPred",
             "sciBet", "singleCellNet", "scID", "scClassify")

# Accuracy
df.plot.acc <- data.frame(stringsAsFactors = F)
for (i in 1:length(sc.dataset)) {
    sc.data <- sc.dataset[i]
    sc.legend <- sc.legends[i]
    file.res <- paste0(path.res, 'results_', ref.data, '_', sc.data, '.txt')
    sub.res <- read.delim(file = file.res, row.names = 1, stringsAsFactors = F)
    sub.res <- sub.res[sub.res$term == 'Accuracy',]
    sub.res$method <- gsub("singleR", "SingleR", sub.res$method)
    rownames(sub.res) <- sub.res$method
    sub.plot <- data.frame(sub.res[methods, 'value'], row.names = methods)
    names(sub.plot) <- sc.legend
    if (i == 1) {
        df.plot.acc <- sub.plot
    } else {
        df.plot.acc <- cbind(df.plot.acc, sub.plot)
    }
}

df.plot.acc$`Mean Accuracy (High-depth reference)` <- rowMeans(df.plot.acc)
# df.plot.acc$`Mean Accuracy` <- rowMeans(df.plot.acc)


# Balanced accuracy
df.plot.macrof1 <- data.frame(stringsAsFactors = F)
for (i in 1:length(sc.dataset)) {
    sc.data <- sc.dataset[i]
    sc.legend <- sc.legends[i]
    file.res <- paste0(path.res, 'results_', ref.data, '_', sc.data, '.txt')
    sub.res <- read.delim(file = file.res, row.names = 1, stringsAsFactors = F)
    sub.res <- sub.res[sub.res$term == 'Balanced accuracy',]
    sub.res$method <- gsub("singleR", "SingleR", sub.res$method)
    rownames(sub.res) <- sub.res$method
    sub.plot <- data.frame(sub.res[methods, 'value'], row.names = methods)
    names(sub.plot) <- sc.legend
    if (i == 1) {
        df.plot.macrof1 <- sub.plot
    } else {
        df.plot.macrof1 <- cbind(df.plot.macrof1, sub.plot)
    }
}

df.plot.macrof1$`Mean Balanced Accuracy (High-depth reference)` <- rowMeans(df.plot.macrof1)
# df.plot.macrof1$`Mean Balanced Accuracy` <- rowMeans(df.plot.macrof1)

mat.plot.Tasic <- cbind(df.plot.acc, df.plot.macrof1)
# mat.plot <- mat.plot[order(mat.plot$`Mean Accuracy`, decreasing = T),]
# 
# df.plot <- reshape2::melt(as.matrix(mat.plot))
# names(df.plot) <- c('Method', 'Dataset', 'Accuracy')
# vec.metrics <- c(rep('Tasic_Accuracy', 50), rep('Tasic_Balanced Accuracy', 50))
# df.plot$Type <- factor(vec.metrics, levels = c('Accuracy', 'Balanced Accuracy'))
# df.mean <- df.plot[df.plot$Dataset == 'Mean Accuracy', ]
# df.plot$Method <- factor(df.plot$Method, 
#                          levels = df.mean$Method[order(df.mean$Accuracy, decreasing = F)])
# sort.methods <- df.mean$Method[order(df.mean$Accuracy, decreasing = F)]
# color_type <- rep('0', nrow(df.plot))
# color_type[(df.plot$Dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy')) & 
#                (df.plot$Accuracy < 0.1)] <- '1'
# color_type[(df.plot$Dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy')) & 
#                (df.plot$Accuracy > 0.9)] <- '1'
# color_type[(df.plot$Dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy')) & 
#                (df.plot$Accuracy > 0.1) & (df.plot$Accuracy < 0.9)] <- '2'
# df.plot$Color_Type <- color_type

ref.data <- 'MCA'
sc.dataset <- c('Tasic2018', 'Campbell', 'HochgernerA', 'Mizrak')
sc.legends <- c('Tasic, mouse Neocortex',
                'Campbell, mouse HArc-ME',
                'Hochgerner, mouse DGH',
                'Mizrak, mouse V-SVZ')

for (i in 1:length(sc.dataset)) {
    sc.data <- sc.dataset[i]
    sc.legend <- sc.legends[i]
    file.res <- paste0(path.res, 'results_', ref.data, '_', sc.data, '.txt')
    sub.res <- read.delim(file = file.res, row.names = 1, stringsAsFactors = F)
    sub.res <- sub.res[sub.res$term == 'Accuracy',]
    sub.res$method <- gsub("singleR", "SingleR", sub.res$method)
    rownames(sub.res) <- sub.res$method
    sub.plot <- data.frame(sub.res[methods, 'value'], row.names = methods)
    names(sub.plot) <- sc.legend
    if (i == 1) {
        df.plot.acc <- sub.plot
    } else {
        df.plot.acc <- cbind(df.plot.acc, sub.plot)
    }
}

df.plot.acc$`Mean Accuracy (Low-depth reference)` <- rowMeans(df.plot.acc)
# df.plot.acc$`Mean Accuracy` <- rowMeans(df.plot.acc)


# Balanced accuracy
df.plot.macrof1 <- data.frame(stringsAsFactors = F)
for (i in 1:length(sc.dataset)) {
    sc.data <- sc.dataset[i]
    sc.legend <- sc.legends[i]
    file.res <- paste0(path.res, 'results_', ref.data, '_', sc.data, '.txt')
    sub.res <- read.delim(file = file.res, row.names = 1, stringsAsFactors = F)
    sub.res <- sub.res[sub.res$term == 'Balanced accuracy',]
    sub.res$method <- gsub("singleR", "SingleR", sub.res$method)
    rownames(sub.res) <- sub.res$method
    sub.plot <- data.frame(sub.res[methods, 'value'], row.names = methods)
    names(sub.plot) <- sc.legend
    if (i == 1) {
        df.plot.macrof1 <- sub.plot
    } else {
        df.plot.macrof1 <- cbind(df.plot.macrof1, sub.plot)
    }
}

df.plot.macrof1$`Mean Balanced Accuracy (Low-depth reference)` <- rowMeans(df.plot.macrof1)
# df.plot.macrof1$`Mean Balanced Accuracy` <- rowMeans(df.plot.macrof1)

mat.plot.MCA <- cbind(df.plot.acc, df.plot.macrof1)
mat.plot <- cbind(mat.plot.Tasic, mat.plot.MCA)
mat.plot$`Mean Accuracy` <- rowMeans(mat.plot[, c(5, 15)])
mat.plot$`Mean Balanced Accuracy` <- rowMeans(mat.plot[, c(10, 20)])
# mat.plot <- cbind(mat.plot[, 1:5], mat.plot[, 6:10], 
#                   mat.plot[, c(11:15)], mat.plot[, c(16:20)], mat.plot[, c(21:22)])
mat.plot <- mat.plot[order(mat.plot$`Mean Accuracy`, decreasing = T),]

df.plot <- reshape2::melt(as.matrix(mat.plot))
names(df.plot) <- c('Method', 'Dataset', 'Accuracy')
vec.metrics <- c(rep('Tasic_Accuracy', 50),  rep('Tasic_Balanced Accuracy', 50),  
                 rep('MCA_Accuracy', 50), rep('MCA_Balanced Accuracy', 50),
                 rep('Mean_Accuracy', 20))
df.plot$Type <- factor(vec.metrics, 
                       levels = c('Tasic_Accuracy', 'Tasic_Balanced Accuracy', 
                                  'MCA_Accuracy', 'MCA_Balanced Accuracy', 
                                  'Mean_Accuracy'))
df.mean <- df.plot[df.plot$Dataset == 'Mean Accuracy', ]
df.plot$Method <- factor(df.plot$Method, 
                         levels = rev(df.mean$Method[order(df.mean$Accuracy, decreasing = F)]))
sort.methods <- df.mean$Method[order(df.mean$Accuracy, decreasing = F)]
color_type <- rep('0', nrow(df.plot))
color_type[(df.plot$Dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy',
                                   'Mean Accuracy (High-depth reference)',
                                   'Mean Balanced Accuracy (High-depth reference)',
                                   'Mean Accuracy (Low-depth reference)',
                                   'Mean Balanced Accuracy (Low-depth reference)')) & 
               (df.plot$Accuracy < 0.1)] <- '1'
color_type[(df.plot$Dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy',
                                   'Mean Accuracy (High-depth reference)',
                                   'Mean Balanced Accuracy (High-depth reference)',
                                   'Mean Accuracy (Low-depth reference)',
                                   'Mean Balanced Accuracy (Low-depth reference)')) & 
               (df.plot$Accuracy > 0.85)] <- '1'
color_type[(df.plot$Dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy',
                                   'Mean Accuracy (High-depth reference)',
                                   'Mean Balanced Accuracy (High-depth reference)',
                                   'Mean Accuracy (Low-depth reference)',
                                   'Mean Balanced Accuracy (Low-depth reference)')) & 
               (df.plot$Accuracy > 0.1) & (df.plot$Accuracy < 0.85)] <- '2'
df.plot$Color_Type <- color_type

df.plot$Dataset[df.plot$Dataset %in% c('Mean Accuracy (High-depth reference)',
                                       'Mean Accuracy (Low-depth reference)')] <- 'Mean Accuracy'
df.plot$Dataset[df.plot$Dataset %in% c('Mean Balanced Accuracy (High-depth reference)',
                                       'Mean Balanced Accuracy (Low-depth reference)')] <- 'Mean Balanced Accuracy'
df.plot$Dataset <- factor(df.plot$Dataset, 
                          levels = c('Mizrak, mouse V-SVZ',
                                     'Hochgerner, mouse DGH',
                                     'Campbell, mouse HArc-ME',
                                     'Tasic, mouse Neocortex', 
                                     'Mean Balanced Accuracy', 'Mean Accuracy'))

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

ggsave(filename = 'heatmap_Tasic_MCA.png', 
       path = path.fig, plot = plot.heatmap,
       units = 'cm', height = 16, width = 20)

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
          axis.text.y = element_text(size = 9, color = 'black', family = 'Arial'),
          axis.text.x = element_text(size = 14, color = 'black', family = 'Arial',
                                     angle = 315, vjust = 1, hjust = 0),
          axis.title = element_text(size = 10, family = 'Arial'),
          legend.title = element_text(
              size = 12, color = "black", family = 'Arial'),
          legend.text = element_text(size = 11, family = 'Arial'),
          # legend.key.size = unit(1, 'cm'),
          legend.key = element_blank()) +
    guides(color = guide_legend(ncol = 1))
ggsave(filename = 'boxplot_Tasic_MCA.png', 
       path = path.fig, plot = plot.box,
       units = 'cm', height = 7, width = 20)
