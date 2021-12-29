library(ggplot2)
library(RColorBrewer)

# data
path.fig <- '/home/zy/scRef/figure/Tasic/'
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
                'Campbell, mouse hArc-ME',
                'Hochgerner, mouse DGH',
                'Mizrak, mouse V-SVZ')

path.res <- '/home/zy/scRef/Benchmark/mouse_brain/'
path.fig <- '/home/zy/scRef/figure/mouse_brain/'

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

df.plot.acc$`Mean Accuracy` <- rowMeans(df.plot.acc)


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

df.plot.macrof1$`Mean Balanced Accuracy` <- rowMeans(df.plot.macrof1)

mat.plot <- cbind(df.plot.acc, df.plot.macrof1)
mat.plot <- mat.plot[order(mat.plot$`Mean Accuracy`, decreasing = T),]

df.plot <- reshape2::melt(as.matrix(mat.plot))
names(df.plot) <- c('Method', 'Dataset', 'Accuracy')
vec.metrics <- c(rep('Accuracy', 50), rep('Balanced Accuracy', 50))
df.plot$Type <- factor(vec.metrics, levels = c('Accuracy', 'Balanced Accuracy'))
df.mean <- df.plot[df.plot$Dataset == 'Mean Accuracy', ]
df.plot$Method <- factor(df.plot$Method, 
                         levels = df.mean$Method[order(df.mean$Accuracy, decreasing = F)])
color_type <- rep('0', nrow(df.plot))
color_type[(df.plot$Dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy')) & 
               (df.plot$Accuracy < 0.2)] <- '1'
color_type[(df.plot$Dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy')) & 
               (df.plot$Accuracy > 0.8)] <- '1'
color_type[(df.plot$Dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy')) & 
               (df.plot$Accuracy > 0.2) & (df.plot$Accuracy < 0.8)] <- '2'
df.plot$Color_Type <- color_type

plot.heatmap <- 
    ggplot(data = df.plot, aes(x = Dataset, y = Method)) + 
    geom_tile(aes(fill = Accuracy)) + 
    facet_grid( ~ Type, scales = 'free', space = 'free') +
    scale_fill_gradientn(
        colors = colorRampPalette(rev(brewer.pal(n = 9, name =  "RdYlBu"))[1:9])(100)) + 
    scale_y_discrete(position = 'right') + 
    labs(x = '', y = '', fill = '') +
    geom_text(aes(label = round(Accuracy, 3), color = Color_Type), 
              family = "Arial", size = 3.0) +
    scale_color_manual(breaks = c('0', '1', '2'), 
                       values = c('transparent', 'white', 'black')) +
    theme(panel.background = element_rect(color = 'transparent',
                                          fill = 'transparent'), 
          axis.text.x = element_text(size = 13, color = "black",
                                     angle = 315, vjust = 1, hjust = 0),
          axis.text.y = element_text(
              size = 13, color = "black", face = 'bold', vjust = 0.5, hjust = 0.5),
          strip.text.x = element_text(
              size = 12, color = "black", face = "bold"), 
          strip.background = element_rect(
              color = "transparent", fill = "transparent"),
          legend.title = element_text(
              size = 12, color = "black", face = "bold"),
          legend.position = 'left') + 
    guides(color = F)

ggsave(filename = 'heatmap_Tasic.png', 
       path = path.fig, plot = plot.heatmap,
       units = 'cm', height = 14, width = 16)


# boxplot
df.box <- df.plot

plot.box <- 
    ggplot(df.box, aes(x = Method, y = Accuracy, color = Type)) +
    geom_boxplot(width = 0.7, outlier.size = 1) +
    scale_color_manual(breaks = c('Accuracy', 'Balanced Accuracy'),
                       values = c("#1F78B4", "#FF7F00")) +
    coord_flip() + 
    labs(title = "", y = 'Accuracy / Balanced Accuracy', x = '', color = '', fill = '') + 
    theme(panel.background = element_rect(color = 'black',
                                          fill = 'transparent'),
          axis.text.y = element_blank(), 
          panel.grid.major.y = element_line(color = 'lightgray', 
                                            size = 0.15, linetype = 'dashed'),
          legend.position = 'bottom',
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 12),
          legend.key = element_blank())
ggsave(filename = 'boxplot_Tasic.png', 
       path = path.fig, plot = plot.box,
       units = 'cm', height = 11, width = 10)


# labeled
# Accuracy
df.plot.acc <- data.frame(stringsAsFactors = F)
for (i in 1:length(sc.dataset)) {
    sc.data <- sc.dataset[i]
    sc.legend <- sc.legends[i]
    file.res <- paste0(path.res, 'results_', ref.data, '_', sc.data, '.txt')
    sub.res <- read.delim(file = file.res, row.names = 1, stringsAsFactors = F)
    sub.res <- sub.res[sub.res$term == 'labeled-Accuracy',]
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

df.plot.acc$`Mean Accuracy` <- rowMeans(df.plot.acc)

# Balanced accuracy
df.plot.macrof1 <- data.frame(stringsAsFactors = F)
for (i in 1:length(sc.dataset)) {
    sc.data <- sc.dataset[i]
    sc.legend <- sc.legends[i]
    file.res <- paste0(path.res, 'results_', ref.data, '_', sc.data, '.txt')
    sub.res <- read.delim(file = file.res, row.names = 1, stringsAsFactors = F)
    sub.res <- sub.res[sub.res$term == 'labeled-Balanced accuracy',]
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

df.plot.macrof1$`Mean Balanced Accuracy` <- rowMeans(df.plot.macrof1)

mat.plot <- cbind(df.plot.acc, df.plot.macrof1)
mat.plot <- mat.plot[order(mat.plot$`Mean Accuracy`, decreasing = T),]

df.plot <- reshape2::melt(as.matrix(mat.plot))
names(df.plot) <- c('Method', 'Dataset', 'Accuracy')
vec.metrics <- c(rep('Accuracy', 50), rep('Balanced Accuracy', 50))
df.plot$Type <- factor(vec.metrics, levels = c('Accuracy', 'Balanced Accuracy'))
df.mean <- df.plot[df.plot$Dataset == 'Mean Accuracy', ]
df.plot$Method <- factor(df.plot$Method, 
                         levels = df.mean$Method[order(df.mean$Accuracy, decreasing = F)])
color_type <- rep('0', nrow(df.plot))
color_type[(df.plot$Dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy')) & 
               (df.plot$Accuracy < 0.2)] <- '1'
color_type[(df.plot$Dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy')) & 
               (df.plot$Accuracy > 0.8)] <- '1'
color_type[(df.plot$Dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy')) & 
               (df.plot$Accuracy > 0.2) & (df.plot$Accuracy < 0.8)] <- '2'
df.plot$Color_Type <- color_type

plot.heatmap <- 
    ggplot(data = df.plot, aes(x = Dataset, y = Method)) + 
    geom_tile(aes(fill = Accuracy)) + 
    facet_grid( ~ Type, scales = 'free', space = 'free') +
    scale_fill_gradientn(
        colors = colorRampPalette(rev(brewer.pal(n = 9, name =  "RdYlBu"))[1:9])(100)) + 
    scale_y_discrete(position = 'right') + 
    labs(x = '', y = '', fill = '') +
    geom_text(aes(label = round(Accuracy, 3), color = Color_Type), 
              family = "Arial", size = 3.0) +
    scale_color_manual(breaks = c('0', '1', '2'), 
                       values = c('transparent', 'white', 'black')) +
    theme(panel.background = element_rect(color = 'transparent',
                                          fill = 'transparent'), 
          axis.text.x = element_text(size = 13, color = "black",
                                     angle = 315, vjust = 1, hjust = 0),
          axis.text.y = element_text(
              size = 13, color = "black", face = 'bold', vjust = 0.5, hjust = 0.5),
          strip.text.x = element_text(
              size = 12, color = "black", face = "bold"), 
          strip.background = element_rect(
              color = "transparent", fill = "transparent"),
          legend.title = element_text(
              size = 12, color = "black", face = "bold"),
          legend.position = 'left') + 
    guides(color = F)

ggsave(filename = 'heatmap_Tasic.png', 
       path = path.fig, plot = plot.heatmap,
       units = 'cm', height = 14, width = 16)


# boxplot
df.box <- data.frame()
for (i in 1:length(sc.dataset)) {
    sc.data <- sc.dataset[i]
    sc.legend <- sc.legends[i]
    file.res <- paste0(path.res, 'results_', ref.data, '_', sc.data, '.txt')
    sub.res <- read.delim(file = file.res, row.names = 1, stringsAsFactors = F)
    sub.res <- sub.res[sub.res$term == 'Percent of unassigned',]
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


plot.box <- 
    ggplot(df.box, aes(x = Method, y = Accuracy, color = Type)) +
    geom_boxplot(width = 0.7, outlier.size = 1) +
    scale_color_manual(breaks = c('Accuracy', 'Balanced Accuracy'),
                       values = c("#1F78B4", "#FF7F00")) +
    coord_flip() + 
    labs(title = "", y = 'Accuracy / Balanced Accuracy', x = '', color = '', fill = '') + 
    theme(panel.background = element_rect(color = 'black',
                                          fill = 'transparent'),
          axis.text.y = element_blank(), 
          panel.grid.major.y = element_line(color = 'lightgray', 
                                            size = 0.15, linetype = 'dashed'),
          legend.position = 'bottom',
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 12),
          legend.key = element_blank())
ggsave(filename = 'boxplot_Tasic.png', 
       path = path.fig, plot = plot.box,
       units = 'cm', height = 11, width = 10)








