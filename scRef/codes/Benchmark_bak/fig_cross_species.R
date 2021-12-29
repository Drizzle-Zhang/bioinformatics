library(ggplot2)
library(RColorBrewer)

# data
ref.dataset <- c('BaronM', 'BaronM', 'MCA', 'MCA')
sc.dataset <- c('panc8_celseq2', 'panc8_smartseq2', 'pbmc3k', 'pbmcsca_10Xv2')
sc.legends <- c('Baron(inDrop) -> Muraro(CEL-seq2), Pancreas',
                'Baron(inDrop) -> Segerstolpe(SMART-seq2), Pancreas',
                'Han(Microwell-seq) -> Butler(10X v1), PBMC',
                'Han(Microwell-seq) -> Ding(10X v2), PBMC')

path.res <- '/home/zy/scRef/Benchmark/cross_species/'
path.fig <- '/home/zy/scRef/figure/cross_species/'

methods <- c("scMAGIC", "singleR", "scmap-cell", "scmap-cluster", "CHETAH", "scPred",
             "sciBet", "singleCellNet", "scID", "scClassify")

# Accuracy
df.plot.acc <- data.frame(stringsAsFactors = F)
for (i in 1:length(sc.dataset)) {
    ref.data <- ref.dataset[i]
    sc.data <- sc.dataset[i]
    sc.legend <- sc.legends[i]
    file.res <- paste0(path.res, 'results_', ref.data, '_', sc.data, '.txt')
    sub.res <- read.delim(file = file.res, row.names = 1, stringsAsFactors = F)
    sub.res <- sub.res[sub.res$term == 'Accuracy',]
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


# Macro F1
df.plot.macrof1 <- data.frame(stringsAsFactors = F)
for (i in 1:length(sc.dataset)) {
    ref.data <- ref.dataset[i]
    sc.data <- sc.dataset[i]
    sc.legend <- sc.legends[i]
    file.res <- paste0(path.res, 'results_', ref.data, '_', sc.data, '.txt')
    sub.res <- read.delim(file = file.res, row.names = 1, stringsAsFactors = F)
    sub.res <- sub.res[sub.res$term == 'Macro F1',]
    rownames(sub.res) <- sub.res$method
    sub.plot <- data.frame(sub.res[methods, 'value'], row.names = methods)
    names(sub.plot) <- sc.legend
    if (i == 1) {
        df.plot.macrof1 <- sub.plot
    } else {
        df.plot.macrof1 <- cbind(df.plot.macrof1, sub.plot)
    }
}

df.plot.macrof1$`Mean Micro F1` <- rowMeans(df.plot.macrof1)

df.plot <- cbind(df.plot.acc, df.plot.macrof1)
df.plot <- df.plot[order(df.plot$`Mean Accuracy`, decreasing = T),]

pheatmap::pheatmap(df.plot,
                   color = colorRampPalette(rev(brewer.pal(n = 9, name =  "RdYlBu"))[1:9])(100),
                   cluster_rows = F, cluster_cols = F, scale = "none",
                   display_numbers = T, number_format = "%.2f", fontsize_number = 15, number_color = 'black',
                   show_rownames = T, show_colnames = T, fontsize_row = 18, fontsize_col = 15,
                   legend = T, drop_levels = F,
                   gaps_col = c(5, 10), angle_col = '315',
                   filename = paste0(path.fig, 'heatmap_CS.png'), width = 11, height = 12
)

# boxplot
df.box <- data.frame(stringsAsFactors = F)
mat.box1 <- df.plot[, c(1:4)]
for (method in rownames(mat.box1)) {
    for (dataset in colnames(mat.box1)) {
        df.sub <- data.frame(Accuracy = mat.box1[method, dataset], 
                             Method = method, Dataset = dataset, Type = 'Accuracy')
        df.box <- rbind(df.box, df.sub)
    }
}
mat.box2 <- df.plot[, c(6:9)]
for (method in rownames(mat.box2)) {
    for (dataset in colnames(mat.box2)) {
        df.sub <- data.frame(Accuracy = mat.box2[method, dataset], 
                             Method = method, Dataset = dataset, Type = 'Macro F1')
        df.box <- rbind(df.box, df.sub)
    }
}
df.box$Method <- factor(df.box$Method, levels = rev(rownames(df.plot)))

plot.box <- 
    ggplot(df.box, aes(x = Method, y = Accuracy, color = Type)) +
    geom_boxplot() +
    theme_classic() + coord_flip() + 
    # facet_wrap(~ Evaluator, scales = 'free_x', ncol = 2) +
    labs(title = "", y = '', x = '', fill = 'Metrics') + 
    theme(axis.text.y = element_blank())
# theme(axis.text = element_text(size = 9),
#       panel.grid = element_blank(),
#       panel.grid.major.y = element_line(color = 'grey', size = 0.2),
#       axis.text.x = element_text(angle = 45, hjust = 1),
#       axis.title = element_text(size = 12))
ggsave(filename = 'boxplot_CS.png', 
       path = path.fig, plot = plot.box,
       units = 'cm', height = 11, width = 14)
