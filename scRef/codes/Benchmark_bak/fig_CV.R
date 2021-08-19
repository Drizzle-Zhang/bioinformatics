library(ggplot2)
require("RColorBrewer")
path.fig <- '/home/zy/scRef/figure/equal_celltype/'

# CV
datasets <- c('Campbell', 'panc8_indrop', 'pbmcsca_10Xv2', 'Tasic')
# GSE_ids <- c('GSE93374', 'GSE84133', 'GSE132044', 'GSE71585')
GSE_ids <- c('Campbell, Mouse hypothalamic Arc-ME', 
             'Baron, Human pancreas', 
             'Ding, PBMC', 'Tasic, Primary visual cortex')
path.output <- '/home/zy/scRef/cross_validation/train4_test1/'

methods <- c("scRef", "SingleR", "scmap-cell", "scmap-cluster", "CHETAH", "scPred",
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

df.plot.CV$`Mean Accuracy (Cross validation)` <- rowMeans(df.plot.CV)


# cross platform
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
    names(sub.plot) <- paste0(ref.GSE, '(', ref.plat, ')', ' -> ', sc.GSE, '(', sc.plat, ')', ', ', tissue)
    if (i == 1) {
        df.plot.CP <- sub.plot
    } else {
        df.plot.CP <- cbind(df.plot.CP, sub.plot)
    }
}

df.plot.CP$`Mean Accuracy (Cross platform)` <- rowMeans(df.plot.CP)

df.plot <- cbind(df.plot.CV, df.plot.CP)
df.plot$`Mean Accuracy` <- (df.plot$`Mean Accuracy (Cross validation)` + 
                                df.plot$`Mean Accuracy (Cross platform)`)/2
rownames(df.plot) <- gsub("scRef", "scMAGIC", rownames(df.plot))
df.plot <- df.plot[order(df.plot$`Mean Accuracy`, decreasing = T),]

anno_col <- data.frame(Type = c(rep('Cross validation', 5), 
                                rep('Cross platform', 5), 'Average'), 
                       row.names = colnames(df.plot))
col_colors <- list(Type = c(`Cross validation` = 'dimgray',
                            `Cross platform` = 'gainsboro',
                            `Average` = 'white'))
# col_colors <- list(Type = c('dimgray', 'gainsboro', 'white'))

pheatmap::pheatmap(df.plot,
                   color = colorRampPalette(rev(brewer.pal(n = 9, name =  "RdYlBu"))[1:9])(100),
                   border_color = NA, 
                   cluster_rows = F, cluster_cols = F, scale = "none",
                   display_numbers = T, number_format = "%.3f", fontsize_number = 10, number_color = 'black',
                   show_rownames = T, show_colnames = T, legend = F, fontsize_row = 14, fontsize_col = 10,
                   annotation_col = anno_col, 
                   annotation_legend = F, annotation_names_col = F,
                   annotation_colors = col_colors,
                   gaps_col = c(5, 10), angle_col = '315',
                   filename = paste0(path.fig, 'heatmap_accuracy.png'), width = 7, height = 9
)

# boxplot
df.box <- data.frame(stringsAsFactors = F)
mat.box1 <- df.plot[, c(1:4)]
for (method in rownames(mat.box1)) {
    for (dataset in colnames(mat.box1)) {
        df.sub <- data.frame(Accuracy = mat.box1[method, dataset], 
                             Method = method, Dataset = dataset, Type = 'Cross validation')
        df.box <- rbind(df.box, df.sub)
    }
}
mat.box2 <- df.plot[, c(6:9)]
for (method in rownames(mat.box2)) {
    for (dataset in colnames(mat.box2)) {
        df.sub <- data.frame(Accuracy = mat.box2[method, dataset], 
                             Method = method, Dataset = dataset, Type = 'Cross platform')
        df.box <- rbind(df.box, df.sub)
    }
}
df.box$Method <- factor(df.box$Method, levels = rev(rownames(df.plot)))
df.box$Type <- factor(df.box$Type, levels = rev(c('Cross validation', 'Cross platform')))

plot.box <- 
    ggplot(df.box,
                   aes(x = Method, y = Accuracy, color = Method, fill = Type)) +
    geom_boxplot(width = 0.5, outlier.size = 1) +
    scale_fill_manual(breaks = c('Cross validation', 'Cross platform'),
                      values = c('dimgray', 'gainsboro')) + 
    scale_color_manual(breaks = rev(rownames(df.plot)),
                       values = brewer.pal(10,"Paired")) + 
    theme_classic() + coord_flip() + 
    # facet_wrap(~ Evaluator, scales = 'free_x', ncol = 2) +
    labs(title = "", y = 'Accuracy', x = '', color = '', fill = '') + 
    theme(axis.text.y = element_blank(), 
          panel.grid.major.y = element_line(color = 'lightgray', 
                                            size = 0.15, linetype = 'dashed'),
          legend.position = 'bottom',
          axis.text.x = element_text(size = 7),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 7),
          legend.key.size = unit(0.5, 'cm')) + 
    guides(color = guide_legend(
              ncol = 3, 
              byrow = TRUE,
              reverse = T),
           fill = guide_legend(
               ncol = 1, 
               byrow = TRUE,
               reverse = T))
    # theme(axis.text = element_text(size = 9),
    #       panel.grid = element_blank(),
    #       panel.grid.major.y = element_line(color = 'grey', size = 0.2),
    #       axis.text.x = element_text(angle = 45, hjust = 1),
    #       axis.title = element_text(size = 12))
ggsave(filename = 'boxplot_Accuracy.png', 
       path = path.fig, plot = plot.box,
       units = 'cm', height = 14, width = 11)









