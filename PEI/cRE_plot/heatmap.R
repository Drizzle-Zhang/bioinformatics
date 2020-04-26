setwd('C:\\Users\\Rain_Zhang\\Desktop\\PEI\\test_data')
library(ggplot2)

df.mtx <- read.delim('./corr.txt', sep = '\t', stringsAsFactors = F, 
                     row.names = 'X')
rows <- c('score_dhs_enhancer', 
          'score_h3k4me3_enhancer', 'pval_h3k4me3_enhancer',
          'score_h3k27ac_enhancer', 'pval_h3k27ac_enhancer',
          'score_dhs_promoter',
          'score_h3k4me3_promoter', 'pval_h3k4me3_promoter',
          'score_h3k27ac_promoter', 'pval_h3k27ac_promoter',
          'gene_expression', 
          'distance', 'score_dhs_insulator', 'score_ctcf_insulator',
          'pcHi-C_ng2019', '3DIV', 'Thurman')
cols <- c('score_dhs_enhancer', 
          'score_h3k4me3_enhancer', 'pval_h3k4me3_enhancer',
          'score_h3k27ac_enhancer', 'pval_h3k27ac_enhancer',
          'score_dhs_promoter',
          'score_h3k4me3_promoter', 'pval_h3k4me3_promoter',
          'score_h3k27ac_promoter', 'pval_h3k27ac_promoter',
          'gene_expression', 
          'distance', 'score_dhs_insulator', 'score_ctcf_insulator',
          'pcHi.C_ng2019', 'X3DIV', 'Thurman')
df.mtx <- df.mtx[rows, cols]
df.heatmap <- data.frame()
labels <- row.names(df.mtx)
for (i in 1:length(labels)) {
    for (j in 1:length(labels)) {
        df.heatmap <- rbind(
            df.heatmap, data.frame(x = labels[i], y = labels[j], 
                                   value = df.mtx[i, j]))
    }
}
plot.heatmap <- ggplot(data = df.heatmap, aes(x, y)) + 
    geom_tile(aes(fill = value)) + 
    scale_fill_continuous(low = "#FFFAFA", high = "#FF0000") + 
    theme_bw() +
    theme(
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90)
    ) + 
    geom_text(aes(label = round(value, 2)), family = "Arial", size = 2)

ggsave(
    plot = plot.heatmap, path = './', filename = "Heatmap_features.png",
    units = 'cm', width = 15, height = 10)


#################################################
df.mtx <- read.delim('./heatmap_label4.txt', sep = '\t', stringsAsFactors = F)

plot.heatmap <- ggplot(data = df.mtx, aes(label1, label2)) + 
    geom_tile(aes(fill = score)) + 
    scale_fill_continuous(low = "#FFFAFA", high = "#FF0000") + 
    theme_bw() +
    theme(
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90)
    ) + 
    geom_text(aes(label = round(score, 2)), size = 2)

ggsave(
    plot = plot.heatmap, path = './', filename = "Heatmap_label4.png",
    units = 'cm', width = 35, height = 25)
