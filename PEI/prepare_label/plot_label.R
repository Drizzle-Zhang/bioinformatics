setwd('C:\\Users\\Rain_Zhang\\Desktop\\PEI\\test_data')
library(ggplot2)

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




########################################################
df.num.pairs <- data.frame(
    CellLine = factor(c('GM12878', 'K562', 'MCF-7', 'HCT116', 'HeLa-S3', 'IMR-90')),
    NumPairs = c(93792, 195780, 199072, 177, 15383, 2807))

plot.num.pairs <- ggplot(df.num.pairs, aes(x = CellLine, y = NumPairs)) + 
    geom_bar(stat = "identity") + coord_flip() + 
    theme(axis.ticks = element_blank())+
    theme(axis.title = element_text(size = 17)) +   
    theme(axis.text = element_text(size = 12)) +   
    theme(panel.background=element_blank())   ## 去掉背景颜色

ggsave(
    plot = plot.num.pairs, path = './', filename = "num_pairs.png",
    units = 'cm', width = 15, height = 8)
