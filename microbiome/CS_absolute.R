library(tidyverse)

# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.out.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)

level = 'class'
################################### class
file = paste0('/home/drizzle_zhang/microbiome/result/3.Community_Structure/abundance/absolute_abundance/', 
              level, '.csv')
class <- read.table(file = file, header = TRUE, row.names = 1,
                    sep = ',', stringsAsFactors = FALSE,
                    check.names = FALSE)

# time series
path.plot <- paste0('/home/drizzle_zhang/microbiome/result/3.Community_Structure/barplot_',
                    level)
series.time <- unique(df.meta$Time)

for (sub.time in series.time) {
    # select meta
    sel.meta <- df.meta[df.meta$Time == sub.time,]
    sel.meta <- sel.meta[sel.meta$Dose %in% c(0, 3),]
    row.names(sel.meta) <- sel.meta$Sample
    
    # select sample
    use.sample <- sel.meta$Sample
    sub.class <- class[, use.sample]
    
    # 挑选 top10 门类群，并将 top10 外的类群合并为“Others”
    class_top10 <- sub.class[1:10, ]
    class_top10['Others', ] <- colSums(sub.class) - colSums(class_top10)
    class_top10.diff <- class_top10
    class_top10$Taxonomy <- fct_inorder(rownames(class_top10))
    
    # 创建分组信息数据集
    group <- data.frame(variable = sel.meta$Sample,
                        group = sel.meta$Group)
    
    class_top10 <- pivot_longer(data = class_top10, cols = -Taxonomy,
                                names_to = "variable", values_to = "value")
    # 按同类项进行合并
    class_top10 <- merge(class_top10, group, by = 'variable')
    
    # 绘制带分面的柱状图
    plot.species <- 
        ggplot(data = class_top10, aes(variable, 100 * value, fill = Taxonomy)) + geom_col(position = 'stack', width = 0.6) +
        # 利用facet_wrap 按组分面并排
        facet_wrap(~ group, scales = 'free_x', ncol = 2) +
        # 用 scale_fill_brewer 的默认11色配色板
        scale_fill_brewer(palette = "Set3") +
        labs(x = '', y = 'Absolute Abundance(%)') +
        # 去掉背景网格，但保留横网线利于比较各组数据，一目了然
        theme(panel.grid.minor.y = element_line(colour = "black"),
              panel.background = element_rect(color = 'black', 
                                              fill = 'transparent'))
    ggsave(plot = plot.species, path = path.plot, 
           filename = paste0('absolute_species', sub.time, '_0vs3.png'))
    
    # degree of difference
    # class_top10.diff <- class[1:10, ]
    diff.in <- cbind(t(class_top10.diff), sel.meta)
    types <- row.names(class_top10.diff)
    df.out <- data.frame()
    interval <- 0.02
    for (type in types) {
        sub.in <- diff.in[, c(type, 'Group')]
        fold.change <- log10(mean(sub.in[sub.in$Group == 'Treat', type]) / 
                                 mean(sub.in[sub.in$Group == 'Control', type]))
        out.wilcox <- 
            wilcox.test(as.formula(paste0(type, ' ~ Group')), data = sub.in)
        p.value <- round(out.wilcox$p.value, 4)
        if (fold.change < 0) {
            position <- fold.change - interval
        } else {
            position <- fold.change + interval
        }
        df.out <- rbind(df.out, 
                        data.frame(type = type, fold.change = fold.change, 
                                   p.value = p.value, position = position))
    }
    plot.fc <- 
        ggplot(aes(x = type, y = fold.change), data = df.out) + 
        geom_bar(stat = 'identity', width = 0.5) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
        geom_text(aes(label = paste0('p = ', p.value), y = position), 
                  size = 3) + 
        labs(x = '', y = 'log10(fold change)')
    
    ggsave(plot = plot.fc, path = path.plot, 
           filename = paste0('absolute_foldchange_', sub.time, '_0vs3.png'))
    
}
