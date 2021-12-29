library(tidyverse)
library(RColorBrewer)
library(stringr)
library(ggplot2)

# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.out.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)

# gender
gender = 'male'
df.meta <- df.meta[df.meta$Gender == gender, ]
# vec.dose <- c(0, 1, 2, 3)
vec.dose <- c(0, 3)

# level = 'family'
level = 'genus'
################################### class
file = paste0('/home/drizzle_zhang/microbiome/result/3.Community_Structure/barplot/Top30_', 
              level, '.csv')
class <- read.table(file = file, header = TRUE, row.names = 1,
                    sep = '\t', stringsAsFactors = FALSE,
                    check.names = FALSE)
row.names(class) <- str_replace_all(row.names(class), '-', '_')
    
# time series
path.plot <- paste0('/home/drizzle_zhang/microbiome/result/3.Community_Structure/barplot_',
                    level)
if (!file.exists(path.plot)) {
    dir.create(path.plot)
}
series.time <- unique(df.meta$Time)
series.time <- c(-1, 1,  5,  9, 17, 21, 25, 29, 33, 41, 49, 60, 68, 84)
mod <- 'sel'
df.plot.point <- data.frame()
for (sub.time in series.time) {
    # select meta
    sel.meta <- df.meta[df.meta$Time == sub.time,]
    sel.meta <- sel.meta[sel.meta$Dose %in% vec.dose,]
    row.names(sel.meta) <- sel.meta$Sample
    
    # select sample
    use.sample <- sel.meta$SampleName
    sub.class <- class[row.names(class) != 'Other', use.sample]
    
    # 挑选 top10 门类群，并将 top10 外的类群合并为“Others”
    class_top10 <- sub.class[1:20, ]
    names(class_top10) <- sel.meta$Sample
    class_top10['Others', ] <- 1 - colSums(class_top10)
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
    colourCount = length(unique(class_top10$Taxonomy))
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))
    plot.species <- 
        ggplot(data = class_top10, aes(variable, 100 * value, fill = Taxonomy)) + geom_col(position = 'stack', width = 0.6) +
        # 利用facet_wrap 按组分面并排
        # facet_wrap(~ group, scales = 'free_x', ncol = 2) +
        facet_grid(~ group, scales = 'free_x', space = 'free') +
        # 用 scale_fill_brewer 的默认11色配色板
        # scale_fill_brewer(palette = "Set3") +
        # more than 11 kinds of color
        scale_fill_manual(values = getPalette(colourCount)) + 
        labs(x = '', y = 'Relative Abundance(%)') +
        # 去掉背景网格，但保留横网线利于比较各组数据，一目了然
        theme(axis.text.x = element_text(size = 6, angle = 45, 
                                         vjust = 1, hjust = 1),
              legend.text = element_text(size = 5),
              legend.key.size = unit(0.5, 'cm'),
              panel.grid.minor.y = element_line(colour = "black"),
              panel.background = element_rect(color = 'black', 
                                              fill = 'transparent')) + 
        guides(fill = guide_legend(ncol = 2))
        
    ggsave(plot = plot.species, path = path.plot, 
           # filename = paste0('gender_diff_relative_species_',
           #                   paste0(as.character(vec.dose), collapse = ''), 
           #                   '_', sub.time, '.png'),
           filename = paste0(gender, '_relative_species_',
                             paste0(as.character(vec.dose), collapse = ''), 
                             '_', sub.time, '.png'),
           height = 10, width = 20, units = 'cm')
    
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
        # a = sub.in[sub.in$Group == 'Treat', type]
        # b = sub.in[sub.in$Group == 'Control',type]
        # out.wilcox <- ks.test(a, b, data = sub.in)
        p.value <- out.wilcox$p.value
        p.value.round <- round(out.wilcox$p.value, 4)
        if (fold.change < 0) {
            position <- fold.change - interval
        } else {
            position <- fold.change + interval
        }
        df.out <- rbind(df.out, 
                        data.frame(type = type, fold.change = fold.change, 
                                   p.value = p.value, position = position,
                                   time = sub.time, p.value.round = p.value.round))
    }
    plot.fc <- 
        ggplot(aes(x = type, y = fold.change), data = df.out) + 
        geom_bar(stat = 'identity', width = 0.5) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
        geom_text(aes(label = paste0('p = ', p.value.round), y = position), 
                  size = 2) + 
        labs(x = '', y = 'log10(fold change)')
    
    ggsave(plot = plot.fc, path = path.plot, 
           filename = paste0(gender, '_relative_foldchange_',
                             paste0(as.character(vec.dose), collapse = ''), 
                             '_', sub.time, '.png'),
           height = 10, width = 20, units = 'cm')
    
    df.plot.point <- rbind(df.plot.point, df.out)
    
}

df.plot.point$log.pval <- -log10(df.plot.point$p.value)
# df.plot.point[df.plot.point$log.pval == Inf, 'log.pval'] <- 4
file.bubble <- paste0(path.plot, '/', gender, '_Bubble_', mod, '_',
                      paste0(as.character(vec.dose), collapse = ''), '.txt')
write.table(df.plot.point, sep = '\t', file = file.bubble, quote = F, row.names = F)
plot.bubble <- 
    ggplot(data = df.plot.point, 
           aes(x = time, y = type, size = log.pval, color = fold.change)) + 
    geom_point(fill = 'cornsilk') + 
    scale_colour_gradient2(low = 'blue', mid = "white", high = 'red') + 
    labs(x = 'Time', y = 'Family', color = 'Fold Change', size = 'P Value') + 
    theme(panel.background = element_rect(color = 'gray',
                                          fill = 'transparent'))
ggsave(plot = plot.bubble, path = path.plot, 
       filename = paste0(gender, '_Bubble_', mod, '_',
                         paste0(as.character(vec.dose), collapse = ''), '.png'),
       height = 15, width = 20, units = 'cm')

