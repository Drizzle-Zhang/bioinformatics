library(tidyverse)
library(RColorBrewer)
library(stringr)
library(ggplot2)

# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)
df.meta$Gender[df.meta$Gender == 'female'] <- 'Female'
df.meta$Gender[df.meta$Gender == 'male'] <- 'Male'
df.meta$Group[df.meta$Group == 'Treat'] <- 'Treatment'

level = 'family'
file = paste0('/home/drizzle_zhang/microbiome/result/3.Community_Structure/barplot/Top30_', 
              level, '.csv')
class <- read.table(file = file, header = TRUE, row.names = 1,
                    sep = '\t', stringsAsFactors = FALSE,
                    check.names = FALSE)
row.names(class) <- str_replace_all(row.names(class), '-', '_')
series.time <- unique(df.meta$Time)
vec.group <- c('Control', 'Treatment')
vec.gender <- c('Male', 'Female')
path.plot <- '/home/drizzle_zhang/microbiome/result/Figs/'

list.vec.dose <- list(c(0, 1, 2, 3), c(0, 1), c(0, 2), c(0, 3))
for (vec.dose in list.vec.dose) {
    # average 
    df.average <- data.frame()
    df.group <- data.frame(stringsAsFactors = F)
    for (sub.time in series.time) {
        # select meta
        sel.meta <- df.meta[df.meta$Time == sub.time,]
        sel.meta <- sel.meta[sel.meta$Dose %in% vec.dose,]
        row.names(sel.meta) <- sel.meta$Sample
        
        # select sample
        use.sample <- sel.meta$SampleName
        sub.class <- class[row.names(class) != 'Other', use.sample]
        class_top10 <- sub.class[1:15, ]
        names(class_top10) <- use.sample
        class_top10['Others', ] <- 1 - colSums(class_top10)
        for (gender in vec.gender) {
            sub.average <- data.frame()
            for (group in vec.group) {
                sub.sample <- 
                    intersect(sel.meta$SampleName[
                        (sel.meta$Group == group) & (sel.meta$Gender == gender)],
                                        colnames(class_top10))
                sub.mtx <- class_top10[,sub.sample]
                sub.average <- rbind(sub.average, rowMeans(sub.mtx))
            }
            rownames(sub.average) <- paste(sub.time, vec.group, gender, sep = '_')
            colnames(sub.average) <- rownames(class_top10)
            df.average <- rbind(df.average, sub.average)
            df.group <- rbind(df.group, 
                              data.frame(ColName = rownames(sub.average),
                                         Time = rep(sub.time, nrow(sub.average)),
                                         Gender = rep(gender, nrow(sub.average)),
                                         Group = vec.group, 
                                         RowName = paste(sub.time, vec.group, sep = '_')))
        }
    }
    df.plot <- as.data.frame(t(df.average))
    df.plot$Taxonomy <- fct_inorder(colnames(df.average))
    df.plot <- pivot_longer(data = df.plot, cols = -Taxonomy,
                                names_to = "ColName", values_to = "value")
    # 按同类项进行合并
    df.plot <- merge(df.plot, df.group, by = 'ColName')
    # 绘制带分面的柱状图
    colourCount = length(unique(df.plot$Taxonomy))
    # brewer.pal(10, "Paired")
    color_vec <- c("#6A3D9A", "#CAB2D6", "#33A02C", "#B2DF8A", 
                   "#1F78B4", "#A6CEE3", "#E31A1C", "#FB9A99",
                   "#FF7F00", "#FDBF6F")
    getPalette = colorRampPalette(color_vec)
    # set.seed(123)
    # vec.color <- sample(getPalette(colourCount), colourCount)
    vec.color <- c(getPalette(colourCount - 1), 'gray')
    plot.species <- 
        ggplot(data = df.plot, aes(Group, 100 * value, fill = Taxonomy)) + 
        geom_col(position = 'stack', width = 0.6) +
        # 利用facet_wrap 按组分面并排
        # facet_wrap(~ group, scales = 'free_x', ncol = 2) +
            facet_grid(Gender ~ Time, scales = 'free', space = 'free') +
        # 用 scale_fill_brewer 的默认11色配色板
        # scale_fill_brewer(palette = "Set3") +
        # more than 11 kinds of color
        scale_fill_manual(values = vec.color) + 
        labs(x = '', y = 'Relative Abundance(%)') +
        # 去掉背景网格，但保留横网线利于比较各组数据，一目了然
        theme(axis.text.x = element_text(size = 8, angle = 45, color = 'black', 
                                         vjust = 1, hjust = 1),
              axis.text.y = element_text(size = 8, color = 'black'),
              legend.text = element_text(size = 8),
              legend.key.size = unit(0.5, 'cm'),
              strip.background = element_rect(
                  color = 'black', fill = 'transparent'),
              strip.text.y = element_text(angle = 270, color = 'black', size = 10),
              # panel.grid.major.y = element_line(colour = "gray"),
              panel.background = element_rect(color = 'transparent', 
                                              fill = 'transparent')) + 
        guides(fill = guide_legend(ncol = 1))
    ggsave(plot = plot.species, path = path.plot, 
           filename = paste0('CS_', level, '_',
                             paste0(as.character(vec.dose), collapse = ''), 
                             '.png'),
           height = 15, width = 25, units = 'cm')
}
