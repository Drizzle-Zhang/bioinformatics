setwd('/home/drizzle_zhang/driver_mutation/cRE_plot/HiC')
library(ggplot2)

df.compare <- read.delim('./compare_results.Cortex500.txt', 
                         sep = '\t', stringsAsFactors = F)

# interaction level
df.interaction <- data.frame()
col.names <- c("total", "count")
for (col.name in col.names) {
    df.sub <- df.compare[, c('label', col.name)]
    names(df.sub) <- c('Source', 'Number')
    df.sub$Label <- rep(col.name, dim(df.compare)[1])
    df.interaction <- rbind(df.interaction, df.sub)
}
df.ratio <- df.compare[, c('label', 'ratio')]
names(df.ratio) <- c('Source', 'Proportion')
df.interaction <- merge(df.interaction, df.ratio)
df.interaction <- df.interaction[df.interaction$Number != 0,]

plot.interaction <- ggplot(df.interaction) + 
    geom_bar(aes(x = Source, y = Number, fill = Label), stat = "identity",
             position = 'dodge') +
    scale_fill_discrete(
        name = 'Legend',
        breaks = c('total', 'count'),
        labels = c(
            'pcHi-C interations', 'pcHi-C interations annotated by cREs')) +
    geom_line(aes(x = Source, y = Proportion * 200000, linetype = 'Proportion', 
                  group = 1), color = '#ca3e1c') + 
    geom_point(aes(x = Source, y = Proportion * 200000), size = 3, shape = 1) + 
    geom_text(aes(x = Source, y = Proportion * 200000), label = 
                  paste(round(df.interaction$Proportion, 4)*100, '%', sep = ''),
              vjust = -1) + 
    theme(legend.title = element_blank(), 
          panel.background = element_rect(fill = 'white', colour = 'white'))
ggsave(
    plot = plot.interaction, path = './', filename = "interation.png",
    units = 'cm', width = 25, height = 12)

# egene-cRE level
df.egene <- data.frame()
col.names <- c("num_overlap_egene", "overlap_num")
col.names2 <- c("num_egene", "num_egene_pairs")
col.names3 <- c(1959, 22781)
for (i in 1:2) {
    df.sub <- df.compare[, c('label', col.names[i])]
    names(df.sub) <- c('Source', 'Number')
    df.sub$Label <- rep(col.names[i], dim(df.compare)[1])
    df.sub$Precision <- df.sub$Number / df.compare[, col.names2[i]]
    df.sub$Recall <- df.sub$Number / col.names3[i]
    df.egene <- rbind(df.egene, df.sub)
}

plot.egene <- ggplot(df.egene) + 
    geom_bar(aes(x = Source, y = Number, fill = Label), stat = "identity",
             position = 'dodge') +
    scale_fill_discrete(
        name = 'Legend',
        breaks = c('num_overlap_egene', 'overlap_num'),
        labels = c(
            'Expressed genes in eGene-cRE pairs overlapped with eQTL', 
            'eGene-cRE pairs overlapped with eQTL')) +
    theme(legend.title = element_blank(), 
          panel.background = element_rect(fill = 'white', colour = 'white'))
ggsave(
    plot = plot.egene, path = './', filename = "egene-cRE_pair.png",
    units = 'cm', width = 30, height = 13)

plot.egene <- ggplot(df.egene) + 
    geom_bar(aes(x = Source, y = Precision, fill = Label), stat = "identity",
             position = 'dodge') +
    scale_fill_discrete(
        name = 'Legend',
        breaks = c('num_overlap_egene', 'overlap_num'),
        labels = c(
            'Expressed genes in eGene-cRE pairs overlapped with eQTL /\n Expressed genes in eGene-cRE pairs', 
            'eGene-cRE pairs overlapped with eQTL /\n eGene-cRE pairs')) +
    theme(legend.title = element_blank(), 
          panel.background = element_rect(fill = 'white', colour = 'white'))
ggsave(
    plot = plot.egene, path = './', filename = "Precision_egene-cRE_pair.png",
    units = 'cm', width = 30, height = 13)


# compare with each other
df.mtx <- read.delim('./overlap.mtx', sep = '\t', stringsAsFactors = F, 
                     row.names = 'X')
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
    geom_text(aes(label = round(value,2)), family = "Arial", size = 2)
ggsave(
    plot = plot.heatmap, path = './', filename = "Heatmap_Hi-C_overlap.png",
    units = 'cm', width = 15, height = 10)

