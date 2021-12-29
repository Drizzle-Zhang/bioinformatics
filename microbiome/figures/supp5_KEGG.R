library(ggplot2)
library(scales)
library(patchwork)

# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.out.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)

# KEGG
file.KEGG.L3 <- '/home/drizzle_zhang/microbiome/result/9.PICRUSt/origin_data/KO_KEGG_L3.txt'
df.db.KEGG <- read.delim(file.KEGG.L3, row.names = 1, header = F, 
                         stringsAsFactors = F)
names.KEGG.L3 <- as.character(df.db.KEGG[1,])
names(df.db.KEGG) <- names.KEGG.L3
df.db.KEGG <- df.db.KEGG[-1,]

# cutoff
type.cutoff <- 'fdr'

# dose
# vec.dose <- c(0, 1, 2, 3)
vec.dose <- c(0, 1)

# time series
# series.time <- unique(df.meta$Time)
# mod <- 'all'
series.time <- c(-1, 1,  5,  9, 17, 21, 25, 29, 33, 41, 49, 60, 68, 84)
mod <- 'sel'
# series.time <- c(-1, 1,  5,  9, 17, 21, 25, 29, 33, 37, 41, 45, 49)
# mod <- 'old'
############################# GSEA
level <- 'L3'
# male
gender <- 'male'
df.meta.gender <- df.meta[df.meta$Gender == gender, ]
path.plot.male <- paste0('/home/drizzle_zhang/microbiome/result/9.PICRUSt/heatmap_',
                         gender, '_', type.cutoff)
if (!file.exists(path.plot.male)) {
    dir.create(path.plot.male)
}
path.out <- '/home/drizzle_zhang/microbiome/result/Figs/'

df.GSEA <- data.frame(ID = names.KEGG.L3)
for (sub.time in series.time) {
    file.GSEA <- paste0(path.plot.male, "/GSEA_", level, "_", sub.time, "_",
                        paste0(as.character(vec.dose), collapse = ''),
                        ".txt")
    # file.GSEA <- paste0(path.plot.male, "/GSEA_", level, "_", sub.time, 
    #                     paste0(as.character(vec.dose), collapse = ''), 
    #                     ".txt")
    sub.GSEA <- read.delim(file.GSEA, row.names = 1)
    sub.GSEA$logPval <- -log10(sub.GSEA$pvalue) * 
        (sub.GSEA$enrichmentScore / abs(sub.GSEA$enrichmentScore))
    # sub.GSEA$order <- rank(sub.GSEA$logPval)
    sub.GSEA <- sub.GSEA[, c("Description", "logPval")]
    names(sub.GSEA) <- c("ID", sub.time)
    df.GSEA <- merge(df.GSEA, sub.GSEA, by = 'ID', all = T)
}
row.names(df.GSEA) <- df.GSEA$ID
df.GSEA$ID <- NULL
df.GSEA[is.na(df.GSEA)] <- 0
df.GSEA.score <- df.GSEA
df.GSEA <- as.data.frame(apply(df.GSEA.score, 2, rank))

# sort
df.sort <- data.frame(stringsAsFactors = F)
for (row in row.names(df.GSEA)) {
    for (col in as.numeric(names(df.GSEA))) {
        if ((col > 0) & (col < 42)) {
            df.sort <- rbind(df.sort, data.frame(pathway = row, time = col,
                                                 value = df.GSEA[row, as.character(col)],
                                                 stringsAsFactors = F))
        }
    }
}
df.sort$ID <- paste(df.sort$pathway, df.sort$time, sep = '_')
df.sort <- df.sort[order(df.sort$value),]
sort.value <- df.sort$value
df.ks <- data.frame(stringsAsFactors = F)
for (pathway in names.KEGG.L3) {
    sub.sort <- df.sort[df.sort$pathway == pathway, 'value']
    enrich.control <- ks.test(sub.sort, sort.value, alternative = 'greater')
    enrich.treat <- ks.test(sub.sort, sort.value, alternative = 'less')
    df.ks <- rbind(df.ks, data.frame(pathway = pathway, 
                                     pvalue.control = enrich.control$p.value,
                                     pvalue.treat = enrich.treat$p.value))
}
df.ks$qvalue.control <- p.adjust(df.ks$pvalue.control, method = 'fdr')
df.ks$qvalue.treat <- p.adjust(df.ks$pvalue.treat, method = 'fdr')

# use ks score to plot
df.ks.male <- df.ks
# df.ks.male.filter <- df.ks.male[
#     df.ks.male$qvalue.control < 0.05 | df.ks.male$qvalue.treat < 0.05,]
df.ks.male.filter <- df.ks.male[
    df.ks.male$pvalue.control < 0.012 | df.ks.male$pvalue.treat < 0.012,]
# df.ks.male.filter <- df.ks.male[
#     df.ks.male$pvalue.control < 0.1 | df.ks.male$pvalue.treat < 0.03,]
# df.ks.male.filter <- df.ks.male[
#     df.ks.male$pvalue.control < 0.01 | df.ks.male$pvalue.treat < 0.01,]
df.ks.male.filter <- 
    df.ks.male.filter[!(df.ks.male.filter$pathway %in% 
                            c('ABC transporters', 'Lipopolysaccharide biosynthesis')),]
log10Pval <- c()
for (i in row.names(df.ks.male.filter)) {
    pvalue.control <- -log10(df.ks.male.filter[i, 'pvalue.control'])
    pvalue.treat <- -log10(df.ks.male.filter[i, 'pvalue.treat'])
    if (pvalue.control > pvalue.treat) {
        log10Pval <- c(log10Pval, -pvalue.control)
    } else {
        log10Pval <- c(log10Pval, pvalue.treat)
    }
}
df.ks.male.filter$log10Pval <- log10Pval
df.ks.male.filter <- df.ks.male.filter[
    order(df.ks.male.filter$log10Pval, decreasing = T), ]
vec.color <- c()
for (pval in df.ks.male.filter$log10Pval) {
    if (pval < 0) {
        vec.color <- c(vec.color, 'Enrich in Control')
    } else {
        vec.color <- c(vec.color, 'Enrich in Treatment')
    }
}
df.ks.male.filter$color <- factor(vec.color, 
                                  levels = c('Enrich in Treatment', 'Enrich in Control'))
plot.male <- 
    ggplot(data = df.ks.male.filter, aes(x = reorder(pathway, X = log10Pval), 
                                         y = log10Pval, fill = color)) + 
    geom_bar(stat = 'identity') + 
    labs(x = 'Pathway', y = expression(paste("-log"[10], "(adj", italic("P"), "-value)")), 
         fill = '') + 
    scale_fill_manual(values = c(muted("red"), muted("blue"))) + 
    coord_flip() + 
    theme_bw() +
    theme(panel.background = element_rect(color = 'black', size = 1.5,
                                          fill = 'transparent'),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(colour = "gray", size = 0.1,
                                            linetype = 2),
          axis.text.x = element_text(size = 9, color = "black", family = 'Arial'),
          axis.text.y = element_text(size = 10, color = "black", family = 'Arial'))
# ggsave(filename = paste0("/Male_Combine_Sum_GSEA_", mod, '_',
#                          paste0(as.character(vec.dose), collapse = ''), ".png"),
#        path = path.out, plot = plot.male,
#        height = 15, width = 20, units = 'cm')

# heatmap
df.GSEA <- df.GSEA.score
df.GSEA.male <- df.GSEA[as.character(df.ks.male.filter$pathway),]
df.heatmap.male <- data.frame(stringsAsFactors = F)
for (pathway in row.names(df.GSEA.male)) {
    for (sub.time in names(df.GSEA.male)) {
        df.heatmap.male <- 
            rbind(df.heatmap.male, 
                  data.frame(pathway = pathway, time = sub.time,
                             score = df.GSEA.male[pathway, sub.time],
                             stringsAsFactors = F))
    }
}
df.heatmap.male$pathway <- 
    factor(df.heatmap.male$pathway, 
           levels = rev(as.character(df.ks.male.filter$pathway)), ordered = T)
df.heatmap.male$time <- 
    factor(df.heatmap.male$time, 
           levels = as.character(series.time), ordered = T)
plot.heatmap.male <- 
    ggplot(data = df.heatmap.male, 
           aes(x = time, y = pathway, fill = score)) + 
    geom_tile() + 
    scale_fill_gradient2(low = muted("blue"), high = muted("red"), mid = "#F5F5F5") + 
    labs(x = 'Time', y = 'Pathway', fill = 'Enrichment Score') + 
    theme_bw() + 
    theme(
        panel.border = element_blank(),
        panel.background = element_rect(color = 'transparent', size = 0,
                                        fill = 'transparent'),
        panel.grid  = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 9, color = "black", family = 'Arial'), 
        legend.text = element_text(size = 9)
    )
# ggsave(filename = paste0("/Male_Combine_Heatmap_GSEA_",  mod, '_',
#                          paste0(as.character(vec.dose), collapse = ''), ".png"),
#        path = path.out, plot = plot.heatmap.male,
#        height = 15, width = 20, units = 'cm')

plot.final.male <- plot.male + plot.heatmap.male + plot_layout(widths = c(1, 1.6),
                                                      guides = 'collect')

ggsave(plot = plot.final.male, path = path.out, 
       filename = paste0("/Male_GSEA_",  
                         paste0(as.character(vec.dose), collapse = ''), ".png"),
       height = 11, width = 25, units = 'cm')





# female
gender <- 'female'
df.meta.gender <- df.meta[df.meta$Gender == gender, ]
path.plot.female <- paste0('/home/drizzle_zhang/microbiome/result/9.PICRUSt/heatmap_',
                           gender, '_', type.cutoff)
df.GSEA <- data.frame(ID = names.KEGG.L3)
for (sub.time in series.time) {
    file.GSEA <- paste0(path.plot.female, "/GSEA_", level, "_", sub.time,  "_",
                        paste0(as.character(vec.dose), collapse = ''), 
                        ".txt")
    sub.GSEA <- read.delim(file.GSEA, row.names = 1)
    sub.GSEA$logPval <- -log10(sub.GSEA$pvalue) * 
        (sub.GSEA$enrichmentScore / abs(sub.GSEA$enrichmentScore))
    sub.GSEA <- sub.GSEA[, c("Description", "logPval")]
    names(sub.GSEA) <- c("ID", sub.time)
    df.GSEA <- merge(df.GSEA, sub.GSEA, by = 'ID', all = T)
}
row.names(df.GSEA) <- df.GSEA$ID
df.GSEA$ID <- NULL
df.GSEA[is.na(df.GSEA)] <- 0

# sort
df.sort <- data.frame(stringsAsFactors = F)
for (row in row.names(df.GSEA)) {
    for (col in as.numeric(names(df.GSEA))) {
        if ((col == 1) | (col == 25) | ((col >= 33) & (col < 85))) {
        # if ((col == 1) | ((col >= 1) & (col < 85))) {
            df.sort <- rbind(df.sort, data.frame(pathway = row, time = col,
                                                 value = df.GSEA[row, as.character(col)],
                                                 stringsAsFactors = F))
        }
    }
}
df.sort$ID <- paste(df.sort$pathway, df.sort$time, sep = '_')
df.sort <- df.sort[order(df.sort$value),]
sort.value <- df.sort$value
df.ks <- data.frame(stringsAsFactors = F)
for (pathway in names.KEGG.L3) {
    sub.sort <- df.sort[df.sort$pathway == pathway, 'value']
    enrich.control <- ks.test(sub.sort, sort.value, alternative = 'greater')
    enrich.treat <- ks.test(sub.sort, sort.value, alternative = 'less')
    df.ks <- rbind(df.ks, data.frame(pathway = pathway, 
                                     pvalue.control = enrich.control$p.value,
                                     pvalue.treat = enrich.treat$p.value))
}
df.ks$qvalue.control <- p.adjust(df.ks$pvalue.control, method = 'fdr')
df.ks$qvalue.treat <- p.adjust(df.ks$pvalue.treat, method = 'fdr')

# use ks score to plot
df.ks.female <- df.ks
# df.ks.female.filter <- df.ks.female[
#     df.ks.female$qvalue.control < 0.15 | df.ks.female$qvalue.treat < 0.4,]
# df.ks.female.filter <- df.ks.female[
#     df.ks.female$pvalue.control < 0.035 | df.ks.female$pvalue.treat < 0.043,]
# df.ks.female.filter <- df.ks.female[
#     df.ks.female$pvalue.control < 0.05 | df.ks.female$pvalue.treat < 0.07,]
# df.ks.female.filter <- df.ks.female[
#     df.ks.female$pvalue.control < 0.01 | df.ks.female$pvalue.treat < 0.003,]
df.ks.female.filter <- df.ks.female[
    df.ks.female$pvalue.control < 0.2 | df.ks.female$pvalue.treat < 0.1,]
df.ks.female.filter <- df.ks.female.filter[
    df.ks.female.filter$pathway != 'ABC transporters',]
sel.pathway <- c('Lipopolysaccharide biosynthesis proteins',
                 'Fatty acid biosynthesis', 'Nitrogen metabolism',
                 'Peptidases', 'Phosphotransferase system (PTS)',
                 'Benzoate degradation',
                 'Tryptophan metabolism', 'Starch and sucrose metabolism',
                 'Peroxisome', 'Galactose metabolism', 'Lysine biosynthesis',
                 'Terpenoid backbone biosynthesis')
df.ks.female.filter <- df.ks.female.filter[df.ks.female.filter$pathway %in% sel.pathway,]
log10Pval <- c()
for (i in row.names(df.ks.female.filter)) {
    pvalue.control <- -log10(df.ks.female.filter[i, 'pvalue.control'])
    pvalue.treat <- -log10(df.ks.female.filter[i, 'pvalue.treat'])
    if (pvalue.control > pvalue.treat) {
        log10Pval <- c(log10Pval, -pvalue.control)
    } else {
        log10Pval <- c(log10Pval, pvalue.treat)
    }
}
df.ks.female.filter$log10Pval <- log10Pval
df.ks.female.filter <- df.ks.female.filter[
    order(df.ks.female.filter$log10Pval, decreasing = T), ]
vec.color <- c()
for (pval in df.ks.female.filter$log10Pval) {
    if (pval > 0) {
        vec.color <- c(vec.color, 'Enrich in Treatment')
    } else {
        vec.color <- c(vec.color, 'Enrich in Control')
    }
}
df.ks.female.filter$color <- factor(vec.color, 
                                    levels = c('Enrich in Treatment', 'Enrich in Control'))
plot.female <- 
    ggplot(data = df.ks.female.filter, aes(x = reorder(pathway, X = log10Pval), 
                                           y = log10Pval, fill = color)) + 
    geom_bar(stat = 'identity') + 
    labs(x = 'Pathway', y = expression(paste("-log"[10], "(adj", italic("P"), "-value)")), 
         fill = '') + 
    scale_fill_manual(values = c(muted("red"), muted("blue"))) + 
    coord_flip() + 
    theme_bw() +
    theme(panel.background = element_rect(color = 'black', size = 1.5,
                                          fill = 'transparent'),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(colour = "gray", size = 0.1,
                                            linetype = 2),
          axis.text.x = element_text(size = 9, color = "black", family = 'Arial'),
          axis.text.y = element_text(size = 10, color = "black", family = 'Arial'))
# ggsave(filename = paste0("/Female_Combine_Sum_GSEA_",  mod, '_',
#                          paste0(as.character(vec.dose), collapse = ''), ".png"),
#        path = path.out, plot = plot.female,
#        height = 15, width = 20, units = 'cm')


# heatmap
df.GSEA.female <- df.GSEA[as.character(df.ks.female.filter$pathway),]
df.heatmap.female <- data.frame(stringsAsFactors = F)
for (pathway in row.names(df.GSEA.female)) {
    for (sub.time in names(df.GSEA.female)) {
        df.heatmap.female <- 
            rbind(df.heatmap.female, 
                  data.frame(pathway = pathway, time = sub.time,
                             score = df.GSEA.female[pathway, sub.time],
                             stringsAsFactors = F))
    }
}
df.heatmap.female$pathway <- 
    factor(df.heatmap.female$pathway, 
           levels = rev(as.character(df.ks.female.filter$pathway)), ordered = T)
df.heatmap.female$time <- 
    factor(df.heatmap.female$time, 
           levels = as.character(series.time), ordered = T)
plot.heatmap.female <- 
    ggplot(data = df.heatmap.female, 
           aes(x = time, y = pathway, fill = score)) + 
    geom_tile() + 
    scale_fill_gradient2(low = muted("blue"), high = muted("red"), mid = "#F5F5F5") + 
    labs(x = 'Time', y = 'Pathway', fill = 'Enrichment Score') + 
    theme_bw() +
    theme(
        panel.border = element_blank(),
        panel.background = element_rect(color = 'transparent', size = 1.5,
                                        fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 9, color = "black", family = 'Arial'), 
        legend.text = element_text(size = 9))
# ggsave(filename = paste0("/Female_Combine_Heatmap_GSEA_",  mod, '_',
#                          paste0(as.character(vec.dose), collapse = ''), ".png"),
#        path = path.out, plot = plot.heatmap.female,
#        height = 15, width = 20, units = 'cm')

plot.final.female <- plot.female + plot.heatmap.female + 
    plot_layout(widths = c(1, 1.6), guides = 'collect')

ggsave(plot = plot.final.female, path = path.out, 
       filename = paste0("/Female_GSEA_",  
                         paste0(as.character(vec.dose), collapse = ''), ".png"),
       height = 8, width = 25, units = 'cm')
