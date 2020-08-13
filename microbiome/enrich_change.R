library(ggplot2)

# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.out.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)

# KEGG
file.KEGG.L3 <- '/home/drizzle_zhang/microbiome/result/9.PICRUSt/origin_data/KO_KEGG_L3.txt'
df.db.KEGG <- read.delim(file.KEGG.L3, row.names = 1)

# cutoff
type.cutoff <- 'fdr'

# dose
# vec.dose <- c(0, 1, 2, 3)
vec.dose <- c(0, 3)

# time series
series.time <- unique(df.meta$Time)

############################# GSEA
level <- 'L3'
# male
gender <- 'male'
df.meta.gender <- df.meta[df.meta$Gender == gender, ]
path.plot <- paste0('/home/drizzle_zhang/microbiome/result/9.PICRUSt/heatmap_',
                    gender, '_', type.cutoff)

for (sub.time in series.time) {
    file.GSEA <- paste0(path.plot, "/GSEA_", level, "_", sub.time, 
                        paste0(as.character(vec.dose), collapse = ''), 
                        ".txt")
    sub.GSEA <- read.delim(file.GSEA, row.names = 1)
    sub.GSEA$logPval <- log10(sub.GSEA$pvalue) * 
        (sub.GSEA$enrichmentScore / abs(sub.GSEA$enrichmentScore))
    sub.GSEA <- sub.GSEA[, c("Description", "logPval")]
    names(sub.GSEA) <- c("ID", sub.time)
    if (sub.time == 'A') {
        df.GSEA <- sub.GSEA
    } else {
        df.GSEA <- merge(df.GSEA, sub.GSEA, by = 'ID', all = T)
    }
}
row.names(df.GSEA) <- df.GSEA$ID
df.GSEA$ID <- NULL
df.GSEA[is.na(df.GSEA)] <- 0
df.plot.GSEA <- df.GSEA[order(abs(rowSums(df.GSEA[,names(df.GSEA) != 'A'])), decreasing = T)[c(1:7, 9:10)],]
df.plot <- data.frame()
for (row in row.names(df.plot.GSEA)) {
    for (col in 1:length(names(df.plot.GSEA))) {
        df.plot <- rbind(df.plot, data.frame(pathway = row, time = col,
                                             value = df.plot.GSEA[row, col]))
    }
}
plot.male <- ggplot(data = df.plot, aes(x = time, y = value, color = pathway)) + 
    geom_point() + 
    geom_smooth(method = lm, formula = y ~ poly(x, 4), se = F) +
    labs(x = 'Time', y = '-log10(Pvalue)', color = 'Pathway')
ggsave(filename = paste0("/Combine_GSEA_", 
                         paste0(as.character(vec.dose), collapse = ''), ".png"),
       path = path.plot, plot = plot.male,
       height = 12, width = 20, units = 'cm')
# sum of scores
df.sum.male <- rowSums(df.plot.GSEA[,names(df.plot.GSEA) != 'A'])
df.sum.male <- data.frame(pathway = names(df.sum.male), logPval = df.sum.male)
vec.color <- c()
for (pval in df.sum.male$logPval) {
    if (pval > 0) {
        vec.color <- c(vec.color, 'Enrich in Control')
    } else {
        vec.color <- c(vec.color, 'Enrich in Treatment')
    }
}
df.sum.male$color <- vec.color
df.sum.male <- df.sum.male[order(df.sum.male$logPval, decreasing = T), ]
plot.sum.male <- ggplot(data = df.sum.male, aes(x = reorder(pathway, X = logPval), 
                                                y = logPval, fill = color)) + 
    geom_bar(stat = 'identity') + 
    labs(x = 'Pathway', y = '-log10(Pvalue)', fill = '') + 
    coord_flip() + 
    theme(panel.background = element_rect(color = 'gray', size = 2,
                                          fill = 'transparent'),
          panel.grid.major.y = element_line(colour = "gray", size = 0.2,
                                            linetype = 2))
ggsave(filename = paste0("/Combine_Sum_GSEA_", 
                         paste0(as.character(vec.dose), collapse = ''), ".png"),
       path = path.plot, plot = plot.sum.male,
       height = 10, width = 20, units = 'cm')

# female
gender <- 'female'
df.meta.gender <- df.meta[df.meta$Gender == gender, ]
path.plot <- paste0('/home/drizzle_zhang/microbiome/result/9.PICRUSt/heatmap_',
                    gender, '_', type.cutoff)

for (sub.time in series.time) {
    file.GSEA <- paste0(path.plot, "/GSEA_", level, "_", sub.time, 
                        paste0(as.character(vec.dose), collapse = ''), 
                        ".txt")
    sub.GSEA <- read.delim(file.GSEA, row.names = 1)
    sub.GSEA$logPval <- log10(sub.GSEA$pvalue) * 
        (sub.GSEA$enrichmentScore / abs(sub.GSEA$enrichmentScore))
    sub.GSEA <- sub.GSEA[, c("Description", "logPval")]
    names(sub.GSEA) <- c("ID", sub.time)
    if (sub.time == 'A') {
        df.GSEA <- sub.GSEA
    } else {
        df.GSEA <- merge(df.GSEA, sub.GSEA, by = 'ID', all = T)
    }
}
row.names(df.GSEA) <- df.GSEA$ID
df.GSEA$ID <- NULL
df.GSEA[is.na(df.GSEA)] <- 0
df.plot.GSEA <- df.GSEA[order(abs(rowSums(df.GSEA[,names(df.GSEA) != 'A'])), decreasing = T)[1:9],]
df.plot.GSEA <- df.plot.GSEA[
    !(row.names(df.plot.GSEA) %in% c("Two-component system", "Protein kinases")),]
df.plot <- data.frame()
for (row in row.names(df.plot.GSEA)) {
    for (col in 1:length(names(df.plot.GSEA))) {
        df.plot <- rbind(df.plot, data.frame(pathway = row, time = col,
                                             value = df.plot.GSEA[row, col]))
    }
}
plot.female <- ggplot(data = df.plot, aes(x = time, y = value, color = pathway)) + 
    geom_point() + 
    geom_smooth(method = lm, formula = y ~ poly(x, 4), se = F) +
    labs(x = 'Time', y = '-log10(Pvalue)', color = 'Pathway')
ggsave(filename = paste0("/Combine_GSEA_", 
                         paste0(as.character(vec.dose), collapse = ''), ".png"),
       path = path.plot, plot = plot.female,
       height = 12, width = 20, units = 'cm')
# sum of scores
df.sum.female <- rowSums(df.plot.GSEA[,names(df.plot.GSEA) != 'A'])
df.sum.female <- data.frame(pathway = names(df.sum.female), logPval = df.sum.female)
vec.color <- c()
for (pval in df.sum.female$logPval) {
    if (pval > 0) {
        vec.color <- c(vec.color, 'Enrich in Control')
    } else {
        vec.color <- c(vec.color, 'Enrich in Treatment')
    }
}
df.sum.female$color <- vec.color
df.sum.female <- df.sum.female[order(df.sum.female$logPval, decreasing = T), ]
plot.sum.female <- ggplot(data = df.sum.female, aes(x = reorder(pathway, X = logPval), 
                                                y = logPval, fill = color)) + 
    geom_bar(stat = 'identity') + 
    labs(x = 'Pathway', y = '-log10(Pvalue)', fill = '') + 
    coord_flip() + 
    theme(panel.background = element_rect(color = 'gray', size = 2,
                                          fill = 'transparent'),
          panel.grid.major.y = element_line(colour = "gray", size = 0.2,
                                            linetype = 2))
ggsave(filename = paste0("/Combine_Sum_GSEA_", 
                         paste0(as.character(vec.dose), collapse = ''), ".png"),
       path = path.plot, plot = plot.sum.female,
       height = 10, width = 20, units = 'cm')


################################ KEGG enrichment analysis
# gender
gender <- 'male'
df.meta.gender <- df.meta[df.meta$Gender == gender, ]
path.plot <- paste0('/home/drizzle_zhang/microbiome/result/9.PICRUSt/heatmap_',
                    gender, '_', type.cutoff)

# pathway-pvalue matrix
for (sub.time in series.time) {
    file.control <- paste0(path.plot, "/Enrich_control_", sub.time, 
                           paste0(as.character(vec.dose), collapse = ''), 
                           ".txt")
    if (file.exists(file.control)) {
        df.enrich <- read.delim(file.control, row.names = 1)
    } else {
        df.enrich <- data.frame(pathway = names(df.db.KEGG),
                                pvalue = rep(1, length(names(df.db.KEGG))))
    }
    sub.enrich <- data.frame(df.enrich[, c('pathway', 'pvalue')])
    names(sub.enrich) <- c('pathway', sub.time)
    if (sub.time == 'A') {
        mat.enrich <- sub.enrich
    } else {
        mat.enrich <- merge(mat.enrich, sub.enrich, by = 'pathway')
    }
}
row.names(mat.enrich) <- mat.enrich$pathway
mat.enrich$pathway <- NULL

sel.pathway <- 
    c('ABC.transporters', 'Transporters', 'Phosphotransferase.system..PTS.',
      'Galactose.metabolism', 'Lipid.metabolism', 'Starch.and.sucrose.metabolism',
      'Terpenoid.backbone.biosynthesis', 'Synthesis.and.degradation.of.ketone.bodies',
      'Peptidoglycan.biosynthesis')
sel.mat.enrich <- mat.enrich[sel.pathway,]
sel.mat.enrich <- -log10(sel.mat.enrich)
df.plot <- data.frame()
for (row in row.names(sel.mat.enrich)) {
    for (col in 1:length(names(sel.mat.enrich))) {
        df.plot <- rbind(df.plot, data.frame(pathway = row, time = col,
                                             value = sel.mat.enrich[row, col]))
    }
}
plot.male <- ggplot(data = df.plot, aes(x = time, y = value, color = pathway)) + 
    geom_point() + 
    geom_smooth(method = lm, formula = y ~ poly(x, 4), se = F) +
    labs(x = 'Time', y = '-log10(Pvalue)', color = 'Pathway')
ggsave(filename = paste0("/Combine_pathway_control_enrich_", 
                         paste0(as.character(vec.dose), collapse = ''), ".png"),
       path = path.plot, plot = plot.male,
       height = 12, width = 20, units = 'cm')

# ggplot(data = df.plot, aes(x = time, y = value, color = pathway)) + 
#     geom_point() + 
#     geom_line()

# gender
gender <- 'female'
df.meta.gender <- df.meta[df.meta$Gender == gender, ]
path.plot <- paste0('/home/drizzle_zhang/microbiome/result/9.PICRUSt/heatmap_',
                    gender, '_', type.cutoff)

# pathway-pvalue matrix
for (sub.time in series.time) {
    file.control <- paste0(path.plot, "/Enrich_control_", sub.time, 
                           paste0(as.character(vec.dose), collapse = ''), 
                           ".txt")
    if (file.exists(file.control)) {
        df.enrich <- read.delim(file.control, row.names = 1)
    } else {
        df.enrich <- data.frame(pathway = names(df.db.KEGG),
                                pvalue = rep(1, length(names(df.db.KEGG))))
    }
    sub.enrich <- data.frame(df.enrich[, c('pathway', 'pvalue')])
    names(sub.enrich) <- c('pathway', sub.time)
    if (sub.time == 'A') {
        mat.enrich <- sub.enrich
    } else {
        mat.enrich <- merge(mat.enrich, sub.enrich, by = 'pathway')
    }
}
row.names(mat.enrich) <- mat.enrich$pathway
mat.enrich$pathway <- NULL

sel.pathway <- 
    c('Transporters', 'Phosphotransferase.system..PTS.',
      'Galactose.metabolism')
sel.mat.enrich <- mat.enrich[sel.pathway,]
sel.mat.enrich <- -log10(sel.mat.enrich)
df.plot <- data.frame()
for (row in row.names(sel.mat.enrich)) {
    for (col in 1:length(names(sel.mat.enrich))) {
        df.plot <- rbind(df.plot, data.frame(pathway = row, time = col,
                                             value = sel.mat.enrich[row, col]))
    }
}
plot.male <- ggplot(data = df.plot, aes(x = time, y = value, color = pathway)) + 
    geom_point() + 
    geom_smooth(method = lm, formula = y ~ poly(x, 4), se = F) +
    labs(x = 'Time', y = '-log10(Pvalue)', color = 'Pathway')
ggsave(filename = paste0("/Combine_pathway_control_enrich_", 
                         paste0(as.character(vec.dose), collapse = ''), ".png"),
       path = path.plot, plot = plot.male,
       height = 12, width = 20, units = 'cm')
