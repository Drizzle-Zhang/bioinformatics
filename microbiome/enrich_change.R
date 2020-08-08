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
