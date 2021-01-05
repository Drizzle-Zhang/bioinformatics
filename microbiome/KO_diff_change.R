library(ggplot2)
library(scales)

# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.out.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)
# KO-KEGG
file.meta.KEGG <- '/home/drizzle_zhang/microbiome/result/9.PICRUSt/origin_data/metadata_KEGG.txt'
df.meta.KEGG <- read.delim(file.meta.KEGG, row.names = 1)
df.meta.KEGG$KO_id <- row.names(df.meta.KEGG)

# cutoff
type.cutoff <- 'fdr'

# dose
# vec.dose <- c(0, 1, 2, 3)
vec.dose <- c(0, 3)

# time series
series.time <- unique(df.meta$Time)

########################## KO diff
# male
gender <- 'male'
df.meta.gender <- df.meta[df.meta$Gender == gender, ]

# path
path.plot <- paste0('/home/drizzle_zhang/microbiome/result/9.PICRUSt/heatmap_',
                    gender, '_', type.cutoff)

for (sub.time in series.time) {
    file.res <- paste0(path.plot, "/OUT_edgeR_", sub.time, 
                       paste0(as.character(vec.dose), collapse = ''), 
                       ".txt")
    res.edgeR <- read.delim(file.res, row.names = 1)
    res.edgeR$KO_id <- row.names(res.edgeR)
    logPval <- -log10(res.edgeR$PValue) * 
        (res.edgeR$logFC / abs(res.edgeR$logFC))
    logPval[is.na(logPval)] <- 0
    df.pval <- data.frame(KO_id = res.edgeR$KO_id, logPval = logPval)
    names(df.pval) <- c('KO_id', sub.time)
    if (sub.time == 'A') {
        mat.pval <- df.pval
    } else {
        mat.pval <- merge(mat.pval, df.pval, by = 'KO_id')
    }
    
}
row.names(mat.pval) <- mat.pval$KO_id
mat.pval$KO_id <- NULL

# sort
df.sort <- data.frame(stringsAsFactors = F)
for (row in row.names(mat.pval)) {
    for (col in 1:length(names(mat.pval))) {
        if (!col %in% c(1)) {
            df.sort <- rbind(df.sort, data.frame(KO_id = row, time = col,
                                             value = mat.pval[row, col],
                                             stringsAsFactors = F))
        }
    }
}
df.sort$ID <- paste(df.sort$KO_id, df.sort$time, sep = '_')
df.sort <- df.sort[order(df.sort$value),]
sort.value <- df.sort$value
df.ks <- data.frame(stringsAsFactors = F)
for (KO_id in row.names(mat.pval)) {
    sub.sort <- df.sort[df.sort$KO_id == KO_id, 'value']
    enrich.control <- ks.test(sub.sort, sort.value, alternative = 'less')
    enrich.treat <- ks.test(sub.sort, sort.value, alternative = 'greater')
    df.ks <- rbind(df.ks, data.frame(KO_id = KO_id, 
                                     pvalue.control = enrich.control$p.value,
                                     pvalue.treat = enrich.treat$p.value))
}
df.ks$qvalue.control <- p.adjust(df.ks$pvalue.control, method = 'BH')
df.ks$qvalue.treat <- p.adjust(df.ks$pvalue.treat, method = 'BH')

df.ks.male <- df.ks
df.ks.male.filter <- df.ks.male[
    df.ks.male$qvalue.control < 0.05 | df.ks.male$qvalue.treat < 0.05,]
log10Pval <- c()
for (i in row.names(df.ks.male.filter)) {
    pvalue.control <- -log10(df.ks.male.filter[i, 'pvalue.control'])
    pvalue.treat <- -log10(df.ks.male.filter[i, 'pvalue.treat'])
    if (pvalue.control > pvalue.treat) {
        log10Pval <- c(log10Pval, pvalue.control)
    } else {
        log10Pval <- c(log10Pval, -pvalue.treat)
    }
}
df.ks.male.filter$log10Pval <- log10Pval
df.ks.male.filter <- df.ks.male.filter[
    order(df.ks.male.filter$log10Pval, decreasing = T), ]

# heatmap
df.mat.male <- mat.pval[as.character(df.ks.male.filter$KO_id),]
df.heatmap.male <- data.frame(stringsAsFactors = F)
for (KO_id in row.names(df.mat.male)) {
    for (sub.time in names(df.mat.male)) {
        df.heatmap.male <- 
            rbind(df.heatmap.male, 
                  data.frame(KO_id = KO_id, time = sub.time,
                             score = df.mat.male[KO_id, sub.time],
                             stringsAsFactors = F))
    }
}
df.heatmap.male$KO_id <- 
    factor(df.heatmap.male$KO_id, 
           levels = as.character(df.ks.male.filter$KO_id), ordered = T)
plot.heatmap.male <- 
    ggplot(data = df.heatmap.male, 
           aes(x = time, y = KO_id, fill = score)) + 
    geom_tile() + 
    scale_fill_gradient2(low = muted("blue"), high = muted("red")) + 
    labs(x = '', y = '', fill = '-log10(Pvalue)') + 
    theme(panel.background = element_rect(color = 'white', size = 1.5,
                                          fill = 'transparent'),
          axis.ticks = element_blank(),
          axis.text.y = element_text(size = 6), 
          legend.text = element_text(size = 10))
ggsave(filename = paste0("/Combine_Heatmap_KO_", 
                         paste0(as.character(vec.dose), collapse = ''), ".png"),
       path = path.plot, plot = plot.heatmap.male,
       height = 35, width = 15, units = 'cm')

# record diff KO
df.record.male <- merge(df.ks.male.filter, df.meta.KEGG, by = 'KO_id')
df.record.male <- 
    df.record.male[order(df.record.male$log10Pval, decreasing = T), 
                   c('KO_id', 'log10Pval', 
                     "metadata_KEGG_Description", "metadata_KEGG_Pathways")]
file.record.male <- paste0(path.plot, '/Combine_records_KO_', 
                      paste0(as.character(vec.dose), collapse = ''), ".txt")
write.table(df.record.male, file = file.record.male, sep = '\t', quote = F,
            row.names = F, col.names = T)


# female
gender <- 'female'
df.meta.gender <- df.meta[df.meta$Gender == gender, ]

# path
path.plot <- paste0('/home/drizzle_zhang/microbiome/result/9.PICRUSt/heatmap_',
                    gender, '_', type.cutoff)

for (sub.time in series.time) {
    file.res <- paste0(path.plot, "/OUT_edgeR_", sub.time, 
                       paste0(as.character(vec.dose), collapse = ''), 
                       ".txt")
    res.edgeR <- read.delim(file.res, row.names = 1)
    res.edgeR$KO_id <- row.names(res.edgeR)
    logPval <- -log10(res.edgeR$PValue) * 
        (res.edgeR$logFC / abs(res.edgeR$logFC))
    logPval[is.na(logPval)] <- 0
    df.pval <- data.frame(KO_id = res.edgeR$KO_id, logPval = logPval)
    names(df.pval) <- c('KO_id', sub.time)
    if (sub.time == 'A') {
        mat.pval <- df.pval
    } else {
        mat.pval <- merge(mat.pval, df.pval, by = 'KO_id')
    }
    
}
row.names(mat.pval) <- mat.pval$KO_id
mat.pval$KO_id <- NULL

# sort
df.sort <- data.frame(stringsAsFactors = F)
for (row in row.names(mat.pval)) {
    for (col in 1:length(names(mat.pval))) {
        if (!col %in% c(1)) {
            df.sort <- rbind(df.sort, data.frame(KO_id = row, time = col,
                                                 value = mat.pval[row, col],
                                                 stringsAsFactors = F))
        }
    }
}
df.sort$ID <- paste(df.sort$KO_id, df.sort$time, sep = '_')
df.sort <- df.sort[order(df.sort$value),]
sort.value <- df.sort$value
df.ks <- data.frame(stringsAsFactors = F)
for (KO_id in row.names(mat.pval)) {
    sub.sort <- df.sort[df.sort$KO_id == KO_id, 'value']
    enrich.control <- ks.test(sub.sort, sort.value, alternative = 'less')
    enrich.treat <- ks.test(sub.sort, sort.value, alternative = 'greater')
    df.ks <- rbind(df.ks, data.frame(KO_id = KO_id, 
                                     pvalue.control = enrich.control$p.value,
                                     pvalue.treat = enrich.treat$p.value))
}
df.ks$qvalue.control <- p.adjust(df.ks$pvalue.control, method = 'BH')
df.ks$qvalue.treat <- p.adjust(df.ks$pvalue.treat, method = 'BH')

df.ks.female <- df.ks
df.ks.female.filter <- df.ks.female[
    df.ks.female$qvalue.control < 0.05 | df.ks.female$qvalue.treat < 0.05,]
log10Pval <- c()
for (i in row.names(df.ks.female.filter)) {
    pvalue.control <- -log10(df.ks.female.filter[i, 'pvalue.control'])
    pvalue.treat <- -log10(df.ks.female.filter[i, 'pvalue.treat'])
    if (pvalue.control > pvalue.treat) {
        log10Pval <- c(log10Pval, pvalue.control)
    } else {
        log10Pval <- c(log10Pval, -pvalue.treat)
    }
}
df.ks.female.filter$log10Pval <- log10Pval
df.ks.female.filter <- df.ks.female.filter[
    order(df.ks.female.filter$log10Pval, decreasing = T), ]

# heatmap
df.mat.female <- mat.pval[as.character(df.ks.female.filter$KO_id),]
df.heatmap.female <- data.frame(stringsAsFactors = F)
for (KO_id in row.names(df.mat.female)) {
    for (sub.time in names(df.mat.female)) {
        df.heatmap.female <- 
            rbind(df.heatmap.female, 
                  data.frame(KO_id = KO_id, time = sub.time,
                             score = df.mat.female[KO_id, sub.time],
                             stringsAsFactors = F))
    }
}
df.heatmap.female$KO_id <- 
    factor(df.heatmap.female$KO_id, 
           levels = as.character(df.ks.female.filter$KO_id), ordered = T)
plot.heatmap.female <- 
    ggplot(data = df.heatmap.female, 
           aes(x = time, y = KO_id, fill = score)) + 
    geom_tile() + 
    scale_fill_gradient2(low = muted("blue"), high = muted("red")) + 
    labs(x = '', y = '', fill = '-log10(Pvalue)') + 
    theme(panel.background = element_rect(color = 'white', size = 1.5,
                                          fill = 'transparent'),
          axis.ticks = element_blank(),
          axis.text.y = element_text(size = 6), 
          legend.text = element_text(size = 10))
ggsave(filename = paste0("/Combine_Heatmap_KO_", 
                         paste0(as.character(vec.dose), collapse = ''), ".png"),
       path = path.plot, plot = plot.heatmap.female,
       height = 35, width = 15, units = 'cm')

# record diff KO
df.record.female <- merge(df.ks.female.filter, df.meta.KEGG, by = 'KO_id')
df.record.female <- 
    df.record.female[order(df.record.female$log10Pval, decreasing = T), 
                   c('KO_id', 'log10Pval', 
                     "metadata_KEGG_Description", "metadata_KEGG_Pathways")]
file.record.female <- paste0(path.plot, '/Combine_records_KO_', 
                           paste0(as.character(vec.dose), collapse = ''), ".txt")
write.table(df.record.female, file = file.record.female, sep = '\t', quote = F,
            row.names = F, col.names = T)




