library(ggplot2)

# KEGG
file.KEGG <- '/home/drizzle_zhang/microbiome/result/9.PICRUSt/ko_predict/ko_predictions.txt'
df.KEGG <- read.delim(file.KEGG, row.names = 1)
df.KEGG <- df.KEGG[ ,paste0('S', 1:(dim(df.KEGG)[2] - 1))]
file.KEGG.L3 <- '/home/drizzle_zhang/microbiome/result/9.PICRUSt/origin_data/KO_KEGG_L3.txt'
df.db.KEGG <- read.delim(file.KEGG.L3, row.names = 1)
file.meta.KEGG <- '/home/drizzle_zhang/microbiome/result/9.PICRUSt/origin_data/metadata_KEGG.txt'
df.meta.KEGG <- read.delim(file.meta.KEGG, row.names = 1)
df.meta.KEGG$KO_id <- row.names(df.meta.KEGG)

# gender
gender <- 'male'

# cutoff
type.cutoff <- 'fdr'

# dose
# vec.dose <- c(0, 1, 2, 3)
vec.dose <- c(0, 3)

# path
path.plot <- paste0('/home/drizzle_zhang/microbiome/result/9.PICRUSt/heatmap_',
                    gender, '_', type.cutoff)
# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.out.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)
series.time <- unique(df.meta$Time)

#### fold change matrix
for (sub.time in series.time) {
    file.res <- paste0(path.plot, "/OUT_edgeR_", sub.time, 
                       paste0(as.character(vec.dose), collapse = ''), 
                       ".txt")
    res.edgeR <- read.delim(file.res, row.names = 1)
    res.edgeR$KO_id <- row.names(res.edgeR)
    df.fc <- data.frame(res.edgeR[, c('KO_id', 'logFC')])
    names(df.fc) <- c('KO_id', sub.time)
    if (sub.time == 'A') {
        mat.fc <- df.fc
    } else {
        mat.fc <- merge(mat.fc, df.fc, by = 'KO_id')
    }
    
}
row.names(mat.fc) <- mat.fc$KO_id
mat.fc$KO_id <- NULL

# abundance
level = 'family'
path.plot <- paste0('/home/drizzle_zhang/microbiome/result/3.Community_Structure/barplot_',
                    level)
file.bubble <- paste0(path.plot, '/', gender, '_Bubble_',
                      paste0(as.character(vec.dose), collapse = ''), '.txt')
df.bubble <- read.delim(file.bubble)
sel.family <- 'Lactobacillaceae'
sel.bubble <- df.bubble[df.bubble$type == sel.family,]
df.corr <- data.frame()
sel.bubble.fc <- sel.bubble$fold.change
for (KO_id in row.names(mat.fc)) {
    corr.spearman <- cor(sel.bubble.fc, as.vector(t(mat.fc[KO_id,])), method = 'spearman')
    df.corr <- rbind(df.corr, data.frame(KO_id = KO_id, corr_spearman = corr.spearman))
}
df.corr.KEGG <- merge(df.corr, df.meta.KEGG, by = 'KO_id')

