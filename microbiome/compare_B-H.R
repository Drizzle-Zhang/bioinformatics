library(ggplot2)

# KEGG
file.KEGG <- '/home/drizzle_zhang/microbiome/result/9.PICRUSt/ko_predict/ko_predictions.txt'
df.KEGG <- read.delim(file.KEGG, row.names = 1)
file.KEGG.L3 <- '/home/drizzle_zhang/microbiome/result/9.PICRUSt/origin_data/KO_KEGG_L3.txt'
df.db.KEGG <- read.delim(file.KEGG.L3, row.names = 1)

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

### sig matrix
for (sub.time in series.time) {
    file.res <- paste0(path.plot, "/OUT_edgeR_", sub.time, 
                       paste0(as.character(vec.dose), collapse = ''), 
                       ".txt")
    res.edgeR <- read.delim(file.res, row.names = 1)
    res.edgeR$KO_id <- row.names(res.edgeR)
    df.sig <- data.frame(res.edgeR[, c('KO_id', 'sig.edger')])
    names(df.sig) <- c('KO_id', sub.time)
    if (sub.time == 'A') {
        mat.sig <- df.sig
    } else {
        mat.sig <- merge(mat.sig, df.sig, by = 'KO_id')
    }
    
}
row.names(mat.sig) <- mat.sig$KO_id
mat.sig$KO_id <- NULL

# male
if (gender == 'male') {
    # sel.time <- c('B', 'C', 'D', 'E', 'F', 'G', 'H')
    sel.time <- c('B', 'C', 'D', 'E', 'F', 'G')
    mat.sig.sel <- mat.sig[, sel.time]
    mat.sig.sig <- mat.sig.sel[
        abs(rowSums(mat.sig.sel)) >= length(sel.time) - 1, ]
    
}
# female
if (gender == 'female') {
    sel.time <- c('B', 'C', 'D', 'E', 'F', 'G')
    mat.sig.sel <- mat.sig[, sel.time]
    mat.sig.sig <- mat.sig.sel[
        abs(rowSums(mat.sig.sel)) >= length(sel.time) - 1, ]
}

########## enrichment analysis
# depleted in treat
enriched <- row.names(mat.sig.sig[rowSums(mat.sig.sig) > 0, ])
vec.sig <- rep(0, dim(df.KEGG)[1])
names(vec.sig) <- row.names(df.KEGG)
for (KO_id in row.names(df.KEGG)) {
    if (KO_id %in% enriched) {
        vec.sig[KO_id] <- 1
    }
}
if (sum(vec.sig) > 0) {
    df.enrich <- data.frame()
    for (pathway in names(df.db.KEGG)) {
        df.fisher <- data.frame(sig = vec.sig, 
                                pathway = df.db.KEGG[,pathway])
        table.fisher <- xtabs(~ sig + pathway, data = df.fisher)
        out.fisher <- fisher.test(table.fisher)
        df.enrich <- rbind(df.enrich, 
                           data.frame(pathway = pathway, 
                                      pvalue = out.fisher$p.value))
    }
    df.enrich$qvalue <- p.adjust(df.enrich$pvalue)
    df.enrich.sig <- df.enrich[df.enrich$pvalue < 0.05,]
    if (dim(df.enrich.sig)[1] > 1) {
        df.enrich.sig$lgp <- -log10(df.enrich.sig$pvalue)
        plot.deplete <- 
            ggplot(data = df.enrich.sig, aes(x = pathway, y = lgp)) + 
            geom_bar(stat = 'identity') +
            labs(y = '-log(Pvalue)', x = 'Pathway', title = 'Enriched in Control') + 
            coord_flip()
        # ggsave(filename = paste0("/Combine_pathway_treat_deplete_", sub.time, 
        #                          paste0(as.character(vec.dose), collapse = ''), ".png"),
        #        path = path.plot, plot = plot.deplete)
    }
}

# enriched in treat
depleted <- row.names(mat.sig.sig[rowSums(mat.sig.sig) < 0, ])
vec.sig <- rep(0, dim(df.KEGG)[1])
names(vec.sig) <- row.names(df.KEGG)
for (KO_id in row.names(df.KEGG)) {
    if (KO_id %in% depleted) {
        vec.sig[KO_id] <- 1
    }
}
if (sum(vec.sig) > 0) {
    df.enrich <- data.frame()
    for (pathway in names(df.db.KEGG)) {
        df.fisher <- data.frame(sig = vec.sig, 
                                pathway = df.db.KEGG[,pathway])
        table.fisher <- xtabs(~ sig + pathway, data = df.fisher)
        out.fisher <- fisher.test(table.fisher)
        df.enrich <- rbind(df.enrich, 
                           data.frame(pathway = pathway, 
                                      pvalue = out.fisher$p.value))
    }
    df.enrich$qvalue <- p.adjust(df.enrich$pvalue)
    df.enrich.sig <- df.enrich[df.enrich$pvalue < 0.05,]
    if (dim(df.enrich.sig)[1] > 1) {
        df.enrich.sig$lgp <- -log10(df.enrich.sig$pvalue)
        plot.enrich <- 
            ggplot(data = df.enrich.sig, aes(x = pathway, y = lgp)) + 
            geom_bar(stat = 'identity') +
            labs(y = '-log(Pvalue)', x = 'Pathway', title = 'Enriched in treat') + 
            coord_flip()
        # ggsave(filename = paste0("/Combine_pathway_treat_enrich_", sub.time, 
        #                          paste0(as.character(vec.dose), collapse = ''), ".png"),
        #        path = path.plot, plot = plot.enrich)
    }
}



################ comparation of degree of diff
# gender
# male
gender <- 'male'
path.plot <- paste0('/home/drizzle_zhang/microbiome/result/9.PICRUSt/heatmap_',
                    gender, '_', type.cutoff)
### sig matrix
for (sub.time in series.time) {
    file.res <- paste0(path.plot, "/OUT_edgeR_", sub.time, 
                       paste0(as.character(vec.dose), collapse = ''), 
                       ".txt")
    res.edgeR <- read.delim(file.res, row.names = 1)
    res.edgeR$KO_id <- row.names(res.edgeR)
    df.sig <- data.frame(res.edgeR[, c('KO_id', 'sig.edger')])
    names(df.sig) <- c('KO_id', sub.time)
    if (sub.time == 'A') {
        mat.sig <- df.sig
    } else {
        mat.sig <- merge(mat.sig, df.sig, by = 'KO_id')
    }
    
}
row.names(mat.sig) <- mat.sig$KO_id
mat.sig$KO_id <- NULL

if (gender == 'male') {
    # sel.time <- c('B', 'C', 'D', 'E', 'F', 'G', 'H')
    sel.time <- c('B', 'C', 'D', 'E', 'F', 'G')
    mat.sig.sel <- mat.sig[, sel.time]
    mat.sig.sig <- mat.sig.sel[
        abs(rowSums(mat.sig.sel)) >= length(sel.time) - 1, ]
    
}
mat.sig.male <- mat.sig.sig
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
mat.fc.male <- mat.fc


# female
gender <- 'female'
path.plot <- paste0('/home/drizzle_zhang/microbiome/result/9.PICRUSt/heatmap_',
                    gender, '_', type.cutoff)
### sig matrix
for (sub.time in series.time) {
    file.res <- paste0(path.plot, "/OUT_edgeR_", sub.time, 
                       paste0(as.character(vec.dose), collapse = ''), 
                       ".txt")
    res.edgeR <- read.delim(file.res, row.names = 1)
    res.edgeR$KO_id <- row.names(res.edgeR)
    df.sig <- data.frame(res.edgeR[, c('KO_id', 'sig.edger')])
    names(df.sig) <- c('KO_id', sub.time)
    if (sub.time == 'A') {
        mat.sig <- df.sig
    } else {
        mat.sig <- merge(mat.sig, df.sig, by = 'KO_id')
    }
    
}
row.names(mat.sig) <- mat.sig$KO_id
mat.sig$KO_id <- NULL

if (gender == 'male') {
    # sel.time <- c('B', 'C', 'D', 'E', 'F', 'G', 'H')
    sel.time <- c('B', 'C', 'D', 'E', 'F', 'G')
    mat.sig.sel <- mat.sig[, sel.time]
    mat.sig.sig <- mat.sig.sel[
        abs(rowSums(mat.sig.sel)) >= length(sel.time) - 1, ]
    
}
mat.sig.female <- mat.sig.sig
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
mat.fc.female <- mat.fc

## compare 
overlap.KO <- intersect(row.names(mat.sig.female), row.names(mat.sig.male))
fc.male <- rowMeans(mat.fc.male[overlap.KO,])
fc.female <- rowMeans(mat.fc.female[overlap.KO,])
wilcox.test(fc.male, fc.female, alternative = 'greater')
# plot
library(ggsignif)
df.plot <- data.frame(FC = fc.male, gender = rep('male', length(fc.male)))
df.plot <- rbind(df.plot, 
                 data.frame(FC = fc.female, gender = rep('female', length(fc.female))))
compaired <- list(c("male", "female"))
ggplot(df.plot, aes(gender, FC, fill = gender)) +
    geom_boxplot(width = 0.5) +
    theme(
        plot.title = element_text(size = 25),
        axis.text.x = element_text(size = 10, angle = 0),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)
    ) + labs(x = 'Gender', y = 'Log10(Fold Change)') + geom_signif(
        comparisons = compaired,
        step_increase = 0.1,
        map_signif_level = F,
        test = 'wilcox.test', alternative = 'greater'
    )


