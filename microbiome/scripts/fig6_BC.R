library(ggplot2)
library(clusterProfiler)
library(foreach)
library(doParallel)
library(scales)
library(patchwork)
registerDoParallel(4)

file.KEGG <- '/home/drizzle_zhang/microbiome/result/9.PICRUSt/ko_predict/ko_predictions_combine.txt'
df.KEGG <- read.delim(file.KEGG, row.names = 1)

# KEGG
file.KEGG.L3 <- '/home/drizzle_zhang/microbiome/result/9.PICRUSt/origin_data/KO_KEGG_L3.txt'
df.db.KEGG.L3 <- read.delim(file.KEGG.L3, row.names = 1)
file.KEGG.L2 <- '/home/drizzle_zhang/microbiome/result/9.PICRUSt/origin_data/KO_KEGG_L2.txt'
df.db.KEGG.L2 <- read.delim(file.KEGG.L2, row.names = 1)
file.meta.KEGG <- '/home/drizzle_zhang/microbiome/result/9.PICRUSt/origin_data/metadata_KEGG.txt'
df.meta.KEGG <- read.delim(file.meta.KEGG, row.names = 1)
df.meta.KEGG$KO_id <- row.names(df.meta.KEGG)
file.gmt.L2 <- '/home/drizzle_zhang/microbiome/result/9.PICRUSt/origin_data/KEGG_KO_L2.gmt'
df.gene.set.L2 <- read.gmt(file.gmt.L2)
file.gmt.L3 <- '/home/drizzle_zhang/microbiome/result/9.PICRUSt/origin_data/KEGG_KO_L3.gmt'
df.gene.set.L3 <- read.gmt(file.gmt.L3)

# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.out.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)
# gender
gender <- 'female'
df.meta.gender <- df.meta[df.meta$Gender == gender, ]

# cutoff
type.cutoff <- 'fdr'

# dose
vec.dose <- c(0, 1, 2, 3)
# vec.dose <- c(0, 3)

# edgeR and heatmap
find.sig.KO <- function(df.meta.gender, df.db.KEGG.L2, df.db.KEGG.L3, 
                        df.gene.set.L2, df.gene.set.L3, 
                        path.plot, vec.dose, type.cutoff, sub.time) {
    library(edgeR)
    library(gplots)
    library(RColorBrewer)
    sel.meta <- df.meta.gender[df.meta.gender$Time == sub.time,]
    sel.meta <- sel.meta[sel.meta$Dose %in% vec.dose,]
    use.sample <- sel.meta$SampleName
    
    sub.KEGG <- df.KEGG[, use.sample]
    
    factor.group <- as.factor(sel.meta$Group)
    d = DGEList(counts=sub.KEGG, group=factor.group)
    d = calcNormFactors(d)
    
    # 生成实验设计矩阵
    design.mat = model.matrix(~ 0 + d$samples$group)
    dimnames(design.mat)[[2]] <- levels(factor.group)
    d2 = estimateGLMCommonDisp(d, design.mat)
    d2 = estimateGLMTagwiseDisp(d2, design.mat)
    fit = glmFit(d2, design.mat)
    
    # 设置比较组
    BvsA <- makeContrasts(
        contrasts = paste(levels(factor.group), collapse = '-'),
        levels=design.mat)
    # 组间比较,统计Fold change, Pvalue
    lrt = glmLRT(fit,contrast=BvsA)
    # FDR检验，控制假阳性率小于5%
    de_lrt = decideTestsDGE(lrt, adjust.method="fdr", p.value=0.05)
    
    # 导出计算结果
    res.edgeR=lrt$table
    res.edgeR$sig.edger=de_lrt
    vec.sig <- rep(0, dim(res.edgeR)[1])
    vec.sig[(res.edgeR$logFC > 1) & (res.edgeR$PValue < 0.05)] <- 1
    vec.sig[(res.edgeR$logFC < -1) & (res.edgeR$PValue < 0.05)] <- -1
    res.edgeR$sig <- vec.sig
    if (type.cutoff == 'fdr') {
        enriched = row.names(subset(res.edgeR, sig.edger==1))
        depleted = row.names(subset(res.edgeR, sig.edger==-1))
    } else {
        enriched = row.names(subset(res.edgeR, sig==1))
        depleted = row.names(subset(res.edgeR, sig==-1))
    }
    
    # write results
    file.res <- paste0(path.plot, "/OUT_edgeR_", sub.time, 
                       paste0(as.character(vec.dose), collapse = ''), 
                       ".txt")
    res.edgeR <- res.edgeR[order(res.edgeR$logFC), ]
    write.table(res.edgeR, file = file.res, quote = F, sep = '\t')
    
    # 热图展示差异OTU
    # Sig OTU in two genotype
    DE <- c(enriched,depleted)
    sub_norm <- as.matrix(sub.KEGG[DE, ])
    
    if (dim(sub_norm)[1] > 1) {
        pdf(file = paste0(path.plot, "/heatmap_", sub.time, 
                          paste0(as.character(vec.dose), collapse = ''), 
                          ".pdf"),
            height = 10, width = 15)
        heatmap.2(sub_norm, scale="row", 
                  Colv=F, Rowv=F, dendrogram="none",
                  col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)),
                  key = T, keysize=1,density.info="none",
                  main=NULL,
                  margins = c(3, 4),
                  trace="none")
        dev.off()
    } else {
        print('Matrix of heatmap have only 1 rows')
        print(dim(sub_norm)[1])
    }
    
    ##### enrichment analysis
    levels <- c('L2', 'L3')
    for (level in levels) {
        if (level == 'L2') {
            df.db.KEGG <- df.db.KEGG.L2
            df.gene.set <- df.gene.set.L2
            p.cutoff <- 0.05
        }
        if (level == 'L3') {
            df.db.KEGG <- df.db.KEGG.L3
            df.gene.set <- df.gene.set.L3
            p.cutoff <- 0.01
        }
        # depleted in treat
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
            # write enrichment results
            file.control <- paste0(path.plot, "/Enrich_", level, "_control_", sub.time, 
                                   paste0(as.character(vec.dose), collapse = ''), 
                                   ".txt")
            write.table(df.enrich, file = file.control, quote = F, sep = '\t')
            
            df.enrich.sig <- df.enrich[df.enrich$pvalue < p.cutoff,]
            if (dim(df.enrich.sig)[1] > 1) {
                df.enrich.sig$lgp <- -log10(df.enrich.sig$pvalue)
                plot.deplete <- 
                    ggplot(data = df.enrich.sig, aes(x = pathway, y = lgp)) + 
                    geom_bar(stat = 'identity') +
                    labs(y = '-log(Pvalue)', x = 'Pathway', title = 'Enriched in Control') + 
                    coord_flip()
                ggsave(filename = paste0("/pathway_", level, "_treat_deplete_", sub.time, 
                                         paste0(as.character(vec.dose), collapse = ''), ".png"),
                       path = path.plot, plot = plot.deplete)
            }
        }
        
        # enriched in treat
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
            # write enrichment results
            file.treat <- paste0(path.plot, "/Enrich_", level, "_treatment_", sub.time, 
                                 paste0(as.character(vec.dose), collapse = ''), 
                                 ".txt")
            write.table(df.enrich, file = file.treat, quote = F, sep = '\t')
            
            df.enrich.sig <- df.enrich[df.enrich$pvalue < p.cutoff,]
            if (dim(df.enrich.sig)[1] > 1) {
                df.enrich.sig$lgp <- -log10(df.enrich.sig$pvalue)
                plot.enrich <- 
                    ggplot(data = df.enrich.sig, aes(x = pathway, y = lgp)) + 
                    geom_bar(stat = 'identity') +
                    labs(y = '-log(Pvalue)', x = 'Pathway', title = 'Enriched in treat') + 
                    coord_flip()
                ggsave(filename = paste0("/pathway_", level, "_treat_enrich_", sub.time, 
                                         paste0(as.character(vec.dose), collapse = ''), ".png"),
                       path = path.plot, plot = plot.enrich)
            }
        }
        # GSEA
        # filter results
        res.edgeR.L <- res.edgeR
        res.edgeR.L$KO_id <- row.names(res.edgeR.L)
        df.counts.mean <- data.frame(KO_id = row.names(sub.KEGG),
                                     mean_counts = rowMeans(sub.KEGG))
        res.edgeR.L <- merge(res.edgeR.L, df.counts.mean, by = 'KO_id')
        row.names(res.edgeR.L) <- res.edgeR.L$KO_id
        res.edgeR.L <- res.edgeR.L[res.edgeR.L$mean_counts > 10,]

        res.edgeR.L$logPval <- log10(res.edgeR.L$PValue) * 
            (res.edgeR.L$logFC / abs(res.edgeR.L$logFC))
        geneList <- res.edgeR.L$logPval
        names(geneList) <- row.names(res.edgeR.L)
        geneList[is.na(geneList)] <- 0
        geneList <- geneList[order(geneList, decreasing = T)]
        egmt <- GSEA(geneList, TERM2GENE = df.gene.set, pvalueCutoff = 0.9)
        res.egmt <- egmt@result
        file.GSEA <- paste0(path.plot, "/GSEA_", level, "_", sub.time, "_",
                            paste0(as.character(vec.dose), collapse = ''),
                            ".txt")
        write.table(res.egmt, file = file.GSEA, quote = F, sep = '\t', row.names = F)
    }
}

# time series
path.plot <- paste0('/home/drizzle_zhang/microbiome/result/9.PICRUSt/heatmap_',
                    gender, '_', type.cutoff)
if (!file.exists(path.plot)) {
    dir.create(path.plot)
}
series.time <- unique(df.meta$Time)
files.res <- foreach(sub.time = series.time, .combine = rbind) %dopar% 
    find.sig.KO(df.meta.gender, df.db.KEGG.L2, df.db.KEGG.L3, 
                df.gene.set.L2, df.gene.set.L3, 
                path.plot, vec.dose, type.cutoff, sub.time)


# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.txt'
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
vec.dose <- c(0, 1, 2, 3)
# vec.dose <- c(0, 1)

# time series
series.time <- unique(df.meta$Time)

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
        df.sort <- rbind(df.sort, data.frame(pathway = row, time = col,
                                             value = df.GSEA[row, as.character(col)],
                                             stringsAsFactors = F))
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
df.ks.male.filter <- df.ks.male[
    df.ks.male$qvalue.control < 0.1 | df.ks.male$qvalue.treat < 0.1,]
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
        panel.grid  = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 9, color = "black", family = 'Arial'), 
        panel.background = element_rect(color = 'white', size = 0,
                                        fill = 'transparent'),
        legend.text = element_text(size = 9)
    )

plot.final.male <- plot.male + plot.heatmap.male + plot_layout(widths = c(1, 1.6),
                                                      guides = 'collect')

ggsave(plot = plot.final.male, path = path.out, 
       filename = paste0("/Male_GSEA_",  
                         paste0(as.character(vec.dose), collapse = ''), ".png"),
       height = 10, width = 25, units = 'cm')





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
        df.sort <- rbind(df.sort, data.frame(pathway = row, time = col,
                                             value = df.GSEA[row, as.character(col)],
                                             stringsAsFactors = F))
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
df.ks.female.filter <- df.ks.female[
    df.ks.female$qvalue.control < 0.1 | df.ks.female$qvalue.treat < 0.1,]
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
    theme(panel.background = element_rect(color = 'transparent', size = 1.5,
                                          fill = 'transparent'),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 9, color = "black", family = 'Arial'), 
          legend.text = element_text(size = 9))

plot.final.female <- plot.female + plot.heatmap.female + 
    plot_layout(widths = c(1, 1.6), guides = 'collect')

ggsave(plot = plot.final.female, path = path.out, 
       filename = paste0("/Female_GSEA_",  
                         paste0(as.character(vec.dose), collapse = ''), ".png"),
       height = 7, width = 25, units = 'cm')

