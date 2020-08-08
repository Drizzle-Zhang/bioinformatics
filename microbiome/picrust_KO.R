# 使用edgeR统计组间差异OTU

library(ggplot2)
library(clusterProfiler)
library(foreach)
library(doParallel)
registerDoParallel(4)
# create DGE list
# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.out.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)

file.KEGG <- '/home/drizzle_zhang/microbiome/result/9.PICRUSt/ko_predict/ko_predictions.txt'
df.KEGG <- read.delim(file.KEGG, row.names = 1)
df.KEGG <- df.KEGG[ ,paste0('S', 1:(dim(df.KEGG)[2] - 1))]
# df.KEGG.scale <- scale(df.KEGG)

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

# gender
gender <- 'female'
df.meta.gender <- df.meta[df.meta$Gender == gender, ]

# cutoff
type.cutoff <- 'fdr'

# dose
# vec.dose <- c(0, 1, 2, 3)
vec.dose <- c(0, 3)

# edgeR and heatmap
find.sig.KO <- function(df.meta.gender, df.db.KEGG.L2, df.db.KEGG.L3, 
                        df.gene.set.L2, df.gene.set.L3, 
                        path.plot, vec.dose, type.cutoff, sub.time) {
    library(edgeR)
    library(gplots)
    library(RColorBrewer)
    sel.meta <- df.meta.gender[df.meta.gender$Time == sub.time,]
    sel.meta <- sel.meta[sel.meta$Dose %in% vec.dose,]
    use.sample <- sel.meta$Sample

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
    # bool.fc <- rep(0, dim(res.edgeR)[1])
    # bool.fc[abs(res.edgeR$logFC) > 1.5] <- 1
    # res.edgeR$bool.fc <- bool.fc
    file.res <- paste0(path.plot, "/OUT_edgeR_", sub.time, 
                       paste0(as.character(vec.dose), collapse = ''), 
                       ".txt")
    # res.edgeR$qvalue <- p.adjust(res.edgeR$PValue, method = 'fdr')
    res.edgeR <- res.edgeR[order(res.edgeR$logFC), ]
    write.table(res.edgeR, file = file.res, quote = F, sep = '\t')

    # 热图展示差异OTU
    # pair_group = subset(sub_design, genotype %in% c("OE", "WT"))
    # Sig OTU in two genotype
    DE <- c(enriched,depleted)
    sub_norm <- as.matrix(sub.KEGG[DE, ])

    if (dim(sub_norm)[1] > 1) {
        #colnames(sub_norm)=gsub("DM","KO",colnames(sub_norm),perl=TRUE) # rename samples ID
        pdf(file = paste0(path.plot, "/heatmap_", sub.time, 
                          paste0(as.character(vec.dose), collapse = ''), 
                          ".pdf"),
            height = 10, width = 15)
        heatmap.2(sub_norm, scale="row", 
                  Colv=F, Rowv=F, dendrogram="none",
                  col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)),
                  # cexCol=1,
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
        # res.edgeR.L$pvalue <- res.edgeR.L$PValue
        # res.edgeR.L$pvalue[res.edgeR.L$mean_counts < 5] <- 1
        res.edgeR.L <- res.edgeR.L[res.edgeR.L$mean_counts > 5,]
        # res.edgeR.L$qvalue <- p.adjust(res.edgeR.L$PValue, method = 'fdr')
        
        res.edgeR.L$logPval <- -log10(res.edgeR.L$PValue) * 
            (res.edgeR.L$logFC / abs(res.edgeR.L$logFC))
        geneList <- res.edgeR.L$logPval
        names(geneList) <- row.names(res.edgeR.L)
        geneList[is.na(geneList)] <- 0
        geneList <- geneList[order(geneList, decreasing = T)]
        egmt <- GSEA(geneList, TERM2GENE = df.gene.set, pvalueCutoff = 0.5)
        res.egmt <- egmt@result
        # vec.KO <- strsplit(
        #     res.egmt['Photosynthesis - antenna proteins', 'core_enrichment'], '/')[[1]]
        # res.edgeR.L[vec.KO,]
        file.GSEA <- paste0(path.plot, "/GSEA_", level, "_", sub.time, 
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

#################### fold change matrix
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
mat.fc.KEGG <- merge(mat.fc, df.meta.KEGG, by = 'KO_id')
row.names(mat.fc) <- mat.fc$KO_id
mat.fc$KO_id <- NULL
# cluster
library(Seurat)
# data preparing
seurat.fc <- CreateSeuratObject(counts = mat.fc, project = "pbmc3k", min.cells = 0, min.features = 0)
seurat.fc <- FindVariableFeatures(seurat.fc, selection.method = "vst", 
                                  nfeatures = 500)
# VariableFeaturePlot(seurat.fc)
all.genes <- rownames(seurat.fc)
seurat.fc <- ScaleData(seurat.fc, features = all.genes)

# add label
seurat.fc@meta.data$original.label <- series.time

# PCA
seurat.fc <- 
    RunPCA(seurat.fc, features = VariableFeatures(object = seurat.fc), 
           verbose = F, npcs = 5)
DimPlot(seurat.fc, reduction = "pca", label = T, group.by = 'original.label')


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
mat.sig.KEGG <- merge(mat.sig, df.meta.KEGG, by = 'KO_id')
row.names(mat.sig) <- mat.sig$KO_id
mat.sig$KO_id <- NULL
# mat.sig.sel <- mat.sig
mat.sig.sel <- mat.sig[rowSums(mat.sig) > 0,]
d.sig <- dist(t(mat.sig.sel), method = 'binary')
fit.sig <- hclust(d.sig, method = 'single')
plot(fit.sig)
pcoa <- cmdscale(d.sig, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points <- as.data.frame(pcoa$points) # get coordinate string, format to dataframme
eig <- pcoa$eig
colnames(points) <- c("x", "y", "z")
points <- cbind(points, series.time)
ggplot(points, aes(x = x, y = z)) +
    labs(x=paste("PCo 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
         y=paste("PCo 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
    geom_point(alpha=.7, size=2) + 
    geom_text(aes(label = series.time)) + 
    theme_classic() + 
    theme(text=element_text(family="sans", size=7))

# male
if (gender == 'male') {
    sel.time <- c('B', 'C', 'D', 'E', 'F', 'G', 'H', 'I')
    mat.sig.sel <- mat.sig[, sel.time]
    mat.sig.sig <- mat.sig.sel[
        abs(rowSums(mat.sig.sel)) >= length(sel.time) - 3, ]
    
}
# female
if (gender == 'female') {
    sel.time <- c('B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M')
    mat.sig.sel <- mat.sig[, sel.time]
    mat.sig.sig <- mat.sig.sel[
        abs(rowSums(mat.sig.sel)) >= length(sel.time) - 5, ]
}


########## heatmap
library(gplots)
library(RColorBrewer)
df.meta.heat <- df.meta.gender[df.meta.gender$Time %in% sel.time, ]
df.meta.heat <- df.meta.heat[df.meta.heat$Dose %in% c(0, 3), ]
# df.meta.heat <- df.meta.heat[order(df.meta.heat$Group),]
use.sample <- df.meta.heat$Sample
sub.KEGG.heat <- df.KEGG[row.names(mat.sig.sig), use.sample]
i = 1
for (sub.time in sel.time) {
    sub.meta <- df.meta.heat[df.meta.heat$Time == sub.time,]
    sub.sample.control <- sub.meta$Sample[sub.meta$Group == 'Control']
    median.time <- apply(sub.KEGG.heat[, sub.sample.control], 1, median)
    sub.sample <- sub.meta$Sample
    sub.heat <- sub.KEGG.heat[, sub.sample] / median.time
    if (i == 1) {
        sub.KEGG.heat.nrom <- sub.heat
    } else {
        sub.KEGG.heat.nrom <- cbind(sub.KEGG.heat.nrom, sub.heat)
    }
    i = i + 1
}
pdf(file=paste0(path.plot, "/Combine_heatmap.pdf"),
    height = 8, width = 15)
heatmap.2(as.matrix(sub.KEGG.heat.nrom), scale="row", 
          Colv=F, Rowv=F, dendrogram="none",
          col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)),
          # cexCol=1,
          key = T, keysize=1,density.info="none",
          main=NULL,
          margins = c(3, 5),
          trace="none")
dev.off()

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
        ggsave(filename = paste0("/Combine_pathway_treat_deplete_", sub.time, 
                                 paste0(as.character(vec.dose), collapse = ''), ".png"),
               path = path.plot, plot = plot.deplete)
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
        ggsave(filename = paste0("/Combine_pathway_treat_enrich_", sub.time, 
                                 paste0(as.character(vec.dose), collapse = ''), ".png"),
               path = path.plot, plot = plot.enrich)
    }
}


########## find diff KO
file.family.KO <- 
    '/home/drizzle_zhang/microbiome/result/9.PICRUSt/ko_predict/select_KOs/mat_family_KO.txt'
df.family.KO <- read.delim(file.family.KO, row.names = 1)
sel.family <- 'f__Lactobacillaceae'
bool.group <- rep('control', dim(df.family.KO)[2])
bool.group[names(df.family.KO) == sel.family] <- 'case'
factor.group <- as.factor(bool.group)
d = DGEList(counts=df.family.KO, group=factor.group)
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
    levels=design.mat)# 组间比较,统计Fold change, Pvalue
lrt = glmLRT(fit,contrast=BvsA)
# FDR检验，控制假阳性率小于5%
de_lrt = decideTestsDGE(lrt, adjust.method="fdr", p.value=0.05)
res.edgeR <- lrt$table
res.edgeR$sig=de_lrt

sig.edgeR <- res.edgeR[(res.edgeR$logFC > 1.5) & (res.edgeR$PValue < 0.01), ]
vec.sig <- rep(0, dim(df.KEGG)[1])
names(vec.sig) <- row.names(df.KEGG)
for (KO_id in row.names(df.KEGG)) {
    if (KO_id %in% row.names(sig.edgeR)) {
        vec.sig[KO_id] <- 1
    }
}
file.KEGG.L3 <- '/home/drizzle_zhang/microbiome/result/9.PICRUSt/origin_data/KO_KEGG_L3.txt'
df.db.KEGG.L3 <- read.delim(file.KEGG.L3, row.names = 1)
df.enrich <- data.frame()
for (pathway in names(df.db.KEGG.L3)) {
    df.fisher <- data.frame(sig = vec.sig, 
                            pathway = df.db.KEGG.L3[,pathway])
    table.fisher <- xtabs(~ sig + pathway, data = df.fisher)
    out.fisher <- fisher.test(table.fisher)
    df.enrich <- rbind(df.enrich, 
                       data.frame(pathway = pathway, 
                                  pvalue = out.fisher$p.value))
}
df.enrich.sig <- df.enrich[df.enrich$pvalue < 0.05,]
df.enrich.sig$lgp <- -log10(df.enrich.sig$pvalue)
ggplot(data = df.enrich.sig, aes(x = pathway, y = lgp)) + 
    geom_bar(stat = 'identity') +
    labs(y = '-log(Pvalue)', x = 'Pathway') + 
    coord_flip()


########## associate bactoria
file.sub.mat <- 
    '/home/drizzle_zhang/microbiome/result/9.PICRUSt/ko_predict/select_KOs/mat_OTU_KO_f__Lactobacillaceae.txt'
df.sub.mat <- read.delim(file.sub.mat, row.names = 1)
vec.sub <- colSums(df.sub.mat)
vec.sig <- vec.sub[row.names(mat.sig.sig)]
wilcox.test(vec.sig, vec.sub, alternative = 'greater')

