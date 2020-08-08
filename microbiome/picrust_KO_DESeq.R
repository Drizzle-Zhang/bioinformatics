# 使用edgeR统计组间差异OTU

library(ggplot2)
library(foreach)
library(doParallel)
registerDoParallel(4)
# create DGE list
# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.out.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)

file.KEGG <- '/home/drizzle_zhang/microbiome/result/9.PICRUSt/ko_predict/ko_predictions.txt'
df.KEGG <- read.delim(file.KEGG, row.names = 1)
# df.KEGG.scale <- scale(df.KEGG)

# KEGG
file.KEGG.L3 <- '/home/drizzle_zhang/microbiome/result/9.PICRUSt/origin_data/KO_KEGG_L3.txt'
df.db.KEGG <- read.delim(file.KEGG.L3, row.names = 1)

# gender
gender <- 'male'
df.meta.gender <- df.meta[df.meta$Gender == gender, ]

# dose
# vec.dose <- c(0, 1, 2, 3)
vec.dose <- c(0, 3)

# edgeR and heatmap
find.sig.KO <- function(df.meta, df.db.KEGG, path.plot, vec.dose, sub.time, p.cutoff = 0.05) {
    library(DESeq2)
    library(gplots)
    library(RColorBrewer)
    sel.meta <- df.meta[df.meta$Time == sub.time,]
    sel.meta <- sel.meta[sel.meta$Dose %in% vec.dose,]
    use.sample <- sel.meta$Sample

    sub.KEGG <- df.KEGG[, use.sample]
    
    # diff analysis by DESeq
    countdata <- as.matrix(sub.KEGG)
    factor.group <- as.factor(sel.meta$Group)
    coldata <- DataFrame(row.names = use.sample, factor.group = factor.group)
    ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                     colData = coldata,
                                     design = ~ 1 + factor.group)
    dds <- DESeq(ddsMat)
    res.dds <- results(dds)
    res.DESeq <- data.frame(res.dds@listData, row.names = res.dds@rownames)
    
    sig.deseq <- rep(0, dim(res.DESeq)[1])
    sig.deseq[(res.DESeq$log2FoldChange > 0) & (res.DESeq$padj < 0.05)] <- 1
    sig.deseq[(res.DESeq$log2FoldChange < 0) & (res.DESeq$padj < 0.05)] <- -1
    res.DESeq$sig.deseq <- sig.deseq
    vec.sig <- rep(0, dim(res.DESeq)[1])
    vec.sig[(res.DESeq$log2FoldChange > 1) & (res.DESeq$pvalue < 0.05)] <- 1
    vec.sig[(res.DESeq$log2FoldChange < -1) & (res.DESeq$pvalue < 0.05)] <- -1
    res.DESeq$sig <- vec.sig
    enriched = row.names(subset(res.DESeq, sig.deseq==1))
    depleted = row.names(subset(res.DESeq, sig.deseq==-1))

    # write results
    file.res <- paste0(path.plot, "/OUT_DESeq_", sub.time, 
                       paste0(as.character(vec.dose), collapse = ''), 
                       ".txt")
    res.DESeq <- res.DESeq[order(res.DESeq$log2FoldChange), ]
    write.table(res.DESeq, file = file.res, quote = F, sep = '\t')

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
        df.enrich.sig <- df.enrich[df.enrich$qvalue < p.cutoff,]
        if (dim(df.enrich.sig)[1] > 1) {
            df.enrich.sig$lgp <- -log10(df.enrich.sig$pvalue)
            plot.deplete <- 
                ggplot(data = df.enrich.sig, aes(x = pathway, y = lgp)) + 
                geom_bar(stat = 'identity') +
                labs(y = '-log(Pvalue)', x = 'Pathway', title = 'Enriched in Treat') + 
                coord_flip()
            ggsave(filename = paste0("/pathway_treat_enrich_", sub.time, 
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
        df.enrich.sig <- df.enrich[df.enrich$qvalue < p.cutoff,]
        if (dim(df.enrich.sig)[1] > 1) {
            df.enrich.sig$lgp <- -log10(df.enrich.sig$pvalue)
            plot.enrich <- 
                ggplot(data = df.enrich.sig, aes(x = pathway, y = lgp)) + 
                geom_bar(stat = 'identity') +
                labs(y = '-log(Pvalue)', x = 'Pathway', title = 'Enriched in Control') + 
                coord_flip()
            ggsave(filename = paste0("/pathway_control_enrich_", sub.time, 
                                     paste0(as.character(vec.dose), collapse = ''), ".png"),
                   path = path.plot, plot = plot.enrich)
        }
    }
    
    return(file.res)

}

# time series
path.plot <- paste0('/home/drizzle_zhang/microbiome/result/9.PICRUSt/DESeq_',
                    gender)
series.time <- unique(df.meta$Time)
p.cutoff = 0.1
files.res <- foreach(sub.time = series.time, .combine = rbind) %dopar% 
    find.sig.KO(df.meta.gender, df.db.KEGG, path.plot, vec.dose, sub.time, p.cutoff)

#################### fold change matrix
for (sub.time in series.time) {
    file.res <- paste0(path.plot, "/OUT_DESeq_", sub.time, 
                       paste0(as.character(vec.dose), collapse = ''), 
                       ".txt")
    res.edgeR <- read.delim(file.res, row.names = 1)
    if (sub.time == 'A') {
        mat.fc <- data.frame(res.edgeR[,'log2FoldChange'],
                              row.names = row.names(res.edgeR))
    } else {
        mat.fc <- cbind(mat.fc,
                         data.frame(res.edgeR[,'log2FoldChange'],
                                    row.names = row.names(res.edgeR)))
    }
    
}
names(mat.fc) <- series.time
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
    file.res <- paste0(path.plot, "/OUT_DESeq_", sub.time, 
                       paste0(as.character(vec.dose), collapse = ''), 
                       ".txt")
    res.edgeR <- read.delim(file.res, row.names = 1)
    if (sub.time == 'A') {
        mat.sig <- data.frame(res.edgeR[,'sig'],
                             row.names = row.names(res.edgeR))
    } else {
        mat.sig <- cbind(mat.sig,
                        data.frame(res.edgeR[,'sig'],
                                   row.names = row.names(res.edgeR)))
    }
    
}
names(mat.sig) <- series.time
mat.sig.sel <- mat.sig[rowSums(mat.sig) > 0,]
d.sig <- dist(t(mat.sig.sel), method = 'binary')
fit.sig <- hclust(d.sig, method = 'complete')
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
    sel.time <- c('B', 'I', 'J', 'K', 'L', 'M')
    mat.sig.sel <- mat.sig[, sel.time]
    mat.sig.sig <- mat.sig.sel[
        abs(rowSums(mat.sig.sel)) >= length(sel.time) - 3, ]
}


########## heatmap
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
pdf(file=paste0(path.plot, "/heatmap_combine.pdf"),
    height = 8, width = 8)
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
# build database of L3
vec.sig <- rep(0, dim(df.KEGG)[1])
names(vec.sig) <- row.names(df.KEGG)
for (KO_id in row.names(df.KEGG)) {
    if (KO_id %in% row.names(mat.sig.sig)) {
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

