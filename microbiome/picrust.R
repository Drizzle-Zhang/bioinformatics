# 使用edgeR统计组间差异OTU

library(edgeR)
library(gplots)
library(RColorBrewer)
# create DGE list
# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.out.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)

file.KEGG <- '/home/drizzle_zhang/microbiome/result/9.PICRUSt/ko_predict/KEGG_L3.tsv.txt'
df.KEGG <- read.delim(file.KEGG, row.names = 1)
# df.KEGG.scale <- scale(df.KEGG)

gender <- 'male'
df.meta <- df.meta[df.meta$Gender == gender, ]

# time series
path.plot <- paste0('/home/drizzle_zhang/microbiome/result/9.PICRUSt/heatmap_',
                    gender)
series.time <- unique(df.meta$Time)
df.fc.pval <- data.frame()
for (sub.time in series.time) {
    # sub.time <- 'D'
    sel.meta <- df.meta[df.meta$Time == sub.time,]
    sel.meta <- sel.meta[sel.meta$Dose %in% c(0, 3),]
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
        res.edgeR$sig=de_lrt
    enriched = row.names(subset(res.edgeR, sig==1))
    depleted = row.names(subset(res.edgeR, sig==-1))
    
    res.edgeR <- res.edgeR[order(res.edgeR$logFC),]
    file.res <- paste0(path.plot, "/pathway_", sub.time, ".txt")
    write.table(res.edgeR, file = file.res, quote = F)
    
    # combine results
    # for (term in row.names(res.edgeR)) {
    #     df.fc.pval <- rbind(df.fc.pval, 
    #                         data.frame(Term = term, Time = sub.time,
    #                                    FoldChange = res.edgeR[term, 'logFC'],
    #                                    PValue = res.edgeR[term, 'PValue']))
    # }
    # 
    # if (sub.time == 'A') {
    #     mat.fc <- data.frame(res.edgeR[,'logFC'], 
    #                          row.names = row.names(res.edgeR))
    # } else {
    #     mat.fc <- cbind(mat.fc,
    #                     data.frame(res.edgeR[,'logFC'], 
    #                                row.names = row.names(res.edgeR)))
    # }
    
    ## 热图展示差异OTU
    # pair_group = subset(sub_design, genotype %in% c("OE", "WT"))
    # Sig OTU in two genotype
    DE=c(enriched,depleted)
    sub_norm = as.matrix(sub.KEGG[DE, ])

    if (dim(sub_norm)[1] < 2) {next()}
    #colnames(sub_norm)=gsub("DM","KO",colnames(sub_norm),perl=TRUE) # rename samples ID
    pdf(file=paste0(path.plot, "/pathway_heatmap_", sub.time, ".pdf"), 
        height = 10, width = 15)
    heatmap.2(sub_norm, scale="row", Colv=F, Rowv=F,
              dendrogram="none",
              col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)),
              cexCol=1, cexRow = 1,
              key = T, keysize=1,density.info="density",
              main=NULL,
              trace="none",margins = c(3, 17))
    dev.off()
}

# names(mat.fc) <- series.time