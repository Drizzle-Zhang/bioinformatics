# 使用edgeR统计组间差异OTU

library(edgeR)
# create DGE list
# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.out.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)

file.KEGG <- '/home/drizzle_zhang/microbiome/result/9.PICRUSt/ko_predict/KEGG_L3.tsv.txt'
df.KEGG <- read.delim(file.KEGG, row.names = 1)

d = DGEList(counts=count, group=sub_design$genotype)
d = calcNormFactors(d)

# 生成实验设计矩阵
design.mat = model.matrix(~ 0 + d$samples$group)
colnames(design.mat)=levels(genotypes)
d2 = estimateGLMCommonDisp(d, design.mat)
d2 = estimateGLMTagwiseDisp(d2, design.mat)
fit = glmFit(d2, design.mat)

# 设置比较组
BvsA <- makeContrasts(contrasts = "OE-WT", levels=design.mat)
# 组间比较,统计Fold change, Pvalue
lrt = glmLRT(fit,contrast=BvsA)
# FDR检验，控制假阳性率小于5%
de_lrt = decideTestsDGE(lrt, adjust.method="fdr", p.value=0.05)

# 导出计算结果
x=lrt$table
x$sig=de_lrt
enriched = row.names(subset(x,sig==1))
depleted = row.names(subset(x,sig==-1))

## 热图展示差异OTU
pair_group = subset(sub_design, genotype %in% c("OE", "WT"))
# Sig OTU in two genotype
DE=c(enriched,depleted)
sub_norm = as.matrix(norm[DE, rownames(pair_group)])
#colnames(sub_norm)=gsub("DM","KO",colnames(sub_norm),perl=TRUE) # rename samples ID
pdf(file=paste("heat_otu_OEvsWT_sig.pdf", sep=""), height = 8, width = 8)
# scale in row, dendrogram only in row, not cluster in column
heatmap.2(sub_norm, scale="row", Colv=FALSE, Rowv=FALSE,dendrogram="none", col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)), cexCol=1,keysize=1,density.info="none",main=NULL,trace="none")
dev.off()