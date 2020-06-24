library("DESeq2")
library(airway)

# import data
data("airway")
se <- airway

# design 参数为 formula，此处为cell和dex两个因素，~ cell + dex表示我们想控制cell研究dex的影响。
dds <- DESeqDataSet(se, design = ~ cell + dex) 
countdata <- assay(se)
class(countdata)
head(countdata)
coldata <- colData(se)
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ cell + dex)
# 将 untrt 定义为dex因素的第一水平，随后的foldchange 将采用 trt/untrt
dds$dex <- relevel(dds$dex, "untrt") 
dds <- DESeq(dds)
# 得到结果，可以根据padj来挑选合适的差异表达基因，log2FoldChange来确定基因上调还是下调
res.deseq <- results(dds)
idx.genes <- which((res.deseq$log2FoldChange > 2) & (res.deseq$padj < 0.01))
diff.genes <- res.deseq@rownames[idx.genes]

