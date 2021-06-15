library(Seurat)

file.MCA <- '/home/disk/scRef/MouseAtlas_SingleCell_Han2018/combinedMCA/MCA_combined_mouse_uniform.txt'
df.MCA <- read.table(file.MCA, header=T, row.names=1, sep='\t', check.name=F)
# a <- t(t(df.MCA)/colSums(df.MCA))*1000000
library(DESeq2)
coldata.MCA <- DataFrame(row.names = names(df.MCA))
obj.DESeq.MCA <- DESeqDataSetFromMatrix(countData = df.MCA, colData = coldata.MCA, 
                                        design = ~ 1)
fpm.MCA <- fpm(obj.DESeq.MCA, robust = T)

file.sample <- '/home/disk/scRef/MouseAtlas_SingleCell_Han2018/combinedMCA/sample.txt'
df.sample <- read.table(file.sample, header=T, row.names=1, sep='\t', check.name=F)

seurat.unlabeled <- CreateSeuratObject(counts = df.sample, project = "MCA", min.cells = 0, min.features = 1000)
seurat.unlabeled <- NormalizeData(seurat.unlabeled, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.unlabeled <- FindVariableFeatures(seurat.unlabeled, selection.method = "mvp", 
                                         mean.cutoff = c(1, Inf),
                                         dispersion.cutoff = c(-0.1, 0.1))
VariableFeaturePlot(seurat.unlabeled)
all.genes <- rownames(seurat.unlabeled)
seurat.unlabeled <- ScaleData(seurat.unlabeled, features = all.genes)
variable.genes <- VariableFeatures(seurat.unlabeled)
a=as.data.frame(seurat.unlabeled@assays$RNA@data[variable.genes,])