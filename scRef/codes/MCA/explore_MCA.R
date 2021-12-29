library(Seurat)

file.MCA <- '/home/disk/scRef/MouseAtlas_SingleCell_Han2018/combinedMCA/MCA_combined_mouse_uniform.txt'
df.MCA <- read.table(file.MCA, header=T, row.names=1, sep='\t', check.name=F)
# a <- t(t(df.MCA)/colSums(df.MCA))*1000000
library(DESeq2)
coldata.MCA <- DataFrame(row.names = names(df.MCA))
obj.DESeq.MCA <- DESeqDataSetFromMatrix(countData = df.MCA, colData = coldata.MCA, 
                                        design = ~ 1)
fpm.MCA <- fpm(obj.DESeq.MCA, robust = T)
fpm.MCA <- as.data.frame(fpm.MCA)

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


# similarity (reference and MCA)
source('/home/zy/my_git/scRef/main/scRef.v5.R')
file.MCA <- '/home/disk/scRef/MouseAtlas_SingleCell_Han2018/combinedMCA/MCA_combined_mouse_uniform.txt'
file.MCA <- '/home/disk/scRef/MouseAtlas_SingleCell_Han2018/combinedMCA/MCA_combined_mouse.txt'
df.MCA <- read.table(file.MCA, header=T, row.names=1, sep='\t', check.name=F)

setwd('/home/zy/scRef/try_data')
# input file
file.ref <- './scRef/Reference/MouseBrain_Bulk_Zhang2014/Reference_expression.txt'
exp_ref_mat <- read.table(file.ref, header=T, row.names=1, sep='\t', check.name=F)

out=.get_cor(df.MCA, exp_ref_mat, method='kendall', CPU=4, print_step=10)

ref.names <- row.names(out)
MCA.names <- dimnames(out)[[2]]
cor.max <- c()
for (i in 1:dim(out)[1]) {
    sub.max <- max(out[i,])
    print('ref.names:')
    print(ref.names[i])
    print('MCA.names:')
    sel.names <- MCA.names[order(out[i,], decreasing = T)[1:5]]
    cor.fk(df.MCA[,sel.names])
    print(sel.names)
    cor.max <- c(cor.max, sub.max)
}

# MCA all tissue
file.MCA.outer <- '/home/disk/scRef/MouseAtlas_SingleCell_Han2018/combinedMCA/MCA_concat_outer.txt'
df.MCA.outer <- read.table(file.MCA.outer, header=T, row.names=1, sep='\t', check.name=F)

out.outer <- .get_cor(df.MCA.outer, exp_ref_mat, method='pearson', CPU=4, print_step=10)

out <- out.outer
ref.names <- row.names(out)
MCA.names <- dimnames(out)[[2]]
cor.max <- c()
for (i in 1:dim(out)[1]) {
    sub.max <- max(out[i,])
    print(ref.names[i])
    print(MCA.names[which(out[i,] == sub.max)])
    cor.max <- c(cor.max, sub.max)
}

file.MCA.inner <- '/home/disk/scRef/MouseAtlas_SingleCell_Han2018/combinedMCA/MCA_concat_inner.txt'
df.MCA.inner <- read.table(file.MCA.inner, header=T, row.names=1, sep='\t', check.name=F)

out.inner <- .get_cor(df.MCA.inner, exp_ref_mat, method='kendall', CPU=10, print_step=50)

library(pcaPP)
out <- out.inner
ref.names <- row.names(out)
MCA.names <- dimnames(out)[[2]]
cor.max <- c()
for (i in 1:dim(out)[1]) {
    sub.max <- max(out[i,])
    sel.names <- MCA.names[order(out[i,], decreasing = T)[1:10]]
    print('ref.names:')
    print(ref.names[i])
    print('MCA.names:')
    print(sel.names)
    cor.fk(df.MCA.inner[,sel.names])
    cor.max <- c(cor.max, sub.max)
}

# 
library(RUVSeq)
library(scMerge)
data("segList", package = "scMerge")
mouse_scSEG <- segList$mouse$mouse_scSEG
use.genes <- intersect(mouse_scSEG, row.names(df.MCA))
RUV.in <- newSeqExpressionSet(as.matrix(fpm.MCA))
seqRUVg <- RUVg(RUV.in, use.genes, k=3)
out.RUV <- normCounts(seqRUVg)

setwd('/home/zy/scRef/try_data')
dataset <- 'Tasic'
file.data.unlabeled <- paste0('./summary/', dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0('./summary/', dataset, '_exp_sc_mat_cluster_merged.txt')
# OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled)
# saveRDS(OUT, file = './Tasic/Tasic.Rdata')
OUT <- readRDS('./Tasic/Tasic.Rdata')
data.filter <- OUT$data.filter
label.filter <- OUT$label.filter
label.in <- data.frame(cell_id = row.names(label.filter), tag = label.filter$label.unlabeled.use.cols...)
exp.Tasic <- .generate_ref(data.filter, label.in, M='mean')
num.cells <- table(label.filter)
num.mean <- mean(num.cells)
exp.Tasic <- exp.Tasic[, names(num.cells)]
exp.Tasic.zoom <- round(exp.Tasic * num.mean)
# exp.Tasic.sum <- .generate_ref(data.filter, label.in, M='SUM')
library(DESeq2)
cell <- "Neuron"
exp.Tasic.zoom <- as.data.frame(exp.Tasic.zoom)
cells <- c(setdiff(names(exp.Tasic.zoom), cell), cell)
coldata <- DataFrame(row.names = cells)
bool.cell <- cells
bool.cell[bool.cell != cell] <- 'other'
bool.cell <- factor(bool.cell, levels = c('other', cell))
coldata$bool.cell <- bool.cell
obj.DESeq.MCA <- DESeqDataSetFromMatrix(countData = exp.Tasic.zoom[,cells], colData = coldata, 
                                        design = ~ 1 + bool.cell)
obj.DESeq.MCA$bool.cell <- relevel(obj.DESeq.MCA$bool.cell, 'other')
obj.DESeq.MCA <- DESeq(obj.DESeq.MCA)
res.DESeq <- results(obj.DESeq.MCA)
df.res <- data.frame(res.DESeq@listData, row.names = res.DESeq@rownames)

library(edgeR)
d = DGEList(counts=exp.Tasic.zoom[,cells], group=bool.cell)
d = calcNormFactors(d)

design.mat = model.matrix(~ 0 + d$samples$group)
dimnames(design.mat)[[2]] <- levels(bool.cell)
d2 = estimateGLMCommonDisp(d, design.mat)
d2 = estimateGLMTagwiseDisp(d2, design.mat)
fit = glmFit(d2, design.mat)

# 设置比较组
BvsA <- makeContrasts(
    contrasts = paste(levels(bool.cell), collapse = '-'),
    levels=design.mat)
# 组间比较,统计Fold change, Pvalue
lrt = glmLRT(fit,contrast=BvsA)
# FDR检验，控制假阳性率小于5%
de_lrt = decideTestsDGE(lrt, adjust.method="fdr", p.value=0.05)

# 导出计算结果
res.edgeR=lrt$table
res.edgeR$sig.edger=de_lrt


# overlap genes
exp_ref_mat=exp.Tasic
df.MCA=df.MCA[order(rownames(df.MCA)),]
exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
gene_MCA=rownames(df.MCA)
gene_ref=rownames(exp_ref_mat)
gene_over= gene_MCA[which(gene_MCA %in% gene_ref)]
df.MCA=df.MCA[which(gene_MCA %in% gene_over),]
exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
print('Number of overlapped genes:')
print(nrow(exp_ref_mat))

cell.MCA <- dimnames(df.MCA)[[2]]
cell.ref <- dimnames(exp_ref_mat)[[2]]
library(RUVSeq)
library(scMerge)
data("segList", package = "scMerge")
mouse_scSEG <- segList$mouse$mouse_scSEG
use.genes <- intersect(mouse_scSEG, row.names(df.MCA))
mtx.in <- cbind(df.MCA, exp_ref_mat)
names(mtx.in) <- c(paste0('MCA.', cell.MCA), paste0('Ref.', cell.ref))
seqRUVg <- RUVg(as.matrix(mtx.in), use.genes, k=10)
out.RUV <- seqRUVg$normalizedCounts


# ratio
mean.HK <- colMeans(mtx.in[use.genes,])
mtx.ratio <- t(t(mtx.in) / mean.HK)
b <- mtx.ratio[,c('MCA.Neuron', 'Ref.Neuron')]
c <- mtx.in[,c('MCA.Neuron', 'Ref.Neuron')]


# housekeeping genes
file.MCA <- '/home/disk/scRef/MouseAtlas_SingleCell_Han2018/combinedMCA/MCA_combined_mouse_uniform.txt'
df.MCA <- read.table(file.MCA, header=T, row.names=1, sep='\t', check.name=F)
# a <- t(t(df.MCA)/colSums(df.MCA))*1000000
library(DESeq2)
coldata.MCA <- DataFrame(row.names = names(df.MCA))
obj.DESeq.MCA <- DESeqDataSetFromMatrix(countData = df.MCA, colData = coldata.MCA, 
                                        design = ~ 1)
fpm.MCA <- fpm(obj.DESeq.MCA, robust = T)
fpm.MCA <- as.data.frame(fpm.MCA)
log2.MCA <- log2(fpm.MCA + 1)

vec.0 <- apply(as.matrix(log2.MCA), 1, function(x) {sum(x == 0)})
genes.0 <- names(vec.0[vec.0 == 0])
vec.sd <- apply(as.matrix(log2.MCA), 1, sd)
genes.sd <- names(vec.sd[vec.sd < 1])
HK.genes <- intersect(genes.0, genes.sd)
mean.cell <- apply(as.matrix(log2.MCA), 2, mean)
vec.fc <- 

a <- t(log2.MCA[HK.genes,])

setwd('/home/zy/scRef/try_data')
# input file
file.ref <- './scRef/Reference/MouseBrain_Bulk_Zhang2014/Reference_expression.txt'
exp_ref_mat <- read.table(file.ref, header=T, row.names=1, sep='\t', check.name=F)
exp_ref_mat <- log1p(exp_ref_mat)
b <- t(exp_ref_mat[HK.genes,])
rowMeans2(b, na.rm = T)

library(Seurat)
seurat.MCA <- CreateSeuratObject(counts = df.MCA, project = "MCA", min.cells = 0, min.features = 1000)
seurat.MCA <- NormalizeData(seurat.MCA, normalization.method = "LogNormalize", scale.factor = 1e6)
mtx.norm <- as.matrix(seurat.MCA@assays$RNA@data)
a <- t(mtx.norm[HK.genes,])
rowMeans2(a, na.rm = T)

library(RUVSeq)
seqRUVg <- RUVg(as.matrix(mtx.in), use.genes, k=5, isLog = T)
out.RUV <- seqRUVg$normalizedCounts


source('/home/zy/my_git/scRef/main/scRef.v5.R')
# two steps
setwd('/home/zy/scRef/try_data')
# input file
file.ref <- './scRef/Reference/MouseBrain_Bulk_Zhang2014/Reference_expression.txt'
exp_ref_mat <- read.table(file.ref, header=T, row.names=1, sep='\t', check.name=F)
exp_ref_mat <- log1p(exp_ref_mat)

file.MCA <- '/home/disk/scRef/MouseAtlas_SingleCell_Han2018/combinedMCA/MCA_combined_mouse_uniform.txt'
df.MCA <- read.table(file.MCA, header=T, row.names=1, sep='\t', check.name=F)
# MCA all tissue
file.MCA.outer <- '/home/disk/scRef/MouseAtlas_SingleCell_Han2018/combinedMCA/MCA_concat_outer.txt'
df.MCA.outer <- read.table(file.MCA.outer, header=T, row.names=1, sep='\t', check.name=F)
df.MCA <- df.MCA.outer

# Tasic
OUT <- readRDS('./Tasic/Tasic.Rdata')
data.filter <- OUT$data.filter
label.filter <- OUT$label.filter
label.in <- data.frame(cell_id = row.names(label.filter), tag = label.filter$label.unlabeled.use.cols...)
exp.Tasic.sum <- .generate_ref(data.filter, label.in, M='SUM')
seurat.Tasic <- CreateSeuratObject(counts = exp.Tasic.sum, project = "MCA", min.cells = 0, min.features = 1000)
seurat.Tasic <- NormalizeData(seurat.Tasic, normalization.method = "LogNormalize", scale.factor = 1e6)
# exp_ref_mat <- as.matrix(seurat.Tasic@assays$RNA@data)
exp.Tasic <- as.matrix(seurat.Tasic@assays$RNA@data)

# overlap genes
exp_ref_mat=exp_ref_mat
df.MCA=df.MCA[order(rownames(df.MCA)),]
exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
gene_MCA=rownames(df.MCA)
gene_ref=rownames(exp_ref_mat)
gene_over= gene_MCA[which(gene_MCA %in% gene_ref)]
df.MCA=df.MCA[which(gene_MCA %in% gene_over),]
exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
print('Number of overlapped genes:')
print(nrow(exp_ref_mat))

library(Seurat)
seurat.MCA <- CreateSeuratObject(counts = df.MCA, project = "MCA", min.cells = 0, min.features = 1000)
seurat.MCA <- NormalizeData(seurat.MCA, normalization.method = "LogNormalize", scale.factor = 1e6)
mtx.norm <- as.matrix(seurat.MCA@assays$RNA@data)


cell <- "Astrocytes"
cells <- c(setdiff(names(exp_ref_mat), cell), cell)
bool.cell <- cells
bool.cell[bool.cell != cell] <- '1'
bool.cell[bool.cell == cell] <- '2'
res.limma <- .getDEgeneF(exp_ref_mat[,cells], bool.cell)
vec.cell <- exp_ref_mat[, cell]
exp.top10 <- quantile(vec.cell, 0.9)
genes.high <- row.names(exp_ref_mat)[vec.cell > exp.top10]
res.limma.high <- res.limma[genes.high,]
res.limma.high <- res.limma.high[order(res.limma.high$adj.P.Val),]

c <- as.data.frame(mtx.norm['Reln',])
top100 <- row.names(res.limma.high)[1:100]
d <- t(mtx.norm[top100,])
b <- t(exp_ref_mat['Hgf',])

out1=.get_cor(mtx.norm, exp_ref_mat, method='pearson', CPU=4, print_step=10)
out1 <- t(out1)
tag1=.get_tag_max(out1)

