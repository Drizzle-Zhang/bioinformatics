library(Seurat)
library(dplyr)
library(ggplot2)

## 读取数据
# 输入路径
path.data <- '/local/hanyang/data/single-cell/'
samples <- list.files(path.data, all.files = T, no.. = T)

list.seurat <- c()
for (sample in samples) {
    sub.mat <- Read10X(paste0(path.data, sample))
    dimnames(sub.mat)[[2]] <- paste(sample, dimnames(sub.mat)[[2]], sep = '_')
    sub.seurat <- CreateSeuratObject(counts = sub.mat)
    sub.seurat@meta.data$sample <- rep(sample, dim(sub.mat)[2])
    sub.seurat[["percent.mt"]] <- PercentageFeatureSet(sub.seurat, pattern = "^MT-")
    if (sample == samples[1]) {
        seurat.C1 <- sub.seurat
    } else {
        list.seurat <- c(list.seurat, sub.seurat)
    }
}
seurat.merge <- merge(seurat.C1, list.seurat)

####################################################
# cutoff median-a*sd ~ median+b*sd
# default a=1 b=2
# nFeature cut 

a=1
b=2
nF_cut <- c()
for (obj in list.seurat) {
    nF_cut <- c(nF_cut, median(obj$nFeature_RNA) + b*sd(obj$nFeature_RNA))
    nF_cut <- c(nF_cut, median(obj$nFeature_RNA) - a*sd(obj$nFeature_RNA))
}
names(nF_cut) <- 
    cutoff_line.df <- 
    data.frame(
        x=rep(c(1:length(list.seurat)-0.3), each=2), 
        xend=rep(c(1:length(list.seurat)+0.3), each=2), 
        y=nF_cut, 
        yend=nF_cut)
p_nf <- VlnPlot(QC.all, features = "nFeature_RNA")
for (i in 1:nrow(cutoff_line.df)) {
    p_nf <- p_nf + geom_segment(x=cutoff_line.df$x[i], 
                            xend=cutoff_line.df$xend[i], 
                            y=cutoff_line.df$y[i], 
                            yend=cutoff_line.df$yend[i], 
                            colour="red")
}
p_nf

######################################
# mt cut 
# cutoff  ~ median+c*sd
# default c=2

c=2
mt_cut <- c()
for (obj in list.seurat) {
    mt_cut <- c(mt_cut, median(obj$percent.mt) + c*sd(obj$percent.mt))
}
cutoff_line.df <- 
    data.frame(
        x=rep(c(1:length(list.seurat)-0.3)), 
        xend=rep(c(1:length(list.seurat)+0.3)), 
        y=mt_cut, 
        yend=mt_cut)
p_mt <- VlnPlot(QC.all, features = "percent.mt")
for (i in 1:nrow(cutoff_line.df)) {
    p_mt <- p_mt + geom_segment(x=cutoff_line.df$x[i], 
                            xend=cutoff_line.df$xend[i], 
                            y=cutoff_line.df$y[i], 
                            yend=cutoff_line.df$yend[i], 
                            colour="red")
}
p_mt

for (i in 1:length(list.seurat)) {
    list.seurat[[i]] <- subset(list.seurat[[i]], 
                               subset = nFeature_RNA < nF_cut[2*i-1] &
                                   nFeature_RNA > nF_cut[2*i] &
                                   percent.mt < mt_cut[i])
}

# 提供质控标准
feature.min <- 800
feature.max <- 6000
UMI.min <- 2000
UMI.max <- 50000
percent.mt.max <- 15
seurat.all_filter <- subset(seurat.merge, subset = nFeature_RNA > feature.min & nFeature_RNA < feature.max & nCount_RNA > UMI.min & nCount_RNA < UMI.max & percent.mt < percent.mt.max)

# 质控结果的展示
# nFeature_RNA
h_df <- data.frame(seurat.merge$sample, seurat.merge$nFeature_RNA)
colnames(h_df) <- c("Sample", "nFeature_RNA")
VlnPlot(seurat.merge, features = "nFeature_RNA", pt.size = 0.1, group.by = 'sample') + 
    geom_hline(aes(yintercept=feature.min), colour="blue", linetype="dashed") +
    geom_hline(aes(yintercept=feature.max), colour="blue", linetype="dashed")
ggplot(h_df, aes(x=nFeature_RNA, color=Sample)) + 
    geom_histogram(binwidth = 20, fill = "white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    geom_vline(aes(xintercept=feature.min), colour="blue", linetype="dashed") + 
    geom_vline(aes(xintercept=feature.max), colour="blue", linetype="dashed") + 
    ylab("Cell Number")
# nCount_RNA
h2_df <- data.frame(seurat.merge$sample, seurat.merge$nCount_RNA)
colnames(h2_df) <- c("Sample", "nCount_RNA")
VlnPlot(seurat.merge, features = "nCount_RNA", pt.size = 0.001) + 
    geom_hline(aes(yintercept=UMI.min), colour="blue", linetype="dashed") +
    geom_hline(aes(yintercept=UMI.max), colour="blue", linetype="dashed")
ggplot(h2_df.M3, aes(x=nCount_RNA, color=Sample)) + 
    geom_histogram(binwidth = 100, fill="white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    geom_vline(aes(xintercept=UMI.min), colour="black", linetype="dashed") + 
    geom_vline(aes(xintercept=UMI.max), colour="black", linetype="dashed") + 
    ylab("Cell Number")
# percent.mt
h3_df.M3 <- data.frame(seurat.merge$sample, seurat.merge$percent.mt)
colnames(h3_df.M3) <- c("Sample", "percent.mt")
VlnPlot(seurat.merge, features = "percent.mt", pt.size = 0.01) + 
    geom_hline(aes(yintercept=percent.mt.max), colour="blue", linetype="dashed")
ggplot(h3_df.M3, aes(x=percent.mt, color=Sample)) + 
    geom_histogram(binwidth = 1, fill="white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    ylab("Cell Number") + 
    geom_vline(aes(xintercept=percent.mt.max), colour="blue", linetype="dashed")
# feature scatter
FeatureScatter(seurat.merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
    geom_vline(aes(xintercept=4000), colour="blue", linetype="dashed") + 
    geom_vline(aes(xintercept=35000), colour="blue", linetype="dashed") + 
    geom_hline(aes(yintercept=1000), colour="blue", linetype="dashed") + 
    geom_hline(aes(yintercept=5000), colour="blue", linetype="dashed")


######
#注：上述质控过程是把所有样本合并后质控的，如果样本间差异过大，需要分别进行质控
#######

### 数据预处理
### 参数列表
num.seed <- 123
num.features <- 4000
num_PCs <- 50
n.neighbors <- 30
min.dist <- 0.3
### 这里也可能需要输入样本信息，比如细胞来自于正常组织还是疾病组织等

Obj.seurat.input <- seurat.all_filter
Obj.seurat.input <- NormalizeData(object = Obj.seurat.input)
Obj.seurat.input <- FindVariableFeatures(object = Obj.seurat.input, 
                                         selection.method = "vst", nfeatures = num.features)
Obj.seurat.input <- ScaleData(object = Obj.seurat.input, 
                              features = VariableFeatures(object = Obj.seurat.input), 
                              split.by = "sample")
Obj.seurat.input <- RunPCA(object = Obj.seurat.input, seed.use = num.seed, npcs=num_PCs, verbose = F)
Obj.seurat.input <- RunTSNE(Obj.seurat.input, dims = 1:num_PCs, 
                            seed.use = num.seed, n.components=2)
Obj.seurat.input <- RunUMAP(Obj.seurat.input, dims = 1:num_PCs, 
                            seed.use = num.seed, n.components=2,
                            n.neighbors = n.neighbors, min.dist = min.dist)

# 输出降维之后的图
DimPlot(Obj.seurat.input, group.by = "sample", reduction = 'tsne')
DimPlot(Obj.seurat.input, group.by = "sample", reduction = 'umap')

# 聚类
# 参数
use.reduction <- "pca"
use.dims <- 50
cluster.k <- 20
cluster.resolution <- 2
Obj.seurat.input <- FindNeighbors(Obj.seurat.input, reduction = use.reduction, 
                                  dims = 1:use.dims, k.param = cluster.k)
Obj.seurat.input <- FindClusters(Obj.seurat.input, resolution = cluster.resolution)

# 展示聚类结果
DimPlot(Obj.seurat.input, group.by = "seurat_clusters", reduction = 'umap')


args <- commandArgs(T)
correct.dhs.score(args[1], args[2])
