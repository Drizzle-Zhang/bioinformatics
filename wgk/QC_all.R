library(Seurat)
library(ggplot2)

path.data <- '/home/yzj/JingMA_ORI/data/'

# read data
C1 <- Read10X(paste0(path.data, 'C1'))
dimnames(C1)[[2]] <- paste('C1', dimnames(C1)[[2]], sep = '_')
seurat.C1 <- CreateSeuratObject(counts = C1)
seurat.C1@meta.data$sample <- rep('C1', dim(C1)[2])
seurat.C1[["percent.mt"]] <- PercentageFeatureSet(seurat.C1, pattern = "^MT-")

C2 <- Read10X(paste0(path.data, 'C2'))
dimnames(C2)[[2]] <- paste('C2', dimnames(C2)[[2]], sep = '_')
seurat.C2 <- CreateSeuratObject(counts = C2)
seurat.C2@meta.data$sample <- rep('C2', dim(C2)[2])
seurat.C2[["percent.mt"]] <- PercentageFeatureSet(seurat.C2, pattern = "^MT-")

C3 <- Read10X(paste0(path.data, 'C3'))
dimnames(C3)[[2]] <- paste('C3', dimnames(C3)[[2]], sep = '_')
seurat.C3 <- CreateSeuratObject(counts = C3)
seurat.C3@meta.data$sample <- rep('C3', dim(C3)[2])
seurat.C3[["percent.mt"]] <- PercentageFeatureSet(seurat.C3, pattern = "^MT-")

C4 <- Read10X(paste0(path.data, 'C4'))
dimnames(C4)[[2]] <- paste('C4', dimnames(C4)[[2]], sep = '_')
seurat.C4 <- CreateSeuratObject(counts = C4)
seurat.C4@meta.data$sample <- rep('C4', dim(C4)[2])
seurat.C4[["percent.mt"]] <- PercentageFeatureSet(seurat.C4, pattern = "^MT-")

C5 <- Read10X(paste0(path.data, 'C5'))
dimnames(C5)[[2]] <- paste('C5', dimnames(C5)[[2]], sep = '_')
seurat.C5 <- CreateSeuratObject(counts = C5)
seurat.C5@meta.data$sample <- rep('C5', dim(C5)[2])
seurat.C5[["percent.mt"]] <- PercentageFeatureSet(seurat.C5, pattern = "^MT-")

C6 <- Read10X(paste0(path.data, 'C6'))
dimnames(C6)[[2]] <- paste('C6', dimnames(C6)[[2]], sep = '_')
seurat.C6 <- CreateSeuratObject(counts = C6)
seurat.C6@meta.data$sample <- rep('C6', dim(C6)[2])
seurat.C6[["percent.mt"]] <- PercentageFeatureSet(seurat.C6, pattern = "^MT-")

M1 <- Read10X(paste0(path.data, 'M1'))
dimnames(M1)[[2]] <- paste('M1', dimnames(M1)[[2]], sep = '_')
seurat.M1 <- CreateSeuratObject(counts = M1)
seurat.M1@meta.data$sample <- rep('M1', dim(M1)[2])
seurat.M1[["percent.mt"]] <- PercentageFeatureSet(seurat.M1, pattern = "^MT-")

M2 <- Read10X(paste0(path.data, 'M2'))
dimnames(M2)[[2]] <- paste('M2', dimnames(M2)[[2]], sep = '_')
seurat.M2 <- CreateSeuratObject(counts = M2)
seurat.M2@meta.data$sample <- rep('M2', dim(M2)[2])
seurat.M2[["percent.mt"]] <- PercentageFeatureSet(seurat.M2, pattern = "^MT-")

M3 <- Read10X(paste0(path.data, 'M3'))
dimnames(M3)[[2]] <- paste('M3', dimnames(M3)[[2]], sep = '_')
seurat.M3 <- CreateSeuratObject(counts = M3)
seurat.M3@meta.data$sample <- rep('M3', dim(M3)[2])
seurat.M3[["percent.mt"]] <- PercentageFeatureSet(seurat.M3, pattern = "^MT-")

# merge
seurat.merge <- merge(seurat.C1, 
                      c(seurat.C2, seurat.C3, seurat.C4, seurat.C5, seurat.C6,
                        seurat.M1, seurat.M2, seurat.M3))

# QC
# nFeature_RNA
h_df <- data.frame(seurat.merge$sample, seurat.merge$nFeature_RNA)
colnames(h_df) <- c("Sample", "nFeature_RNA")
VlnPlot(seurat.merge, features = "nFeature_RNA", pt.size = 0.1, group.by = 'sample') + 
    geom_hline(aes(yintercept=1000), colour="blue", linetype="dashed") +
    geom_hline(aes(yintercept=5000), colour="blue", linetype="dashed")
ggplot(h_df, aes(x=nFeature_RNA, color=Sample)) + 
    geom_histogram(binwidth = 20, fill = "white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    geom_vline(aes(xintercept=1000), colour="blue", linetype="dashed") + 
    geom_vline(aes(xintercept=5500), colour="blue", linetype="dashed") + 
    ylab("Cell Number")
# nCount_RNA
h2_df <- data.frame(seurat.merge$sample, seurat.merge$nCount_RNA)
colnames(h2_df) <- c("Sample", "nCount_RNA")
VlnPlot(seurat.merge, features = "nCount_RNA", pt.size = 0.001) + 
    geom_hline(aes(yintercept=4000), colour="blue", linetype="dashed") +
    geom_hline(aes(yintercept=35000), colour="blue", linetype="dashed")
ggplot(h2_df.M3, aes(x=nCount_RNA, color=Sample)) + 
    geom_histogram(binwidth = 100, fill="white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    geom_vline(aes(xintercept=2500), colour="black", linetype="dashed") + 
    geom_vline(aes(xintercept=30000), colour="black", linetype="dashed") + 
    ylab("Cell Number")
# percent.mt
h3_df.M3 <- data.frame(seurat.merge$sample, seurat.merge$percent.mt)
colnames(h3_df.M3) <- c("Sample", "percent.mt")
VlnPlot(seurat.merge, features = "percent.mt", pt.size = 0.01) + 
    geom_hline(aes(yintercept=10), colour="blue", linetype="dashed")
ggplot(h3_df.M3, aes(x=percent.mt, color=Sample)) + 
    geom_histogram(binwidth = 1, fill="white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    ylab("Cell Number") + geom_vline(aes(xintercept=9), colour="blue", linetype="dashed")
# feature scatter
FeatureScatter(seurat.merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
    geom_vline(aes(xintercept=4000), colour="blue", linetype="dashed") + 
    geom_vline(aes(xintercept=35000), colour="blue", linetype="dashed") + 
    geom_hline(aes(yintercept=1000), colour="blue", linetype="dashed") + 
    geom_hline(aes(yintercept=5000), colour="blue", linetype="dashed")
# filter  
seurat.all_filter <- subset(seurat.merge, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & nCount_RNA > 4000 & nCount_RNA < 35000 & percent.mt < 10)

# status
samples <- seurat.all_filter$sample
vec.st <- rep('0', length(samples))
vec.st[samples %in% c('C1', 'C2', 'C3', 'C4', 'C5', 'C6')] <- 'Normal'
vec.st[samples %in% c('M1', 'M2', 'M3')] <- 'Microtia'
seurat.all_filter$status <- vec.st

file.ear <- '/home/disk/drizzle/wgk/data/seurat_all_filtered.Rdata'
saveRDS(seurat.all_filter, file = file.ear)

# dim(seurat.all_filter)
# 32738 46916




