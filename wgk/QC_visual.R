library(Seurat)
library(ggplot2)

path.data <- '/home/yzj/JingMA/data/'

############ C1
C1 <- Read10X(paste0(path.data, 'C1'))
dimnames(C1)[[2]] <- paste('C1', dimnames(C1)[[2]], sep = '_')
seurat.C1 <- CreateSeuratObject(counts = C1)
seurat.C1@meta.data$sample <- rep('C1', dim(C1)[2])
seurat.C1[["percent.mt"]] <- PercentageFeatureSet(seurat.C1, pattern = "^MT-")
# nFeature_RNA
h_df.C1 <- data.frame(seurat.C1$sample, seurat.C1$nFeature_RNA)
colnames(h_df.C1) <- c("Sample", "nFeature_RNA")
VlnPlot(seurat.C1, features = "nFeature_RNA", pt.size = 0.1) + 
    geom_hline(aes(yintercept=1500), colour="blue", linetype="dashed") +
    geom_hline(aes(yintercept=5500), colour="blue", linetype="dashed")
ggplot(h_df.C1, aes(x=nFeature_RNA, color=Sample)) + 
    geom_histogram(binwidth = 50, fill = "white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    geom_vline(aes(xintercept=1500), colour="blue", linetype="dashed") + 
    geom_vline(aes(xintercept=5500), colour="blue", linetype="dashed") + 
    ylab("Cell Number")
# nCount_RNA
h2_df.C1 <- data.frame(seurat.C1$sample, seurat.C1$nCount_RNA)
colnames(h2_df.C1) <- c("Sample", "nCount_RNA")
VlnPlot(seurat.C1, features = "nCount_RNA", pt.size = 0.1) + 
    geom_hline(aes(yintercept=4000), colour="blue", linetype="dashed") +
    geom_hline(aes(yintercept=37000), colour="blue", linetype="dashed")
ggplot(h2_df.C1, aes(x=nCount_RNA, color=Sample)) + 
    geom_histogram(binwidth = 100, fill="white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    geom_vline(aes(xintercept=4000), colour="blue", linetype="dashed") + 
    geom_vline(aes(xintercept=37000), colour="blue", linetype="dashed") + 
    ylab("Cell Number")
# percent.mt
h3_df.C1 <- data.frame(seurat.C1$sample, seurat.C1$percent.mt)
colnames(h3_df.C1) <- c("Sample", "percent.mt")
VlnPlot(seurat.C1, features = "percent.mt", pt.size = 0.1) + 
    geom_hline(aes(yintercept=14), colour="blue", linetype="dashed")
ggplot(h3_df.C1, aes(x=percent.mt, color=Sample)) + 
    geom_histogram(binwidth = 2, fill="white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    ylab("Cell Number") + geom_vline(aes(xintercept=14), colour="red", linetype="dashed")
# feature scatter
FeatureScatter(seurat.C1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
    geom_vline(aes(xintercept=4000), colour="blue", linetype="dashed") + 
    geom_vline(aes(xintercept=37000), colour="blue", linetype="dashed") + 
    geom_hline(aes(yintercept=1500), colour="blue", linetype="dashed") + 
    geom_hline(aes(yintercept=5500), colour="blue", linetype="dashed")
# filter  
seurat.C1_filter <- subset(seurat.C1, subset = nFeature_RNA > 1500 & nFeature_RNA < 5500 & nCount_RNA > 4000 & nCount_RNA < 37000 & percent.mt < 14)


############ C2
C2 <- Read10X(paste0(path.data, 'C2'))
dimnames(C2)[[2]] <- paste('C2', dimnames(C2)[[2]], sep = '_')
seurat.C2 <- CreateSeuratObject(counts = C2)
seurat.C2@meta.data$sample <- rep('C2', dim(C2)[2])
seurat.C2[["percent.mt"]] <- PercentageFeatureSet(seurat.C2, pattern = "^MT-")
# nFeature_RNA
h_df.C2 <- data.frame(seurat.C2$sample, seurat.C2$nFeature_RNA)
colnames(h_df.C2) <- c("Sample", "nFeature_RNA")
VlnPlot(seurat.C2, features = "nFeature_RNA", pt.size = 0.1) + 
    geom_hline(aes(yintercept=1200), colour="blue", linetype="dashed") +
    geom_hline(aes(yintercept=5500), colour="blue", linetype="dashed")
ggplot(h_df.C2, aes(x=nFeature_RNA, color=Sample)) + 
    geom_histogram(binwidth = 50, fill = "white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    geom_vline(aes(xintercept=1200), colour="blue", linetype="dashed") + 
    geom_vline(aes(xintercept=5500), colour="blue", linetype="dashed") + 
    ylab("Cell Number")
# nCount_RNA
h2_df.C2 <- data.frame(seurat.C2$sample, seurat.C2$nCount_RNA)
colnames(h2_df.C2) <- c("Sample", "nCount_RNA")
VlnPlot(seurat.C2, features = "nCount_RNA", pt.size = 0.1) + 
    geom_hline(aes(yintercept=4200), colour="blue", linetype="dashed") +
    geom_hline(aes(yintercept=40000), colour="blue", linetype="dashed")
ggplot(h2_df.C2, aes(x=nCount_RNA, color=Sample)) + 
    geom_histogram(binwidth = 100, fill="white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    geom_vline(aes(xintercept=4200), colour="black", linetype="dashed") + 
    geom_vline(aes(xintercept=40000), colour="black", linetype="dashed") + 
    ylab("Cell Number")
# percent.mt
h3_df.C2 <- data.frame(seurat.C2$sample, seurat.C2$percent.mt)
colnames(h3_df.C2) <- c("Sample", "percent.mt")
VlnPlot(seurat.C2, features = "percent.mt", pt.size = 0.1) + 
    geom_hline(aes(yintercept=13), colour="blue", linetype="dashed")
ggplot(h3_df.C2, aes(x=percent.mt, color=Sample)) + 
    geom_histogram(binwidth = 2, fill="white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    ylab("Cell Number") + geom_vline(aes(xintercept=13), colour="red", linetype="dashed")
# feature scatter
FeatureScatter(seurat.C2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
    geom_vline(aes(xintercept=4200), colour="blue", linetype="dashed") + 
    geom_vline(aes(xintercept=40000), colour="blue", linetype="dashed") + 
    geom_hline(aes(yintercept=1200), colour="blue", linetype="dashed") + 
    geom_hline(aes(yintercept=5500), colour="blue", linetype="dashed")
# filter  
seurat.C2_filter <- subset(seurat.C2, subset = nFeature_RNA > 1200 & nFeature_RNA < 5500 & nCount_RNA > 4200 & nCount_RNA < 40000 & percent.mt < 13)


############ C3
C3 <- Read10X(paste0(path.data, 'C3'))
dimnames(C3)[[2]] <- paste('C3', dimnames(C3)[[2]], sep = '_')
seurat.C3 <- CreateSeuratObject(counts = C3)
seurat.C3@meta.data$sample <- rep('C3', dim(C3)[2])
seurat.C3[["percent.mt"]] <- PercentageFeatureSet(seurat.C3, pattern = "^MT-")
# nFeature_RNA
h_df.C3 <- data.frame(seurat.C3$sample, seurat.C3$nFeature_RNA)
colnames(h_df.C3) <- c("Sample", "nFeature_RNA")
VlnPlot(seurat.C3, features = "nFeature_RNA", pt.size = 0.1) + 
    geom_hline(aes(yintercept=1600), colour="blue", linetype="dashed") +
    geom_hline(aes(yintercept=5500), colour="blue", linetype="dashed")
ggplot(h_df.C3, aes(x=nFeature_RNA, color=Sample)) + 
    geom_histogram(binwidth = 50, fill = "white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    geom_vline(aes(xintercept=1600), colour="blue", linetype="dashed") + 
    geom_vline(aes(xintercept=5500), colour="blue", linetype="dashed") + 
    ylab("Cell Number")
# nCount_RNA
h2_df.C3 <- data.frame(seurat.C3$sample, seurat.C3$nCount_RNA)
colnames(h2_df.C3) <- c("Sample", "nCount_RNA")
VlnPlot(seurat.C3, features = "nCount_RNA", pt.size = 0.1) + 
    geom_hline(aes(yintercept=5000), colour="blue", linetype="dashed") +
    geom_hline(aes(yintercept=40000), colour="blue", linetype="dashed")
ggplot(h2_df.C3, aes(x=nCount_RNA, color=Sample)) + 
    geom_histogram(binwidth = 100, fill="white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    geom_vline(aes(xintercept=5000), colour="black", linetype="dashed") + 
    geom_vline(aes(xintercept=40000), colour="black", linetype="dashed") + 
    ylab("Cell Number")
# percent.mt
h3_df.C3 <- data.frame(seurat.C3$sample, seurat.C3$percent.mt)
colnames(h3_df.C3) <- c("Sample", "percent.mt")
VlnPlot(seurat.C3, features = "percent.mt", pt.size = 0.1) + 
    geom_hline(aes(yintercept=13), colour="blue", linetype="dashed")
ggplot(h3_df.C3, aes(x=percent.mt, color=Sample)) + 
    geom_histogram(binwidth = 1, fill="white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    ylab("Cell Number") + geom_vline(aes(xintercept=13), colour="blue", linetype="dashed")
# feature scatter
FeatureScatter(seurat.C3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
    geom_vline(aes(xintercept=5000), colour="blue", linetype="dashed") + 
    geom_vline(aes(xintercept=40000), colour="blue", linetype="dashed") + 
    geom_hline(aes(yintercept=1600), colour="blue", linetype="dashed") + 
    geom_hline(aes(yintercept=5500), colour="blue", linetype="dashed")
# filter  
seurat.C3_filter <- subset(seurat.C3, subset = nFeature_RNA > 1600 & nFeature_RNA < 5500 & nCount_RNA > 5000 & nCount_RNA < 40000 & percent.mt < 13)


############ C4
C4 <- Read10X(paste0(path.data, 'C4'))
dimnames(C4)[[2]] <- paste('C4', dimnames(C4)[[2]], sep = '_')
seurat.C4 <- CreateSeuratObject(counts = C4)
seurat.C4@meta.data$sample <- rep('C4', dim(C4)[2])
seurat.C4[["percent.mt"]] <- PercentageFeatureSet(seurat.C4, pattern = "^MT-")
# nFeature_RNA
h_df.C4 <- data.frame(seurat.C4$sample, seurat.C4$nFeature_RNA)
colnames(h_df.C4) <- c("Sample", "nFeature_RNA")
VlnPlot(seurat.C4, features = "nFeature_RNA", pt.size = 0.1) + 
    geom_hline(aes(yintercept=2000), colour="blue", linetype="dashed") +
    geom_hline(aes(yintercept=5000), colour="blue", linetype="dashed")
ggplot(h_df.C4, aes(x=nFeature_RNA, color=Sample)) + 
    geom_histogram(binwidth = 20, fill = "white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    geom_vline(aes(xintercept=2000), colour="blue", linetype="dashed") + 
    geom_vline(aes(xintercept=5000), colour="blue", linetype="dashed") + 
    ylab("Cell Number")
# nCount_RNA
h2_df.C4 <- data.frame(seurat.C4$sample, seurat.C4$nCount_RNA)
colnames(h2_df.C4) <- c("Sample", "nCount_RNA")
VlnPlot(seurat.C4, features = "nCount_RNA", pt.size = 0.1) + 
    geom_hline(aes(yintercept=6000), colour="blue", linetype="dashed") +
    geom_hline(aes(yintercept=50000), colour="blue", linetype="dashed")
ggplot(h2_df.C4, aes(x=nCount_RNA, color=Sample)) + 
    geom_histogram(binwidth = 100, fill="white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    geom_vline(aes(xintercept=6000), colour="black", linetype="dashed") + 
    geom_vline(aes(xintercept=50000), colour="black", linetype="dashed") + 
    ylab("Cell Number")
# percent.mt
h3_df.C4 <- data.frame(seurat.C4$sample, seurat.C4$percent.mt)
colnames(h3_df.C4) <- c("Sample", "percent.mt")
VlnPlot(seurat.C4, features = "percent.mt", pt.size = 0.1) + 
    geom_hline(aes(yintercept=5), colour="blue", linetype="dashed")
ggplot(h3_df.C4, aes(x=percent.mt, color=Sample)) + 
    geom_histogram(binwidth = 1, fill="white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    ylab("Cell Number") + geom_vline(aes(xintercept=5), colour="blue", linetype="dashed")
# feature scatter
FeatureScatter(seurat.C4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
    geom_vline(aes(xintercept=6000), colour="blue", linetype="dashed") + 
    geom_vline(aes(xintercept=50000), colour="blue", linetype="dashed") + 
    geom_hline(aes(yintercept=2000), colour="blue", linetype="dashed") + 
    geom_hline(aes(yintercept=5000), colour="blue", linetype="dashed")
# filter  
seurat.C4_filter <- subset(seurat.C4, subset = nFeature_RNA > 2000 & nFeature_RNA < 5000 & nCount_RNA > 6000 & nCount_RNA < 50000 & percent.mt < 5)


############ M1
M1 <- Read10X(paste0(path.data, 'M1'))
dimnames(M1)[[2]] <- paste('M1', dimnames(M1)[[2]], sep = '_')
seurat.M1 <- CreateSeuratObject(counts = M1)
seurat.M1@meta.data$sample <- rep('M1', dim(M1)[2])
seurat.M1[["percent.mt"]] <- PercentageFeatureSet(seurat.M1, pattern = "^MT-")
# nFeature_RNA
h_df.M1 <- data.frame(seurat.M1$sample, seurat.M1$nFeature_RNA)
colnames(h_df.M1) <- c("Sample", "nFeature_RNA")
VlnPlot(seurat.M1, features = "nFeature_RNA", pt.size = 0.1) + 
    geom_hline(aes(yintercept=800), colour="blue", linetype="dashed") +
    geom_hline(aes(yintercept=5500), colour="blue", linetype="dashed")
ggplot(h_df.M1, aes(x=nFeature_RNA, color=Sample)) + 
    geom_histogram(binwidth = 20, fill = "white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    geom_vline(aes(xintercept=800), colour="blue", linetype="dashed") + 
    geom_vline(aes(xintercept=5500), colour="blue", linetype="dashed") + 
    ylab("Cell Number")
# nCount_RNA
h2_df.M1 <- data.frame(seurat.M1$sample, seurat.M1$nCount_RNA)
colnames(h2_df.M1) <- c("Sample", "nCount_RNA")
VlnPlot(seurat.M1, features = "nCount_RNA", pt.size = 0.1) + 
    geom_hline(aes(yintercept=2500), colour="blue", linetype="dashed") +
    geom_hline(aes(yintercept=50000), colour="blue", linetype="dashed")
ggplot(h2_df.M1, aes(x=nCount_RNA, color=Sample)) + 
    geom_histogram(binwidth = 100, fill="white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    geom_vline(aes(xintercept=2500), colour="black", linetype="dashed") + 
    geom_vline(aes(xintercept=50000), colour="black", linetype="dashed") + 
    ylab("Cell Number")
# percent.mt
h3_df.M1 <- data.frame(seurat.M1$sample, seurat.M1$percent.mt)
colnames(h3_df.M1) <- c("Sample", "percent.mt")
VlnPlot(seurat.M1, features = "percent.mt", pt.size = 0.1) + 
    geom_hline(aes(yintercept=13), colour="blue", linetype="dashed")
ggplot(h3_df.M1, aes(x=percent.mt, color=Sample)) + 
    geom_histogram(binwidth = 1, fill="white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    ylab("Cell Number") + geom_vline(aes(xintercept=13), colour="blue", linetype="dashed")
# feature scatter
FeatureScatter(seurat.M1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
    geom_vline(aes(xintercept=2500), colour="blue", linetype="dashed") + 
    geom_vline(aes(xintercept=50000), colour="blue", linetype="dashed") + 
    geom_hline(aes(yintercept=800), colour="blue", linetype="dashed") + 
    geom_hline(aes(yintercept=5500), colour="blue", linetype="dashed")
# filter  
seurat.M1_filter <- subset(seurat.M1, subset = nFeature_RNA > 800 & nFeature_RNA < 5500 & nCount_RNA > 2500 & nCount_RNA < 50000 & percent.mt < 13)


############ M2
M2 <- Read10X(paste0(path.data, 'M2'))
dimnames(M2)[[2]] <- paste('M2', dimnames(M2)[[2]], sep = '_')
seurat.M2 <- CreateSeuratObject(counts = M2)
seurat.M2@meta.data$sample <- rep('M2', dim(M2)[2])
seurat.M2[["percent.mt"]] <- PercentageFeatureSet(seurat.M2, pattern = "^MT-")
# nFeature_RNA
h_df.M2 <- data.frame(seurat.M2$sample, seurat.M2$nFeature_RNA)
colnames(h_df.M2) <- c("Sample", "nFeature_RNA")
VlnPlot(seurat.M2, features = "nFeature_RNA", pt.size = 0.1) + 
    geom_hline(aes(yintercept=800), colour="blue", linetype="dashed") +
    geom_hline(aes(yintercept=5500), colour="blue", linetype="dashed")
ggplot(h_df.M2, aes(x=nFeature_RNA, color=Sample)) + 
    geom_histogram(binwidth = 20, fill = "white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    geom_vline(aes(xintercept=800), colour="blue", linetype="dashed") + 
    geom_vline(aes(xintercept=5500), colour="blue", linetype="dashed") + 
    ylab("Cell Number")
# nCount_RNA
h2_df.M2 <- data.frame(seurat.M2$sample, seurat.M2$nCount_RNA)
colnames(h2_df.M2) <- c("Sample", "nCount_RNA")
VlnPlot(seurat.M2, features = "nCount_RNA", pt.size = 0.1) + 
    geom_hline(aes(yintercept=2500), colour="blue", linetype="dashed") +
    geom_hline(aes(yintercept=50000), colour="blue", linetype="dashed")
ggplot(h2_df.M2, aes(x=nCount_RNA, color=Sample)) + 
    geom_histogram(binwidth = 100, fill="white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    geom_vline(aes(xintercept=2500), colour="black", linetype="dashed") + 
    geom_vline(aes(xintercept=50000), colour="black", linetype="dashed") + 
    ylab("Cell Number")
# percent.mt
h3_df.M2 <- data.frame(seurat.M2$sample, seurat.M2$percent.mt)
colnames(h3_df.M2) <- c("Sample", "percent.mt")
VlnPlot(seurat.M2, features = "percent.mt", pt.size = 0.1) + 
    geom_hline(aes(yintercept=12), colour="blue", linetype="dashed")
ggplot(h3_df.M2, aes(x=percent.mt, color=Sample)) + 
    geom_histogram(binwidth = 1, fill="white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    ylab("Cell Number") + geom_vline(aes(xintercept=12), colour="blue", linetype="dashed")
# feature scatter
FeatureScatter(seurat.M2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
    geom_vline(aes(xintercept=2500), colour="blue", linetype="dashed") + 
    geom_vline(aes(xintercept=50000), colour="blue", linetype="dashed") + 
    geom_hline(aes(yintercept=800), colour="blue", linetype="dashed") + 
    geom_hline(aes(yintercept=5500), colour="blue", linetype="dashed")
# filter  
seurat.M2_filter <- subset(seurat.M2, subset = nFeature_RNA > 800 & nFeature_RNA < 5500 & nCount_RNA > 2500 & nCount_RNA < 50000 & percent.mt < 12)


############ M3
M3 <- Read10X(paste0(path.data, 'M3'))
dimnames(M3)[[2]] <- paste('M3', dimnames(M3)[[2]], sep = '_')
seurat.M3 <- CreateSeuratObject(counts = M3)
seurat.M3@meta.data$sample <- rep('M3', dim(M3)[2])
seurat.M3[["percent.mt"]] <- PercentageFeatureSet(seurat.M3, pattern = "^MT-")
# nFeature_RNA
h_df.M3 <- data.frame(seurat.M3$sample, seurat.M3$nFeature_RNA)
colnames(h_df.M3) <- c("Sample", "nFeature_RNA")
VlnPlot(seurat.M3, features = "nFeature_RNA", pt.size = 0.1) + 
    geom_hline(aes(yintercept=1000), colour="blue", linetype="dashed") +
    geom_hline(aes(yintercept=4900), colour="blue", linetype="dashed")
ggplot(h_df.M3, aes(x=nFeature_RNA, color=Sample)) + 
    geom_histogram(binwidth = 20, fill = "white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    geom_vline(aes(xintercept=1000), colour="blue", linetype="dashed") + 
    geom_vline(aes(xintercept=4900), colour="blue", linetype="dashed") + 
    ylab("Cell Number")
# nCount_RNA
h2_df.M3 <- data.frame(seurat.M3$sample, seurat.M3$nCount_RNA)
colnames(h2_df.M3) <- c("Sample", "nCount_RNA")
VlnPlot(seurat.M3, features = "nCount_RNA", pt.size = 0.1) + 
    geom_hline(aes(yintercept=2500), colour="blue", linetype="dashed") +
    geom_hline(aes(yintercept=30000), colour="blue", linetype="dashed")
ggplot(h2_df.M3, aes(x=nCount_RNA, color=Sample)) + 
    geom_histogram(binwidth = 100, fill="white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    geom_vline(aes(xintercept=2500), colour="black", linetype="dashed") + 
    geom_vline(aes(xintercept=30000), colour="black", linetype="dashed") + 
    ylab("Cell Number")
# percent.mt
h3_df.M3 <- data.frame(seurat.M3$sample, seurat.M3$percent.mt)
colnames(h3_df.M3) <- c("Sample", "percent.mt")
VlnPlot(seurat.M3, features = "percent.mt", pt.size = 0.1) + 
    geom_hline(aes(yintercept=9), colour="blue", linetype="dashed")
ggplot(h3_df.M3, aes(x=percent.mt, color=Sample)) + 
    geom_histogram(binwidth = 1, fill="white") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    ylab("Cell Number") + geom_vline(aes(xintercept=9), colour="blue", linetype="dashed")
# feature scatter
FeatureScatter(seurat.M3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
    geom_vline(aes(xintercept=2500), colour="blue", linetype="dashed") + 
    geom_vline(aes(xintercept=30000), colour="blue", linetype="dashed") + 
    geom_hline(aes(yintercept=1000), colour="blue", linetype="dashed") + 
    geom_hline(aes(yintercept=4900), colour="blue", linetype="dashed")
# filter  
seurat.M3_filter <- subset(seurat.M3, subset = nFeature_RNA > 1000 & nFeature_RNA < 4900 & nCount_RNA > 2500 & nCount_RNA < 30000 & percent.mt < 9)


seurat.merge <- merge(seurat.C1_filter, 
                      c(seurat.C2_filter, seurat.C3_filter, seurat.C4_filter,
                        seurat.M1_filter, seurat.M2_filter, seurat.M3_filter))

# status
samples <- seurat.merge$sample
vec.st <- rep('0', length(samples))
vec.st[samples %in% c('C1', 'C2', 'C3', 'C4')] <- 'Normal'
vec.st[samples %in% c('M1', 'M2', 'M3')] <- 'Microtia'
seurat.merge$status <- vec.st

file.ear <- '/home/disk/drizzle/wgk/data/ear.filtered.Rdata'
saveRDS(seurat.merge, file = file.ear)

###############
# combined_filter for downstream
# dim(combined_filter)
# [1] 31733 37703
###############

###############