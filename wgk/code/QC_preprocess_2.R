# parameter
# path of input data
path.data <- '/local/hanyang/data/test11@qq.com/single-cell/'
# output path
path.QC.PP <- '/mdshare/node9/zy/wgk/out_QC_PP/'
# QC
feature.min <- 800
feature.max <- 6000
UMI.min <- 2000
UMI.max <- 50000
percent.mt.max <- 15
# seed
num.seed <- 123
### dimension reduction
num.features <- 4000
num_PCs <- 50
n.neighbors <- 30
min.dist <- 0.3
# cluster
use.reduction <- "pca"
use.dims <- 50
cluster.k <- 20
cluster.resolution <- 2


QC_preprocess <- function(path.data, path.QC.PP, 
                          feature.min, feature.max, UMI.min, UMI.max, percent.mt.max, 
                          num.seed, num.features, num_PCs, n.neighbors, min.dist,
                          use.reduction, use.dims, cluster.k, cluster.resolution) {
    
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    
    # parameter
    # QC
    feature.min <- as.numeric(feature.min)
    feature.max <- as.numeric(feature.max)
    UMI.min <- as.numeric(UMI.min)
    UMI.max <- as.numeric(UMI.max)
    percent.mt.max <- as.numeric(percent.mt.max)
    # seed
    num.seed <- as.numeric(num.seed)
    ### dimension reduction
    num.features <- as.numeric(num.features)
    num_PCs <- as.numeric(num_PCs)
    n.neighbors <- as.numeric(n.neighbors)
    min.dist <- as.numeric(min.dist)
    # cluster
    use.dims <- as.numeric(use.dims)
    cluster.k <- as.numeric(cluster.k)
    cluster.resolution <- as.numeric(cluster.resolution)
    
    ## read
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
    
    seurat.all_filter <- 
        subset(seurat.merge, 
               subset = nFeature_RNA > feature.min & nFeature_RNA < feature.max & nCount_RNA > UMI.min & nCount_RNA < UMI.max & percent.mt < percent.mt.max)

    # show results of QC
    # nFeature_RNA
    h_df <- data.frame(seurat.merge$sample, seurat.merge$nFeature_RNA)
    colnames(h_df) <- c("Sample", "nFeature_RNA")
    # VlnPlot(seurat.merge, features = "nFeature_RNA", pt.size = 0.1, group.by = 'sample') + 
    #     geom_hline(aes(yintercept=feature.min), colour="blue", linetype="dashed") +
    #     geom_hline(aes(yintercept=feature.max), colour="blue", linetype="dashed")
    plot_feature <- 
        ggplot(h_df, aes(x=nFeature_RNA, color=Sample)) + 
        geom_histogram(binwidth = 20, fill = "white") + 
        theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
        geom_vline(aes(xintercept=feature.min), colour="blue", linetype="dashed") + 
        geom_vline(aes(xintercept=feature.max), colour="blue", linetype="dashed") + 
        ylab("Cell Number")
    ggsave(plot = plot_feature, filename = 'Feature.png', path = path.QC.PP)
    # nCount_RNA
    h2_df <- data.frame(seurat.merge$sample, seurat.merge$nCount_RNA)
    colnames(h2_df) <- c("Sample", "nCount_RNA")
    # VlnPlot(seurat.merge, features = "nCount_RNA", pt.size = 0.001) + 
    #     geom_hline(aes(yintercept=UMI.min), colour="blue", linetype="dashed") +
    #     geom_hline(aes(yintercept=UMI.max), colour="blue", linetype="dashed")
    plot_UMI <- 
        ggplot(h2_df, aes(x=nCount_RNA, color=Sample)) + 
        geom_histogram(binwidth = 100, fill="white") + 
        theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
        geom_vline(aes(xintercept=UMI.min), colour="black", linetype="dashed") + 
        geom_vline(aes(xintercept=UMI.max), colour="black", linetype="dashed") + 
        ylab("Cell Number")
    ggsave(plot = plot_UMI, filename = 'UMI.png', path = path.QC.PP)
    # percent.mt
    h3_df <- data.frame(seurat.merge$sample, seurat.merge$percent.mt)
    colnames(h3_df) <- c("Sample", "percent.mt")
    # VlnPlot(seurat.merge, features = "percent.mt", pt.size = 0.01) + 
    #     geom_hline(aes(yintercept=percent.mt.max), colour="blue", linetype="dashed")
    plot_mt <- 
        ggplot(h3_df, aes(x=percent.mt, color=Sample)) + 
        geom_histogram(binwidth = 1, fill="white") + 
        theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
        ylab("Cell Number") + 
        geom_vline(aes(xintercept=percent.mt.max), colour="blue", linetype="dashed")
    ggsave(plot = plot_mt, filename = 'MT.png', path = path.QC.PP)
    # feature scatter
    # FeatureScatter(seurat.merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
    #     geom_vline(aes(xintercept=4000), colour="blue", linetype="dashed") + 
    #     geom_vline(aes(xintercept=35000), colour="blue", linetype="dashed") + 
    #     geom_hline(aes(yintercept=1000), colour="blue", linetype="dashed") + 
    #     geom_hline(aes(yintercept=5000), colour="blue", linetype="dashed")

    ### normalization and dimension reduction
    set.seed(num.seed)
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
    
    # show results of TNSE and UMAP
    plot.TSNE <- DimPlot(Obj.seurat.input, group.by = "sample", reduction = 'tsne')
    ggsave(plot = plot.TSNE, filename = 'TSNE.png', path = path.QC.PP)
    plot.UMAP <- DimPlot(Obj.seurat.input, group.by = "sample", reduction = 'umap')
    ggsave(plot = plot.UMAP, filename = 'UMAP.png', path = path.QC.PP)
    
    # cluster
    Obj.seurat.input <- FindNeighbors(Obj.seurat.input, reduction = use.reduction, 
                                      dims = 1:use.dims, k.param = cluster.k)
    Obj.seurat.input <- FindClusters(Obj.seurat.input, resolution = cluster.resolution)
    
    # show cluster results
    plot_cluster <- DimPlot(Obj.seurat.input, group.by = "seurat_clusters", reduction = 'umap')
    ggsave(plot = plot_cluster, filename = 'Cluster.png', path = path.QC.PP)
    
}

args <- commandArgs(T)
QC_preprocess(args[1], args[2], args[3], args[4], args[5], args[6], 
                  args[7], args[8], args[9], args[10], args[11], args[12], 
                  args[13], args[14], args[15], args[16])

# Rscript QC_preprocess_2.R '/local/hanyang/data/test11@qq.com/single-cell/' '/mdshare/node9/zy/wgk/out_QC_PP/' 800 6000 2000 50000 15 123 4000 50 30 0.3 "pca" 50 20 2