library(Seurat)

# read data
setwd('/home/zy/scRef/cross_validation/hPBMC')
hpbmc <- readRDS('/home/zy/scRef/cross_validation/hPBMC/pbmc68k_data.rds')
data.filter <- (t(hpbmc$all_data$`17820`$hg19$mat))
dimnames(data.filter)[[1]] <- hpbmc$all_data$`17820`$hg19$gene_symbols
dimnames(data.filter)[[2]] <- hpbmc$all_data$`17820`$hg19$barcodes

label.filter <- read.delim('./68k_pbmc_barcodes_annotation.tsv', row.names = 1)
label.filter <- as.matrix(label.filter)[dimnames(data.filter)[[2]],]

# sample
sub.barcode <- sample(dimnames(data.filter)[[2]], 10000)
data.filter <- data.filter[, sub.barcode]
label.filter <- label.filter[sub.barcode]

# data preparing
seurat.unlabeled <- CreateSeuratObject(counts = data.filter)
seurat.unlabeled <- NormalizeData(seurat.unlabeled, normalization.method = "LogNormalize", 
                                  scale.factor = 10000)
seurat.unlabeled <- FindVariableFeatures(seurat.unlabeled, selection.method = "vst", nfeatures = 2000)
# VariableFeatures(seurat.unlabeled)
# VariableFeaturePlot(seurat.unlabeled)
# seurat.unlabeled[["percent.mt"]] <- PercentageFeatureSet(seurat.unlabeled, pattern = "^mt-")
# seurat.unlabeled <- ScaleData(seurat.unlabeled, vars.to.regress = "percent.mt")
# seurat.unlabeled <- ScaleData(seurat.unlabeled, features = all.genes)
seurat.unlabeled <- ScaleData(seurat.unlabeled)

# add label
seurat.unlabeled@meta.data$original.label <- label.filter
# mtx.tag <- as.matrix(meta.tag)
# seurat.unlabeled@meta.data$scRef.tag <- mtx.tag[use.cells, 'scRef.tag'] 
# seurat.unlabeled@meta.data$new.tag <- new.tag

# PCA
seurat.unlabeled <- RunPCA(seurat.unlabeled, npcs = 75)

# cluster
# seurat.unlabeled <- FindNeighbors(seurat.unlabeled, reduction = "pca", dims = 1:75, nn.eps = 0.5)
# seurat.unlabeled <- FindClusters(seurat.unlabeled, resolution = 3, n.start = 20)

# UMAP
seurat.unlabeled <- RunUMAP(seurat.unlabeled, dims = 1:20, n.neighbors = 30)
# figure1: ture label
DimPlot(seurat.unlabeled, reduction = "umap", label = T, group.by = 'original.label')
