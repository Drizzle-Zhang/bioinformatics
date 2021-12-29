library(Seurat)
library(liger)
library(SeuratWrappers)
library(cowplot)

file.ear <- '/home/disk/drizzle/wgk/data/ear.filtered.Rdata'
seurat.ear.liger <- readRDS(file.ear)
cells <- colnames(seurat.ear.liger@assays$RNA@counts)
# normal
normal.cell <- cells[seurat.ear.liger$status == 'Normal']
seurat.normal <- subset(seurat.ear.liger, cells = normal.cell)
seurat.normal <- NormalizeData(seurat.normal)
seurat.normal <- FindVariableFeatures(seurat.normal, nfeatures = 4000)
# do.center = FALSE, NMF needs non-negative input
seurat.normal <- ScaleData(seurat.normal, split.by = "sample", do.center = FALSE)
# k = 20 for default, large k leads to more aggressive correction
seurat.normal <- RunOptimizeALS(seurat.normal, k = 50, lambda = 5, split.by = "sample")
seurat.normal <- RunQuantileNorm(seurat.normal, split.by = "sample")
# reduction using iNMF
seurat.normal <- FindNeighbors(seurat.normal, reduction = "iNMF_raw", dims = 1:50)
seurat.normal <- FindClusters(seurat.normal, resolution = 1)
seurat.normal <- RunUMAP(seurat.normal, dims = 1:ncol(seurat.normal[["iNMF_raw"]]), 
                            reduction = "iNMF_raw", n.neighbors = 30)

# after batch correction
file.ear <- '/home/disk/drizzle/wgk/data/ear.filtered.Rdata'
seurat.ear.liger <- readRDS(file.ear)
seurat.ear.liger <- NormalizeData(seurat.ear.liger)
seurat.ear.liger <- FindVariableFeatures(seurat.ear.liger, nfeatures = 4000)
# do.center = FALSE, NMF needs non-negative input
seurat.ear.liger <- ScaleData(seurat.ear.liger, split.by = "sample", do.center = FALSE)
# k = 20 for default, large k leads to more aggressive correction
seurat.ear.liger <- RunOptimizeALS(seurat.ear.liger, k = 20, lambda = 5, split.by = "sample")
seurat.ear.liger <- RunQuantileNorm(seurat.ear.liger, split.by = "sample")
# reduction using iNMF
seurat.ear.liger <- FindNeighbors(seurat.ear.liger, reduction = "iNMF_raw", dims = 1:20)
seurat.ear.liger <- FindClusters(seurat.ear.liger, resolution = 1)
seurat.ear.liger <- RunUMAP(seurat.ear.liger, dims = 1:ncol(seurat.ear.liger[["iNMF_raw"]]), 
                            reduction = "iNMF_raw", n.neighbors = 30)
# plot
DimPlot(seurat.ear.liger, group.by = "sample")
DimPlot(seurat.ear.liger, group.by = "seurat_clusters", label = T)
DimPlot(seurat.ear.liger, group.by = "status", label = T)
# normal
cells <- row.names(seurat.ear.liger@reductions$umap@cell.embeddings)
normal.cell <- cells[seurat.ear.liger$status == 'Normal']
DimPlot(SubsetData(seurat.ear.liger, cells = normal.cell), group.by = "sample")
normal.cell <- cells[seurat.ear.liger$sample %in% c('C1', 'C2', 'C3')]
DimPlot(SubsetData(seurat.ear.liger, cells = normal.cell), group.by = "sample")

# microtia
microtia.cell <- cells[seurat.ear.liger$status == 'Microtia']
DimPlot(SubsetData(seurat.ear.liger, cells = microtia.cell), group.by = "sample")

file.liger <- '/home/disk/drizzle/wgk/data/ear.liger_20_5.Rdata'
saveRDS(seurat.ear.liger, file = file.liger)
