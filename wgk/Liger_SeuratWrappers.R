# Liger for batch correction
# With SeuratWrappers, directly run liger on Seurat object
# http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/liger.html

library(Seurat)
library(liger)
library(SeuratWrappers)
library(cowplot)

########################
# before batch correction
file.ear <- '/home/disk/drizzle/wgk/data/ear.filtered.Rdata'
seurat.ear <- readRDS(file.ear)
seurat.ear <- NormalizeData(seurat.ear)
seurat.ear <- FindVariableFeatures(seurat.ear, nfeatures = 5000)
seurat.ear <- ScaleData(seurat.ear, split.by = "sample")
seurat.ear <- RunPCA(seurat.ear, verbose = F, npcs = 150)
seurat.ear <- RunUMAP(seurat.ear, dims = 1:100, n.neighbors = 30)
# plot
P1 <- DimPlot(seurat.ear, group.by = "sample")
DimPlot(seurat.ear, group.by = "sample")
DimPlot(seurat.ear.liger, group.by = "seurat_clusters", label = T)
DimPlot(seurat.ear, group.by = "status", label = T)
# normal
cells <- row.names(seurat.ear@reductions$umap@cell.embeddings)
normal.cell <- cells[seurat.ear$status == 'Normal']
DimPlot(SubsetData(seurat.ear, cells = normal.cell), group.by = "sample")
normal.cell <- cells[seurat.ear$sample %in% c('C1', 'C2', 'C3')]
DimPlot(SubsetData(seurat.ear, cells = normal.cell), group.by = "sample")
# microtia
microtia.cell <- cells[seurat.ear$status == 'Microtia']
DimPlot(SubsetData(seurat.ear, cells = microtia.cell), group.by = "sample")
# first
first.cell <- cells[seurat.ear$sample %in% c('C1', 'C2', 'C3', 'M1', 'M2')]
DimPlot(SubsetData(seurat.ear, cells = first.cell), group.by = "sample")
######################################

file.seurat <- '/home/disk/drizzle/wgk/data/ear.seurat.Rdata'
saveRDS(seurat.ear, file = file.seurat)
seurat.ear <- readRDS(file.seurat)

FeaturePlot(seurat.ear, features = c('COL1A1', 'COL1A2', 'PDGFRB', 'COL'))
FeaturePlot(seurat.ear, features = c('COL2A1', 'COL9A2', 'ACAN', 'SOX9'))

cells <- row.names(seurat.ear@reductions$umap@cell.embeddings)
C4 <- cells[seurat.ear$sample %in% c('C4')]
DimPlot(subset(seurat.ear, cells = C4), group.by = "sample")

# after batch correction
file.ear <- '/home/disk/drizzle/wgk/data/ear.filtered.Rdata'
seurat.ear.liger <- readRDS(file.ear)
seurat.ear.liger <- NormalizeData(seurat.ear.liger)
seurat.ear.liger <- FindVariableFeatures(seurat.ear.liger, nfeatures = 3000)
# do.center = FALSE, NMF needs non-negative input
seurat.ear.liger <- ScaleData(seurat.ear.liger, split.by = "sample", do.center = FALSE)
# k = 20 for default, large k leads to more aggressive correction
k = 20
seurat.ear.liger <- RunOptimizeALS(seurat.ear.liger, k = k, lambda = 5, split.by = "sample")
seurat.ear.liger <- RunQuantileNorm(seurat.ear.liger, split.by = "sample")
# reduction using iNMF
seurat.ear.liger <- FindNeighbors(seurat.ear.liger, reduction = "iNMF", dims = 1:k)
seurat.ear.liger <- FindClusters(seurat.ear.liger, resolution = 0.5)
seurat.ear.liger <- RunUMAP(seurat.ear.liger, dims = 1:ncol(seurat.ear.liger[["iNMF"]]), 
                            reduction = "iNMF", n.neighbors = 30)
file.liger <- '/home/disk/drizzle/wgk/data/ear.liger_3000_20_5.Rdata'
saveRDS(seurat.ear.liger, file = file.liger)

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

# feature plot
FeaturePlot(seurat.ear.liger, features = c('COL1A1', 'COL1A2', 'PDGFRB', 'COL'))
FeaturePlot(seurat.ear.liger, features = c('COL2A1', 'COL9A2', 'ACAN', 'SOX9'))


################################################
# after batch correction
file.ear <- '/home/disk/drizzle/wgk/data/ear.filtered.Rdata'
seurat.ear.liger <- readRDS(file.ear)
seurat.ear.liger <- NormalizeData(seurat.ear.liger)
seurat.ear.liger <- FindVariableFeatures(seurat.ear.liger, nfeatures = 4000)
# do.center = FALSE, NMF needs non-negative input
seurat.ear.liger <- ScaleData(seurat.ear.liger, split.by = "sample", do.center = FALSE)
# k = 20 for default, large k leads to more aggressive correction
seurat.ear.liger <- RunOptimizeALS(seurat.ear.liger, k = 40, lambda = 5, split.by = "sample")
seurat.ear.liger <- RunQuantileNorm(seurat.ear.liger, split.by = "sample")
# reduction using iNMF
seurat.ear.liger <- FindNeighbors(seurat.ear.liger, reduction = "iNMF_raw", dims = 1:40)
seurat.ear.liger <- FindClusters(seurat.ear.liger, resolution = 1)
seurat.ear.liger <- RunUMAP(seurat.ear.liger, dims = 1:ncol(seurat.ear.liger[["iNMF_raw"]]), 
                            reduction = "iNMF_raw", n.neighbors = 30)

file.liger <- '/home/disk/drizzle/wgk/data/ear.liger_40_5.Rdata'
saveRDS(seurat.ear.liger, file = file.liger)

seurat.ear.liger <- readRDS(file.liger)

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

# feature plot
FeaturePlot(seurat.ear.liger, features = c('COL1A1', 'COL1A2', 'PDGFRB', 'COL'))
FeaturePlot(seurat.ear.liger, features = c('COL2A1', 'COL9A2', 'ACAN', 'SOX9'))
