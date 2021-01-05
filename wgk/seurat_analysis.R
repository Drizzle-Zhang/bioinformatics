# Liger for batch correction
# With SeuratWrappers, directly run liger on Seurat object
# http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/liger.html

library(Seurat)
library(liger)
library(SeuratWrappers)
library(cowplot)

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


file.seurat <- '/home/disk/drizzle/wgk/data/ear.seurat.Rdata'
saveRDS(seurat.ear, file = file.seurat)
seurat.ear <- readRDS(file.seurat)

FeaturePlot(seurat.ear, features = c('COL1A1', 'COL1A2', 'PDGFRB', 'COL'))
FeaturePlot(seurat.ear, features = c('COL2A1', 'COL9A2', 'ACAN', 'SOX9'))

cells <- row.names(seurat.ear@reductions$umap@cell.embeddings)
C4 <- cells[seurat.ear$sample %in% c('C4')]
DimPlot(subset(seurat.ear, cells = C4), group.by = "sample")

# first batch
file.ear <- '/home/disk/drizzle/wgk/data/ear.filtered.Rdata'
seurat.ear <- readRDS(file.ear)
seurat.first <- subset(seurat.ear, subset = sample %in% c('C1', 'C2', 'C3', 'M1', 'M2'))
seurat.first <- NormalizeData(seurat.first)
seurat.first <- FindVariableFeatures(seurat.first, nfeatures = 5000)
seurat.first <- ScaleData(seurat.first, split.by = "sample")
seurat.first <- RunPCA(seurat.first, verbose = F, npcs = 100)
seurat.first <- RunUMAP(seurat.first, dims = 1:100, n.neighbors = 50)
seurat.first <- FindNeighbors(seurat.first, reduction = "pca", dims = 1:100)
seurat.first <- FindClusters(seurat.first, resolution = 0.15)

file.seurat <- '/home/disk/drizzle/wgk/data/ear.seurat.first.Rdata'
saveRDS(seurat.first, file = file.seurat)
seurat.first <- readRDS(file.seurat)

# cluster
table(seurat.first$sample, seurat.first$RNA_snn_res.0.8)/as.vector(table(seurat.first$sample))

DimPlot(seurat.first, group.by = "sample")
DimPlot(seurat.first, group.by = "seurat_clusters", label = T)
DimPlot(seurat.first, group.by = "status", label = T)
# normal
cells <- row.names(seurat.first@reductions$umap@cell.embeddings)
normal.cell <- cells[seurat.first$status == 'Normal']
DimPlot(SubsetData(seurat.first, cells = normal.cell), group.by = "sample")
C3.cell <- cells[seurat.first$sample == 'C3']
DimPlot(SubsetData(seurat.first, cells = C3.cell), group.by = "sample")


FeaturePlot(seurat.first, features = c('COL1A1', 'COL1A2', 'PDGFRB', 'COL5A2'))
FeaturePlot(seurat.first, features = c('COL2A1', 'COL9A2', 'ACAN', 'SOX9'))


# SCTranform
file.ear <- '/home/disk/drizzle/wgk/data/ear.filtered.Rdata'
seurat.ear.all <- readRDS(file.ear)
pbmc <- SCTransform(pbmc)