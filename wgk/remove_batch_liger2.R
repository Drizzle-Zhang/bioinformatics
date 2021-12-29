# .libPaths("/home/zy/tools/R-4.0.0/library")
library(Seurat)
library(liger)
library(SeuratWrappers)
library(cowplot)

file.ear <- '/home/disk/drizzle/wgk/data/ear.filtered.Rdata'
seurat.ear.liger <- readRDS(file.ear)
cells <- colnames(seurat.ear.liger@assays$RNA@counts)

# first batch
vec.first.N <- c('C1', 'C2', 'C3')
first.cell.N <- cells[seurat.ear.liger$sample %in% vec.first.N]
seurat.first <- subset(seurat.ear.liger, cells = first.cell.N)
seurat.first <- NormalizeData(seurat.first)
seurat.first <- FindVariableFeatures(seurat.first, nfeatures = 5000)
# do.center = FALSE, NMF needs non-negative input
seurat.first <- ScaleData(seurat.first, split.by = "sample", do.center = FALSE)
# k = 20 for default, large k leads to more aggressive correction
k <- 20
seurat.first <- RunOptimizeALS(seurat.first, k = k, lambda = 5, split.by = "sample")
seurat.first <- RunQuantileNorm(seurat.first, split.by = "sample")
# reduction using iNMF
seurat.first <- FindNeighbors(seurat.first, reduction = "iNMF", dims = 1:k)
seurat.first <- FindClusters(seurat.first, resolution = 0.5)
seurat.first <- RunUMAP(seurat.first, dims = 1:k, 
                         reduction = "iNMF", n.neighbors = 30)
file.liger <- '/home/disk/drizzle/wgk/data/ear.liger_first_normal_20_5.Rdata'
saveRDS(seurat.first, file = file.liger)
seurat.first <- readRDS(file.liger)

DimPlot(seurat.first, group.by = "sample")
DimPlot(seurat.first, group.by = "seurat_clusters", label = T)
DimPlot(seurat.first, group.by = "status", label = T)
scMAGIC.label.1 <- readRDS('/home/zy/scRef/ear/scMAGIC.2.Rdata')
seurat.first@meta.data$scMAGIC.label.1 <- scMAGIC.label.1[first.cell,'scRef.tag']
DimPlot(seurat.first, group.by = "scMAGIC.label.1", label = T)

# normal
cells <- row.names(seurat.first@reductions$umap@cell.embeddings)
normal.cell <- cells[seurat.first$status == 'Normal']
DimPlot(SubsetData(seurat.first, cells = normal.cell), group.by = "sample", pt.size = 0.3)
# microtia
microtia.cell <- cells[seurat.first$status == 'Microtia']
DimPlot(SubsetData(seurat.first, cells = microtia.cell), group.by = "sample")

library(ggplot2)
P1 <- DimPlot(subset(seurat.first, subset = sample == 'C1'), group.by = "sample") + ylim(-10, 15)
P2 <- DimPlot(subset(seurat.first, subset = sample == 'C2'), group.by = "sample") + ylim(-10, 15)
P3 <- DimPlot(subset(seurat.first, subset = sample == 'C3'), group.by = "sample") + ylim(-10, 15)
P1 + P2 + P3

# find marker
all.marker <- FindAllMarkers(seurat.first)
all.marker$diff.pct <- abs(all.marker$pct.2 - all.marker$pct.1)
marker.5 <- all.marker[all.marker$cluster == 5,]
marker.7 <- all.marker[all.marker$cluster == 7,]
marker.14 <- all.marker[all.marker$cluster == 14,]

file.liger <- '/home/disk/drizzle/wgk/data/ear.liger_first_50_1.Rdata'
saveRDS(seurat.first, file = file.liger)
file.liger <- '/home/disk/drizzle/wgk/data/ear.liger_first_50_3.Rdata'
saveRDS(seurat.first, file = file.liger)
file.liger <- '/home/disk/drizzle/wgk/data/ear.liger_first_50_5.Rdata'
saveRDS(seurat.first, file = file.liger)


