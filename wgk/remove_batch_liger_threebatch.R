# .libPaths("/home/zy/tools/R-4.0.0/library")

library(Seurat)
library(liger)
library(SeuratWrappers)
library(cowplot)

file.ear <- '/home/disk/drizzle/wgk/data/ear.filtered.Rdata'
seurat.ear.liger <- readRDS(file.ear)
vec.batch <- rep('0', length(seurat.ear.liger$sample))
vec.batch[seurat.ear.liger$sample %in% c('C1', 'C2', 'C3', 'M1', 'M2')] <- 'batch1'
vec.batch[seurat.ear.liger$sample == 'M3'] <- 'batch2'
vec.batch[seurat.ear.liger$sample == 'C4'] <- 'batch3'
seurat.ear.liger@meta.data$batch <- vec.batch

seurat.ear.liger <- NormalizeData(seurat.ear.liger)
seurat.ear.liger <- FindVariableFeatures(seurat.ear.liger, nfeatures = 4000)
# do.center = FALSE, NMF needs non-negative input
seurat.ear.liger <- ScaleData(seurat.ear.liger, split.by = "sample", do.center = FALSE)
# k = 20 for default, large k leads to more aggressive correction
k = 20
seurat.ear.liger <- RunOptimizeALS(seurat.ear.liger, k = k, lambda = 3, split.by = "batch")
seurat.ear.liger <- RunQuantileNorm(seurat.ear.liger, split.by = "batch")
# reduction using iNMF
seurat.ear.liger <- FindNeighbors(seurat.ear.liger, reduction = "iNMF", dims = 1:k)
seurat.ear.liger <- FindClusters(seurat.ear.liger, resolution = 0.15)
seurat.ear.liger <- RunUMAP(seurat.ear.liger, dims = 1:ncol(seurat.ear.liger[["iNMF"]]), 
                            reduction = "iNMF", n.neighbors = 30)
# file.liger <- '/home/disk/drizzle/wgk/data/ear.liger_3000_20_5.Rdata'
file.liger <- '/home/disk/drizzle/wgk/data/ear.liger_three_20_3.Rdata'
saveRDS(seurat.ear.liger, file = file.liger)
seurat.ear.liger <- readRDS(file.liger)

# cluster
table(seurat.ear.liger$sample, seurat.ear.liger$RNA_snn_res.0.15)/as.vector(table(seurat.ear.liger$sample))
mytable <- table(seurat.ear.liger$sample, seurat.ear.liger$RNA_snn_res.0.15)[,1:4]
mytable/rowSums(mytable)

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

# first
DimPlot(subset(seurat.ear.liger, subset = sample %in% c('C1', 'C2', 'C3', 'M1', 'M2')), 
        group.by = "sample")

library(ggplot2)
P1 <- DimPlot(subset(seurat.ear.liger, subset = sample == 'C1'), group.by = "sample") + ylim(-10, 15)
P2 <- DimPlot(subset(seurat.ear.liger, subset = sample == 'C2'), group.by = "sample") + ylim(-10, 15)
P3 <- DimPlot(subset(seurat.ear.liger, subset = sample == 'C3'), group.by = "sample") + ylim(-10, 15)
P4 <- DimPlot(subset(seurat.ear.liger, subset = sample == 'C4'), group.by = "sample") + ylim(-10, 15)
P1 + P2 + P3 + P4

P1 <- DimPlot(subset(seurat.ear.liger, subset = sample == 'M1'), group.by = "sample") + ylim(-10, 15)
P2 <- DimPlot(subset(seurat.ear.liger, subset = sample == 'M2'), group.by = "sample") + ylim(-10, 15)
P3 <- DimPlot(subset(seurat.ear.liger, subset = sample == 'M3'), group.by = "sample") + ylim(-10, 15)
P1 + P2 + P3

# feature plot
FeaturePlot(seurat.ear.liger, features = c('COL1A1', 'COL1A2', 'PDGFRB', 'COL'))
FeaturePlot(seurat.ear.liger, features = c('COL2A1', 'COL9A2', 'ACAN', 'SOX9'))
