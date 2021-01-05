library(Seurat)
library(ggplot2)
library(patchwork)
options(future.globals.maxSize = 7000 * 1024^2)

file.sct <- '/home/disk/drizzle/wgk/data/ear.sct.Rdata'
list.sct <- readRDS(file.sct)
list.sample <- list.sct$list.sample
sel.features <- list.sct$sel.features

# normal
normal.list <- list.sample[c('C1', 'C2', 'C3', 'C4')]
normal.list <- PrepSCTIntegration(object.list = normal.list, anchor.features = sel.features, 
                                  verbose = FALSE)
normal.anchors <- FindIntegrationAnchors(object.list = normal.list, normalization.method = "SCT", 
                                         anchor.features = sel.features, verbose = FALSE)
normal.integrated <- IntegrateData(anchorset = normal.anchors, normalization.method = "SCT", 
                                   verbose = FALSE)
# plot
normal.integrated <- RunPCA(normal.integrated, npcs = 100, verbose = FALSE)
normal.integrated <- RunUMAP(normal.integrated, dims = 1:50, n.neighbors = 30, min.dist = 0.05)
DimPlot(normal.integrated, group.by = c("sample"))
DimPlot(normal.integrated, split.by = 'sample', ncol = 2)
# save
file.normal <- '/home/disk/drizzle/wgk/data/seurat.normal.Rdata'
saveRDS(normal.integrated, file = file.normal)
normal.integrated <- readRDS(file.normal)


normal.integrated@assays$SCT <- NULL
normal.integrated@assays$RNA <- NULL
list.N.M <- list()
list.N.M$N <- normal.integrated
list.N.M$M <- list.sample$M1
list.N.M <- PrepSCTIntegration(object.list = list.N.M, anchor.features = sel.features, 
                               verbose = FALSE)
merge.anchors <- FindIntegrationAnchors(object.list = list.N.M, normalization.method = "SCT", 
                                        anchor.features = sel.features)
merge.integrated <- IntegrateData(anchorset = merge.anchors, normalization.method = "SCT")
# plot
merge.integrated <- RunPCA(merge.integrated, npcs = 100, verbose = FALSE, features = sel.features)
merge.integrated <- RunUMAP(merge.integrated, dims = 1:100, n.neighbors = 20)
DimPlot(merge.integrated, group.by = c("sample"))
DimPlot(merge.integrated, split.by = 'sample', ncol = 3)
FeaturePlot(merge.integrated, features = c('COL1A1', 'COL1A2', 'PDGFRB', 'COL5A2'))
