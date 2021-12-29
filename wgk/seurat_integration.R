library(Seurat)
library(ggplot2)
library(patchwork)
options(future.globals.maxSize = 7000 * 1024^2)

file.ear <- '/home/disk/drizzle/wgk/data/ear.filtered.Rdata'
seurat.ear <- readRDS(file.ear)

list.sample <- SplitObject(seurat.ear, split.by = 'sample')
for (i in 1:length(list.sample)) {
    list.sample[[i]] <- SCTransform(list.sample[[i]], verbose = FALSE)
}
sel.features <- SelectIntegrationFeatures(object.list = list.sample, nfeatures = 4000)

list.sct <- list()
list.sct$list.sample <- list.sample
list.sct$sel.features <- sel.features
file.sct <- '/home/disk/drizzle/wgk/data/ear.sct.Rdata'
saveRDS(list.sct, file = file.sct)
list.sct <- readRDS(file.sct)
list.sample <- list.sct$list.sample
sel.features <- list.sct$sel.features

# C1, C2, C3, C4
normal.list <- list.sample[c('C1', 'C2', 'C3', 'C4')]
normal.anchors <- FindIntegrationAnchors(object.list = normal.list, normalization.method = "SCT", 
                                         anchor.features = sel.features, verbose = FALSE)
normal.integrated <- IntegrateData(anchorset = normal.anchors, normalization.method = "SCT", 
                                   verbose = FALSE)


# C1, C2, C3 + C4
# C1, C2, C3
normal.list.1 <- list.sample[c('C1', 'C2', 'C3')]
normal.list.1 <- PrepSCTIntegration(object.list = normal.list.1, anchor.features = sel.features, 
                                    verbose = FALSE)
normal.anchors.1 <- FindIntegrationAnchors(object.list = normal.list.1, normalization.method = "SCT", 
                                           anchor.features = sel.features, verbose = FALSE)
normal.integrated.1 <- IntegrateData(anchorset = normal.anchors.1, normalization.method = "SCT", 
                                     verbose = FALSE)

normal.integrated.1 <- RunPCA(normal.integrated.1, npcs = 100, verbose = FALSE)
normal.integrated.1 <- RunUMAP(normal.integrated.1, dims = 1:100, n.neighbors = 30)
DimPlot(normal.integrated.1, group.by = c("sample"))
DimPlot(normal.integrated.1, split.by = 'sample', ncol = 2)

normal.integrated.1@assays$SCT <- NULL
normal.list.2 <- list()
normal.list.2$C123 <- normal.integrated.1
normal.list.2$C4 <- list.sample$C4
normal.list.2 <- PrepSCTIntegration(object.list = normal.list.2, anchor.features = sel.features, 
                                    verbose = FALSE)
normal.anchors.2 <- FindIntegrationAnchors(object.list = normal.list.2, normalization.method = "SCT", 
                                           anchor.features = sel.features, verbose = FALSE)
normal.integrated.2 <- IntegrateData(anchorset = normal.anchors.2, normalization.method = "SCT", 
                                     verbose = FALSE)

normal.integrated.2 <- RunPCA(normal.integrated.2, npcs = 100, verbose = FALSE)
normal.integrated.2 <- RunUMAP(normal.integrated.2, dims = 1:100, n.neighbors = 30)
DimPlot(normal.integrated.2, group.by = c("sample"))
DimPlot(normal.integrated.2, split.by = 'sample', ncol = 2)

normal.integrated.2 <- FindNeighbors(normal.integrated.2, dims = 1:100)
normal.integrated.2 <- FindClusters(normal.integrated.2, resolution = 0.5)
DimPlot(normal.integrated.2, group.by = "seurat_clusters", label = T)

file.seurat.step.C1234 <- '/home/disk/drizzle/wgk/data/seurat.step.C1234.Rdata'
saveRDS(normal.integrated.2, file = file.seurat.step.C1234)

# plot
normal.integrated <- RunPCA(normal.integrated, npcs = 100, verbose = FALSE)
normal.integrated <- RunUMAP(normal.integrated, dims = 1:100, n.neighbors = 30)
DimPlot(normal.integrated, group.by = c("sample"))

# all samples
list.sample <- PrepSCTIntegration(object.list = list.sample, anchor.features = sel.features, 
                                    verbose = FALSE)
all.anchors <- FindIntegrationAnchors(object.list = list.sample, normalization.method = "SCT", 
                                           anchor.features = sel.features, verbose = FALSE)
all.integrated <- IntegrateData(anchorset = all.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)

all.integrated <- RunPCA(all.integrated, npcs = 100, verbose = FALSE)
all.integrated <- RunUMAP(all.integrated, dims = 1:100, n.neighbors = 30)
DimPlot(all.integrated, group.by = c("sample"))
DimPlot(all.integrated, split.by = 'sample', ncol = 2)

file.seurat.all <- '/home/disk/drizzle/wgk/data/seurat.all.Rdata'
saveRDS(all.integrated, file = file.seurat.all)
