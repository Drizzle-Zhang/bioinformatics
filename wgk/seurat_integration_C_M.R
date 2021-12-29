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
normal.integrated <- RunUMAP(normal.integrated, dims = 1:50, n.neighbors = 20)
DimPlot(normal.integrated, group.by = c("sample"))
DimPlot(normal.integrated, split.by = 'sample', ncol = 2)
# save
file.normal <- '/home/disk/drizzle/wgk/data/seurat.normal.Rdata'
saveRDS(normal.integrated, file = file.normal)
normal.integrated <- readRDS(file.normal)

# microtia
microtia.list <- list.sample[c('M1', 'M2', 'M3')]
microtia.list <- PrepSCTIntegration(object.list = microtia.list, anchor.features = sel.features, 
                                    verbose = FALSE)
microtia.anchors <- FindIntegrationAnchors(object.list = microtia.list, normalization.method = "SCT", 
                                         anchor.features = sel.features, verbose = FALSE)
microtia.integrated <- IntegrateData(anchorset = microtia.anchors, normalization.method = "SCT", 
                                   verbose = FALSE)
# plot
microtia.integrated <- RunPCA(microtia.integrated, npcs = 100, verbose = FALSE)
microtia.integrated <- RunUMAP(microtia.integrated, dims = 1:100, n.neighbors = 30)
DimPlot(microtia.integrated, group.by = c("sample"))
DimPlot(microtia.integrated, split.by = 'sample', ncol = 2)
# save
file.microtia <- '/home/disk/drizzle/wgk/data/seurat.microtia.Rdata'
saveRDS(microtia.integrated, file = file.microtia)
microtia.integrated <- readRDS(file.microtia)

# merge directly
normal.integrated@assays$integrated@counts <- normal.integrated@assays$integrated@data
microtia.integrated@assays$integrated@counts <- microtia.integrated@assays$integrated@data
seurat.merge <- merge(normal.integrated, microtia.integrated)
# merge N and M
normal.integrated@assays$SCT <- NULL
normal.integrated@assays$RNA <- NULL
microtia.integrated@assays$SCT <- NULL
microtia.integrated@assays$RNA <- NULL
list.N.M <- list()
list.N.M$N <- normal.integrated
list.N.M$M <- microtia.integrated

# sel.features.2 <- SelectIntegrationFeatures(object.list = list.N.M, nfeatures = 4000)
list.N.M <- PrepSCTIntegration(object.list = list.N.M, anchor.features = sel.features, 
                                    verbose = FALSE)
merge.anchors <- FindIntegrationAnchors(object.list = list.N.M, normalization.method = "SCT", 
                                           anchor.features = sel.features, verbose = FALSE)
merge.integrated <- IntegrateData(anchorset = merge.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)

# pca
list.N.M <- PrepSCTIntegration(object.list = list.N.M, anchor.features = sel.features, 
                               verbose = FALSE)
pca.anchors <- FindIntegrationAnchors(object.list = list.N.M, dims = 1:100, reduction = 'rpca', 
                                      normalization.method = "SCT", anchor.features = sel.features)
pca.integrated <- IntegrateData(anchorset = pca.anchors, normalization.method = "SCT", dims = 1:100)
pca.integrated <- RunPCA(pca.integrated, npcs = 100, verbose = FALSE, features = sel.features)
pca.integrated <- RunUMAP(pca.integrated, dims = 1:100, n.neighbors = 30)
DimPlot(pca.integrated, group.by = c("sample"))


# plot
merge.integrated <- RunPCA(merge.integrated, npcs = 100, verbose = FALSE, features = sel.features)
merge.integrated <- RunUMAP(merge.integrated, dims = 1:100, n.neighbors = 30)
DimPlot(merge.integrated, group.by = c("sample"))
DimPlot(merge.integrated, split.by = 'sample', ncol = 4)
# save
file.merge <- '/home/disk/drizzle/wgk/data/seurat.N.M.Rdata'
saveRDS(merge.integrated, file = file.merge)


