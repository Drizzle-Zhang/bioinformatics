library(Seurat)
library(ggplot2)

file.ear <- '/home/disk/drizzle/wgk/data/ear.filtered.Rdata'
seurat.ear <- readRDS(file.ear)

list.sample <- SplitObject(seurat.ear, split.by = 'sample')
for (i in 1:length(list.sample)) {
    list.sample[[i]] <- NormalizeData(list.sample[[i]], verbose = FALSE)
    list.sample[[i]] <- FindVariableFeatures(list.sample[[i]], selection.method = "vst", 
                                               nfeatures = 4000, verbose = FALSE)
}
file.standard <- '/home/disk/drizzle/wgk/data/ear.standard.Rdata'
saveRDS(list.sample, file = file.standard)
list.sample <- readRDS(file.standard)

features <- SelectIntegrationFeatures(object.list = list.sample, nfeatures = 4000)
list.sample <- lapply(X = list.sample, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, npcs = 100, verbose = FALSE)
})

# normal
normal.list <- list.sample[c('C1', 'C2', 'C3', 'C4')]
normal.anchors <- FindIntegrationAnchors(object.list = normal.list, anchor.features = features,
                                         dims = 1:50, reduction = "rpca")
normal.integrated <- IntegrateData(anchorset = normal.anchors, dims = 1:50)
DefaultAssay(normal.integrated) <- "integrated"
normal.integrated <- ScaleData(normal.integrated, verbose = FALSE)
normal.integrated <- RunPCA(normal.integrated, npcs = 50, verbose = FALSE)
normal.integrated <- RunUMAP(normal.integrated, reduction = "pca", dims = 1:50)
# plot
normal.integrated <- RunPCA(normal.integrated, npcs = 100, verbose = FALSE)
normal.integrated <- RunUMAP(normal.integrated, dims = 1:50, n.neighbors = 30, min.dist = 0.05)
DimPlot(normal.integrated, group.by = c("sample"))
DimPlot(normal.integrated, split.by = 'sample', ncol = 2)



