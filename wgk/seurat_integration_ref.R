library(Seurat)
library(ggplot2)
library(patchwork)
options(future.globals.maxSize = 7000 * 1024^2)
library(future)
plan("multiprocess", workers = 6)

file.ear <- '/home/disk/drizzle/wgk/data/seurat_all_filtered.Rdata'
seurat.ear <- readRDS(file.ear)


# all
list.sample <- SplitObject(seurat.ear, split.by = 'sample')
for (i in 1:length(list.sample)) {
    list.sample[[i]] <- SCTransform(list.sample[[i]], verbose = FALSE)
}
sel.features <- SelectIntegrationFeatures(object.list = list.sample, nfeatures = 5000)

list.sct <- list()
list.sct$list.sample <- list.sample
list.sct$sel.features <- sel.features
file.sct <- '/home/disk/drizzle/wgk/data/SCT_all.Rdata'
saveRDS(list.sct, file = file.sct)
list.sct <- readRDS(file.sct)
list.sample <- list.sct$list.sample
sel.features <- list.sct$sel.features

# first
list.sample <- SplitObject(seurat.ear, split.by = 'sample')
list.sample <- list.sample[c('C1', 'C2', 'C3', 'M1', 'M2')]
for (i in 1:length(list.sample)) {
    list.sample[[i]] <- SCTransform(list.sample[[i]], verbose = FALSE)
}
sel.features <- SelectIntegrationFeatures(object.list = list.sample, nfeatures = 5000)

list.sct <- list()
list.sct$list.sample <- list.sample
list.sct$sel.features <- sel.features
file.sct <- '/home/disk/drizzle/wgk/data/ear.sct.first.Rdata'
saveRDS(list.sct, file = file.sct)
list.sct <- readRDS(file.sct)
list.sample <- list.sct$list.sample
sel.features <- list.sct$sel.features

# ref-based integration
list.sample <- PrepSCTIntegration(object.list = list.sample, anchor.features = sel.features, 
                                    verbose = FALSE)
reference_dataset <- which(names(list.sample) %in% c('C1', 'C2', 'C3'))
ref.anchors <- FindIntegrationAnchors(object.list = list.sample, normalization.method = "SCT", 
                                      anchor.features = sel.features, reference = reference_dataset)
ref.integrated <- IntegrateData(anchorset = ref.anchors, normalization.method = "SCT")

ref.integrated <- RunPCA(ref.integrated, npcs = 100, verbose = FALSE)
ref.integrated <- RunUMAP(ref.integrated, dims = 1:100, n.neighbors = 60)
DimPlot(ref.integrated, group.by = c("sample"))
DimPlot(ref.integrated, split.by = 'sample', ncol = 3)
FeaturePlot(ref.integrated, features = c('COL1A1', 'COL1A2', 'PDGFRB', 'COL5A2'))
FeaturePlot(ref.integrated, features = c('COL2A1', 'COL9A2', 'ACAN', 'SOX9'))
FeaturePlot(ref.integrated, features = c('CDK18', 'TOP2A', 'AURKB'))


file.ref.C123 <- '/home/disk/drizzle/wgk/data/seurat.ref.C123.first.Rdata'
saveRDS(ref.integrated, file = file.ref.C123)


# all
list.sample <- SplitObject(seurat.ear, split.by = 'sample')
for (i in 1:length(list.sample)) {
    list.sample[[i]] <- SCTransform(list.sample[[i]], verbose = FALSE)
}
sel.features <- SelectIntegrationFeatures(object.list = list.sample, nfeatures = 5000)

list.sct <- list()
list.sct$list.sample <- list.sample
list.sct$sel.features <- sel.features
file.sct <- '/home/disk/drizzle/wgk/data/ear.sct.first.Rdata'
saveRDS(list.sct, file = file.sct)
list.sct <- readRDS(file.sct)
list.sample <- list.sct$list.sample
sel.features <- list.sct$sel.features

# reference based
# C1, C2, C3, M1, M2
file.sct <- '/home/disk/drizzle/wgk/data/ear.sct.Rdata'
list.sct <- readRDS(file.sct)
list.sample <- list.sct$list.sample
sel.features <- list.sct$sel.features

list.sample <- PrepSCTIntegration(object.list = list.sample, anchor.features = sel.features, 
                                  verbose = FALSE)
reference_dataset <- which(names(list.sample) %in% c('C1', 'C2', 'C3', 'M1', 'M2'))
ref.anchors <- FindIntegrationAnchors(object.list = list.sample, normalization.method = "SCT", 
                                      anchor.features = sel.features, reference = reference_dataset)
ref.integrated <- IntegrateData(anchorset = ref.anchors, normalization.method = "SCT")

ref.integrated <- RunPCA(ref.integrated, npcs = 100, verbose = FALSE)
ref.integrated <- RunUMAP(ref.integrated, dims = 1:50, n.neighbors = 35)
DimPlot(ref.integrated, group.by = c("sample"))
DimPlot(ref.integrated, split.by = 'sample', ncol = 3)
FeaturePlot(ref.integrated, features = c('COL1A1', 'COL1A2', 'PDGFRB', 'COL5A2'))
FeaturePlot(ref.integrated, features = c('COL2A1', 'COL9A2', 'ACAN', 'SOX9'))
FeaturePlot(ref.integrated, features = c('CDK18', 'TOP2A', 'AURKB'))


file.seurat.all <- '/home/disk/drizzle/wgk/data/seurat.all.Rdata'
saveRDS(all.integrated, file = file.seurat.all)
