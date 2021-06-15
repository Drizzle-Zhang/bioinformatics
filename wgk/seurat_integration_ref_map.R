library(uwot)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
options(future.globals.maxSize = 7000 * 1024^2)
# library(reticulate)
# use_python('/home/zy/tools/anaconda3/bin/python3', required = T)
# library(future)
# plan("multiprocess", workers = 6)

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
# saveRDS(list.sct, file = file.sct)
list.sct <- readRDS(file.sct)
list.sample <- list.sct$list.sample
sel.features <- list.sct$sel.features

# first
seurat.first <- subset(seurat.ear, subset = sample %in% c('C1', 'C2', 'C3', 'M1', 'M2'))
seurat.first <- NormalizeData(seurat.first)
seurat.first <- FindVariableFeatures(seurat.first, nfeatures = 5000)
seurat.first <- ScaleData(seurat.first, split.by = "sample")
seurat.first <- RunPCA(seurat.first, verbose = F, npcs = 100)
seurat.first <- RunUMAP(seurat.first, dims = 1:100, n.neighbors = 50)
seurat.first <- FindNeighbors(seurat.first, reduction = "pca", dims = 1:100)
seurat.first <- FindClusters(seurat.first, resolution = 1.2)

DimPlot(seurat.first, group.by = "sample")
DimPlot(seurat.first, group.by = "seurat_clusters", label = T)
DimPlot(seurat.first, group.by = "status", label = T)

file.seurat <- '/home/disk/drizzle/wgk/data/marker_1.5_merge/seurat_first.Rdata'
seurat.first <- readRDS(file.seurat)
seurat.first <- RunPCA(seurat.first, verbose = F, npcs = 200)
seurat.first <- RunUMAP(seurat.first, dims = 1:75, n.neighbors = 50, umap.method = "uwot")

seurat.C4 <- list.sample$`C4`
seurat.C4 <- NormalizeData(seurat.C4)
seurat.C4 <- FindVariableFeatures(seurat.C4, nfeatures = 5000)
seurat.C4 <- ScaleData(seurat.C4)
reference <- seurat.first
reference@reductions$`umot-learn` <- seurat.first@reductions$umap
DimPlot(reference, group.by = "sample")

overlap.features <- intersect(VariableFeatures(seurat.C4), VariableFeatures(seurat.first))
seurat.C4 <- subset(list.sample$`C4`, features = overlap.features)
reference <- subset(reference, features = overlap.features)
anchors <- FindTransferAnchors(
    reference = reference,
    query = seurat.C4,
    reference.assay = 'RNA', query.assay = 'RNA',
    features = overlap.features,
    reference.reduction = 'pca',
    npcs = NULL, project.query = F,
    normalization.method = "LogNormalize",
    dims = 1:100
)
seurat.C4 <- MapQuery(
    anchorset = anchors,
    query = seurat.C4,
    reference = reference,
    refdata = list(celltype = 'celltype'),
    reference.reduction = 'pca',
    reduction.model = "uwot"
)
DimPlot(seurat.C4, reduction = "ref.umap", group.by = 'celltype')

# try
reference <- LoadH5Seurat("/home/disk/drizzle/wgk/data/pbmc_multimodal.h5seurat")


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
