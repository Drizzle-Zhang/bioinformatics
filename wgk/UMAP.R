setwd('/home/zy/my_git/bioinformatics/wgk')
.libPaths('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)
library(harmony)
library(ggplot2)

file.all <- '/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/PBMC_harmony.RDS'
seurat.all <- readRDS(file.all)
seurat.all.2 <- CreateSeuratObject(seurat.all@assays$RNA@counts[
    setdiff(rownames(seurat.all@assays$RNA@counts), 'CRISPLD1'), ])
seurat.all.2@meta.data <- seurat.all@meta.data

seurat.all.2 <- NormalizeData(seurat.all.2)
seurat.all.2 <- FindVariableFeatures(seurat.all.2, nfeatures = 4000)
seurat.all.2 <- ScaleData(seurat.all.2, split.by = "batch")
seurat.all.2 <- RunPCA(seurat.all.2, verbose = F, npcs = 50)
seurat.all.2 <- RunHarmony(seurat.all.2, "batch", reduction.save = "harmony")
seurat.all.2 <- RunUMAP(seurat.all.2, reduction = "harmony", 
                        dims = 1:50, n_neighbors = 30, min.dist = 0.075)
saveRDS(seurat.all.2, file = '/home/disk/drizzle/wgk/UMAP/seurat_50_30_0.075.Rdata')

DimPlot(seurat.all.2, group.by = "batch", label = T, reduction = 'umap', pt.size = 1)
DimPlot(seurat.all.2, reduction = 'harmony')

FeaturePlot(seurat.all.2, 'EGR1', reduction = 'umap')
FeaturePlot(seurat.all.2, 'COL2A1', reduction = 'umap')


path_umap <- '/home/disk/drizzle/wgk/UMAP'
if (!file.exists(path_tsne)) {
    dir.create(path_tsne)
}

Plot.UMAP <- function(seurat.all.2, path_umap, param) {
    dim <- param$dim
    n <- param$n
    dist <- param$dist
    seurat.all.2 <- RunUMAP(seurat.all.2, reduction = "harmony", 
                                 dims = 1:dim, n_neighbors = n, min.dist = dist)
    plot.umap <- DimPlot(seurat.all.2, group.by = 'batch', 
                         label = T, reduction = 'umap')
    ggsave(plot.umap, path = path_umap, 
           filename = paste0('umap_', dim, '_', n, '_', dist, '.png'),
           height = 20, width = 25, units = 'cm')
    plot.col2 <- FeaturePlot(seurat.all.2, 'COL2A1', reduction = 'umap')
    ggsave(plot.col2, path = path_umap, 
           filename = paste0('COL2A1_', dim, '_', n, '_', dist, '.png'),
           height = 20, width = 25, units = 'cm')
    plot.egr1 <- FeaturePlot(seurat.all.2, 'EGR1', reduction = 'umap')
    ggsave(plot.egr1, path = path_umap, 
           filename = paste0('EGR1_', dim, '_', n, '_', dist, '.png'),
           height = 20, width = 25, units = 'cm')
}
# Plot.UMAP(seurat.all.2, path_umap, param)

dims <- seq(10, 50, 5)
n.neigh <- seq(5, 50, 5)
dists <- c(0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5)
list.param <- list()
idx <- 0
for (dim in dims) {
    for (n.n in n.neigh) {
        for (dist in dists) {
            idx = idx + 1
            list.param[[idx]] <- data.frame(dim = dim, n = n.n, dist = dist)
            # df.param <- rbind(df.param, data.frame(dim = dim, n = n.n, dist = dist))
        }
    }
}
library(foreach)
library(doParallel)
registerDoParallel(cores = 30)
df.res <- foreach(param = list.param) %dopar% Plot.UMAP(seurat.all.2, path_umap, param)
