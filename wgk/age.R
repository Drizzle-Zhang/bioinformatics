library(Seurat)

file.all <- '/home/yzj/JingMA/res/Harmony/ALL/RDS/PBMC_harmony_noUNK_ORI.RDS'
seurat.all <- readRDS(file.all)

DimPlot(seurat.all, group.by = "batch")
status <- rep('0', length(seurat.all$batch))
status[seurat.all$batch %in% c('C1', 'C2', 'C3', 'C4', 'C5')] <- 'Normal'
status[seurat.all$batch %in% c('M1', 'M2', 'M3')] <- 'Microtia'
seurat.all$status <- status
DimPlot(seurat.all, group.by = "status")

seurat.all <- FindNeighbors(seurat.all, reduction = "pca", dims = 1:100)
seurat.all <- FindClusters(seurat.all, resolution = 1.5)
DimPlot(seurat.all, group.by = "seurat_clusters", label = T)

FeaturePlot(seurat.first, features = c('FRZB', 'CTGF', 
                                       'SERPINA1', "SCRG1", 'COL9A3', 'FGFBP2'), ncol = 3)

# 1
# FeaturePlot(seurat.first, features = row.names(list.marker$`1`)[1:20], ncol = 5)
FeaturePlot(seurat.first, features = c('CFD', 'APOE', 'CSF3', 
                                       'MEDAG', "VCAN", "IL33", 
                                       'TYMP', 'VEGFA', 'PTGES'), ncol = 3)

# dashijie
file.CSC.down <- '/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdult/DEG/FC1.5_adjP0.01/Chondral_stem_cell/sigDN_gene.csv'
df.CSC.down <- read.csv(file.CSC.down)
file.CSC.up <- '/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdult/DEG/FC1.5_adjP0.01/Chondral_stem_cell/sigUP_gene.csv'
df.CSC.up <- read.csv(file.CSC.up)
file.CSC.BP.down <- '/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdult/ClusterPro/FC1.5_adjP0.01/Chondral_stem_cell/BP_DN.csv'
df.CSC.BP.down <- read.csv(file.CSC.BP.down)


