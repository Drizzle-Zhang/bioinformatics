library(Seurat)
library(ggplot2)

path.data <- '/homeold/yzj/JingMA_ORI/data/'

# male C2 24
C2 <- Read10X(paste0(path.data, 'C2'))
dimnames(C2)[[2]] <- paste('C2', dimnames(C2)[[2]], sep = '_')
seurat.C2 <- CreateSeuratObject(counts = C2)
seurat.C2@meta.data$sample <- rep('C2', dim(C2)[2])
seurat.C2[["percent.mt"]] <- PercentageFeatureSet(seurat.C2, pattern = "^MT-")
seurat.C2@meta.data$gender <- rep('Male', dim(C2)[2])

# female C3 35
C3 <- Read10X(paste0(path.data, 'C3'))
dimnames(C3)[[2]] <- paste('C3', dimnames(C3)[[2]], sep = '_')
seurat.C3 <- CreateSeuratObject(counts = C3)
seurat.C3@meta.data$sample <- rep('C3', dim(C3)[2])
seurat.C3[["percent.mt"]] <- PercentageFeatureSet(seurat.C3, pattern = "^MT-")
seurat.C3@meta.data$gender <- rep('Female', dim(C3)[2])

# merge
seurat.merge <- merge(seurat.C2, seurat.C3)
# filter  
seurat.all_filter <- 
    subset(seurat.merge, 
           subset = nFeature_RNA > 800 & nFeature_RNA < 6000 & nCount_RNA > 2000 & nCount_RNA < 50000 & percent.mt < 10)

# file.all <- '/homeold/yzj/JingMA_NEW/res/Harmony/ALL/RDS/PBMC_harmony.RDS'
# seurat.all <- readRDS(file.all)
# DimPlot(seurat.all, group.by = "batch")
status <- rep('0', length(seurat.all$batch))
status[seurat.all$batch %in% c('C1', 'C2', 'C3', 'C4', 'C5')] <- 'Normal'
status[seurat.all$batch %in% c('M1', 'M2', 'M3')] <- 'Microtia'
seurat.all$status <- status
DimPlot(seurat.all, group.by = "status")


