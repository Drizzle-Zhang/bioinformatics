###################
## 1. Seurat preprocessing: Normalization, scale and PCA
###################
library(Seurat)
library(dplyr)
library(ggplot2)
pbmc_batch <- readRDS('/home/yzj/JingMA_NEW/res/QC/ALL/RDS/PBMC_QC.RDS')

####
pbmc_batch <- NormalizeData(object = pbmc_batch, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_batch <- FindVariableFeatures(object = pbmc_batch, selection.method = "vst", nfeatures = 4000)
pbmc_batch <- ScaleData(object = pbmc_batch, features = VariableFeatures(object = pbmc_batch))
pbmc_batch <- RunPCA(object = pbmc_batch, seed.use=123, npcs=150,
                     features = VariableFeatures(object = pbmc_batch), ndims.print=1,nfeatures.print=1)
pbmc_batch <- RunTSNE(pbmc_batch, dims = 1:50, seed.use = 123,n.components=2)
pbmc_batch <- RunUMAP(pbmc_batch, dims = 1:50, seed.use = 123,n.components=2)

########
saveRDS(pbmc_batch,'/home/yzj/JingMA_NEW/res/QC/ALL/RDS/PBMC_ORI.RDS')
########


##################
# 2. Harmony
##################
library(Seurat)
library(dplyr)
library(ggplot2)
library(SeuratWrappers)
library(harmony)

pbmc_batch <- readRDS(file='/home/yzj/JingMA_NEW/res/QC/ALL/RDS/PBMC_ORI.RDS')
system('mkdir -p /home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/')

theta=2
block.size = 0.05
pbmc <- RunHarmony(pbmc_batch,group.by.vars = 'batch',theta = theta,block.size = block.size)
dimN=15
neighborN=20
mindistN=0.05
pbmc <- RunTSNE(pbmc, dims = 1:dimN,reduction = "harmony")
pbmc <- RunUMAP(pbmc,dims=1:dimN, reduction = "harmony", seed.use = 123,n.components=2,
                     n.neighbors = neighborN,min.dist = mindistN)
saveRDS(pbmc,'JingMA_NEW/res/Harmony/ALL/RDS/PBMC_harmony.RDS')
