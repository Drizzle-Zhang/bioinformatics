setwd('/home/zy/my_git/bioinformatics/wgk')
.libPaths('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)

path.scenic <- '/home/yzj/JingMA_NEW/res/SCENIC_All_2/'
# path.scenic <- '/home/yzj/JingMA_NEW/res/SCENIC_All/'
path.env <- paste0(path.scenic, 'int/')
if (!dir.exists(path.scenic)) {
    dir.create(path.scenic)
}
if (!dir.exists(path.env)) {
    dir.create(path.env)
}
setwd(path.scenic)

# read seurat
path.data <- '/home/disk/drizzle/wgk/data/AllSample_2_merge/'
file.all <- paste0(path.data, 'seurat_celltype.Rdata')
seurat.all <- readRDS(file.all)
seurat.chon <- subset(seurat.all, subset = celltype %in% 
                        c('Chondrocyte1', 'Chondrocyte2', 'Chondral stem cell', 
                          'Transitional chondrocyte') & batch %in% c('C4', 'C6', 'M1', 'M2', 'M3'))

# cell info
cellInfo <- data.frame(seurat.chon@meta.data)
cellInfo <- cellInfo[, c("batch", "type", "celltype")]
colnames(cellInfo) <- c("sample", "status", "celltype")
file.cellinfo <- paste0(path.env, 'cellInfo.Rds')
saveRDS(cellInfo, file=file.cellinfo)

# expression matrix
subcell <- sample(colnames(seurat.chon), 2000)
scRNAsub <- seurat.chon[,subcell]
file.sub2000 <- paste0(path.scenic, 'scRNAsub.rds')
saveRDS(scRNAsub, file.sub2000)
exprMat <- as.matrix(scRNAsub@assays$RNA@counts)

# environment
mydbDIR <- "/home/yzj/publicData/SCENIC"
mydbs <- c("hg19-500bp-upstream-7species.mc9nr.feather",
           "hg19-tss-centered-10kb-7species.mc9nr.feather")
names(mydbs) <- c("500bp", "10kb")
scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=20,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "HNSCC")
file.options <- paste0(path.env, 'scenicOptions.rds')
saveRDS(scenicOptions, file.options)

##==转录调控网络推断==##
##基因过滤
#过滤标准是基因表达量之和>细胞数*3%，且在1%的细胞中表达
genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]
##计算相关性矩阵
runCorrelation(exprMat_filtered, scenicOptions)
##TF-Targets相关性回归分析
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions, nParts = 5)
#这一步消耗的计算资源非常大，个人电脑需要几个小时的运行时间

##推断共表达模块
runSCENIC_1_coexNetwork2modules(scenicOptions)
##推断转录调控网络（regulon）
scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=5,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "HNSCC")
runSCENIC_2_createRegulons(scenicOptions)
#以上代码可增加参数coexMethod=c("w001", "w005", "top50", "top5perTarget", "top10perTarget", "top50perTarget"))
#默认6种方法的共表达网络都计算，可以少选几种方法以减少计算量

##==regulon活性评分与可视化==##
##regulons计算AUC值并进行下游分析
exprMat_all <- as.matrix(seurat.chon@assays$RNA@counts)
exprMat_all <- log2(exprMat_all+1)
scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=1,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "HNSCC")
runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_all)

runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_all)

# view results
RegulonInfo <- readRDS('int/2.5_regulonTargetsInfo.Rds')
Regulons <- readRDS('int/3.1_regulons_forAUCell.Rds')

##导入原始regulonAUC矩阵
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
# HighConf <- lapply(strsplit(RegulonName_AUC, split='_'), 
#                    function(sub_list) {
#                        if (sub_list[2] != 'extended') 
#                         {return(paste(sub_list[1], sub_list[2], sep = '_'))}
#                     })
# HighConf <- unlist(HighConf)
# seurat.chon@assays$AUC <- CreateAssayObject(data = t(as.matrix(AUCmatrix[,HighConf])))
seurat.chon@assays$AUC <- CreateAssayObject(data = t(as.matrix(AUCmatrix))[, rownames(seurat.chon@meta.data)])

# diff across cell type
clusters <- unique(seurat.chon$celltype)
list.marker <- list()
for (cluster in clusters) {
    sub.markers <- FindMarkers(seurat.chon, assay = 'AUC',
                               ident.1 = cluster, group.by = 'celltype',
                               logfc.threshold = 0, min.diff.pct = 0, only.pos = T)
    list.marker[[cluster]] <- sub.markers
}

celltypes <- unique(seurat.chon$celltype)
list.diff <- list()
for (cell in celltypes) {
    sub.seurat <- subset(seurat.chon, subset = celltype == cell)
    sub.markers <- FindMarkers(sub.seurat, assay = 'AUC', 
                               ident.1 = 'Microtia', group.by = 'type',
                               logfc.threshold = 0, min.diff.pct = 0)
    sub.markers <- sub.markers[sub.markers$p_val_adj < 0.01,]
    list.diff[[as.character(cell)]] <- sub.markers
}

# diff TF
path.M123 <- '/home/disk/drizzle/wgk/microtia_child_M1M2M3/'
file.gsea.marker <- paste0(path.M123, 'marker_GSEA.Rdata')
list.marker.gsea <- readRDS(file.gsea.marker)
TFs <- unique(RegulonInfo$TF)
View(list.marker.gsea$`Chondral stem cell`[TFs, ])





