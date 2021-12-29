library(Seurat)

PBMC_LIGER <- readRDS('/home/yzj/JingMA_Liger3/res/Liger/ALL/RDS/PBMC_LIGER.RDS')
ref.labels <- PBMC_LIGER$CellType
ref.mtx <- PBMC_LIGER@assays$RNA@counts

path.data <- '/home/yzj/JingMA/data/'
C4 <- Read10X(paste0(path.data, 'C4'))
dimnames(C4)[[2]] <- paste('C4', dimnames(C4)[[2]], sep = '_')
seurat.C4 <- CreateSeuratObject(counts = C4)
seurat.C4@meta.data$sample <- rep('C4', dim(C4)[2])
seurat.C4[["percent.mt"]] <- PercentageFeatureSet(seurat.C4, pattern = "^MT-")
seurat.C4_filter <- subset(seurat.C4, subset = nFeature_RNA > 800 & nFeature_RNA < 6000 & 
                               nCount_RNA > 3000 & nCount_RNA < 50000 & percent.mt < 15)
exp_sc_mat <- seurat.C4_filter@assays$RNA@counts

source('/home/zy/my_git/scRef/main/scRef.v20.R')
setwd('~/my_git/scRef')
result.scref <- SCREF(exp_sc_mat, ref.mtx, ref.labels,
                      type_ref = 'sc-counts', use.RUVseq = T, out.group = 'HCA',
                      cluster.speed = F, cluster.cell = 20,
                      method1 = 'spearman', 
                      # corr_use_HVGene1 = 5000, corr_use_HVGene2 = 5000,
                      GMM.ceiling_cutoff = 20,
                      min_cell = 1, CPU = 8)
pred.scRef <- result.scref$final.out
file.label.1 <- '/home/zy/scRef/ear/scMAGIC.C4.Rdata'
saveRDS(pred.scRef, file = file.label.1)

result.scref.2 <- SCREF(exp_sc_mat, ref.mtx, ref.labels,
                        type_ref = 'sc-counts', identify_unassigned = F, single_round = F,
                        method1 = 'spearman', 
                        min_cell = 1, CPU = 8)
pred.scRef.2 <- result.scref.2$final.out
file.label.2 <- '/home/zy/scRef/ear/scMAGIC.C4.2.Rdata'
saveRDS(pred.scRef.2, file = file.label.2)

result.scref.3 <- SCREF(exp_sc_mat, ref.mtx, ref.labels,
                        type_ref = 'sc-counts', identify_unassigned = F, single_round = T,
                        method1 = 'spearman', 
                        min_cell = 1, CPU = 8)
pred.scRef.3 <- result.scref.3$final.out

### sciBet
suppressMessages(library(tidyverse))
suppressMessages(library(scibet))
suppressMessages(library(viridis))
suppressMessages(library(ggsci))
train_set <- as.data.frame(t(ref.mtx))
train_set$label <- ref.labels
test_set <- as.data.frame(t(exp_sc_mat))
sciBet <- SciBet(train_set, test_set)

## plot
seurat.unlabeled <- seurat.ear
seurat.unlabeled@meta.data$scRef.tag <- pred.scRef.3
seurat.unlabeled@meta.data$sciBet <- sciBet

DimPlot(seurat.unlabeled, reduction = "umap", label = T, repel = T, group.by = 'scRef.tag')

DimPlot(seurat.unlabeled, reduction = "umap", label = T, repel = T, group.by = 'sciBet')

plot.umap.scRef <- 
    DimPlot(seurat.unlabeled, reduction = "umap", label = T, repel = T, group.by = 'scRef.tag') + 
    scale_color_manual(values = c(hue_pal()(10), 'gray'),
                       breaks = c(names(table(pred.scRef)))) + 
    theme_bw() + 
    theme(axis.text = element_text(size = 9),
          panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 11))


