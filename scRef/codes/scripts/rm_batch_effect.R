# import data
library(Seurat)
library(SeuratData)
data("panc8")

# delete some cells 
use.cells <- c('acinar', 'activated_stellate', 'alpha', 'beta', 'delta', 'ductal', 
               'endothelial', 'gamma', 'quiescent_stellate')
exp_panc <- as.matrix(panc8@assays$RNA@counts)
label_panc <- as.character(panc8$celltype)
platform <- as.character(panc8$dataset)
exp_panc <- exp_panc[, label_panc %in% use.cells]
platform <- platform[label_panc %in% use.cells]
label_panc <- label_panc[label_panc %in% use.cells]

# reference
ref.dataset <- 'celseq2'
ref.mtx <- exp_panc[, platform %in% c('celseq2')]
ref.labels <- label_panc[platform %in% c('celseq2')]

source('/home/zy/my_git/scRef/main/scRef.v20.R')
setwd('~/my_git/scRef')

# query data
dataset <- 'indrop'
exp_indrop <- exp_panc[, platform %in% c('indrop1', 'indrop2', 'indrop3', 'indrop4')]
label_indrop <- label_panc[platform %in% c('indrop1', 'indrop2', 'indrop3', 'indrop4')]

result.scref_indrop <- SCREF(exp_indrop, ref.mtx, ref.labels,
                      type_ref = 'sc-counts', use.RUVseq = T, 
                      out.group = 'HCA',
                      cluster.speed = F, cluster.cell = 3,
                      GMM.floor_cutoff = 3, GMM.ceiling_cutoff = 20,
                      threshold.recall = 0.5,
                      min_cell = 1, CPU = 8)
pred.scRef_indrop <- result.scref_indrop$final.out$scRef.tag
table(label_indrop, pred.scRef_indrop)
pred.scRef_indrop[pred.scRef_indrop == 'Unassigned'] <- 'Unassigned_indrop'

dataset <- 'smartseq2'
exp_smartseq2 <- exp_panc[, platform %in% c('smartseq2')]
label_smartseq2 <- label_panc[platform %in% c('smartseq2')]

result.scref_smartseq2 <- SCREF(exp_smartseq2, ref.mtx, ref.labels,
                      type_ref = 'sc-counts', use.RUVseq = T, 
                      out.group = 'HCA',
                      cluster.speed = F, cluster.cell = 5,
                      GMM.floor_cutoff = 3, GMM.ceiling_cutoff = 20,
                      threshold.recall = 0.5,
                      min_cell = 1, CPU = 8)
pred.scRef_smartseq2 <- result.scref_smartseq2$final.out$scRef.tag
table(label_smartseq2, pred.scRef_smartseq2)
pred.scRef_smartseq2[pred.scRef_smartseq2 == 'Unassigned'] <- 'Unassigned_smartseq2'

dataset <- 'celseq'
exp_celseq <- exp_panc[, platform %in% c('celseq')]
label_celseq <- label_panc[platform %in% c('celseq')]

result.scref_celseq <- SCREF(exp_celseq, ref.mtx, ref.labels,
                                type_ref = 'sc-counts', use.RUVseq = T, 
                                out.group = 'HCA',
                                cluster.speed = F, cluster.cell = 3,
                                GMM.floor_cutoff = 3, GMM.ceiling_cutoff = 20,
                                threshold.recall = 0.5,
                                min_cell = 1, CPU = 8)
pred.scRef_celseq <- result.scref_celseq$final.out$scRef.tag
table(label_celseq, pred.scRef_celseq)
pred.scRef_celseq[pred.scRef_celseq == 'Unassigned'] <- 'Unassigned_celseq'

dataset <- 'fluidigmc1'
exp_fluidigmc1 <- exp_panc[, platform %in% c('fluidigmc1')]
label_fluidigmc1 <- label_panc[platform %in% c('fluidigmc1')]

result.scref_fluidigmc1 <- SCREF(exp_fluidigmc1, ref.mtx, ref.labels,
                                type_ref = 'sc-counts', use.RUVseq = T, 
                                out.group = 'HCA',
                                cluster.speed = F, cluster.cell = 3,
                                GMM.floor_cutoff = 3, GMM.ceiling_cutoff = 20,
                                threshold.recall = 0.5,
                                min_cell = 1, CPU = 8)
pred.scRef_fluidigmc1 <- result.scref_fluidigmc1$final.out$scRef.tag
table(label_fluidigmc1, pred.scRef_fluidigmc1)
pred.scRef_fluidigmc1[pred.scRef_fluidigmc1 == 'Unassigned'] <- 'Unassigned_fluidigmc1'

# data preparing
exp_sc_mat <- cbind(exp_indrop, exp_smartseq2, exp_celseq, exp_fluidigmc1)
pred.label <- c(pred.scRef_indrop, pred.scRef_smartseq2, pred.scRef_celseq, pred.scRef_fluidigmc1)
true.label <- c(label_indrop, label_smartseq2, label_celseq, label_fluidigmc1)
type.platform <- c(rep('indrop', length(pred.scRef_indrop)),
                   rep('smartseq2', length(pred.scRef_smartseq2)),
                   rep('celseq', length(pred.scRef_celseq)),
                   rep('fluidigmc1', length(pred.scRef_fluidigmc1)))
seurat.unlabeled <- CreateSeuratObject(counts = exp_sc_mat)
seurat.unlabeled <- NormalizeData(seurat.unlabeled, normalization.method = "LogNormalize", 
                                  scale.factor = 10000)
seurat.unlabeled <- FindVariableFeatures(seurat.unlabeled, selection.method = "vst", nfeatures = 2000)
seurat.unlabeled <- ScaleData(seurat.unlabeled)

# add label
# use.cells <- dimnames(seurat.unlabeled@assays$RNA@counts)[[2]]
seurat.unlabeled@meta.data$original.label <- true.label
seurat.unlabeled@meta.data$scRef.tag <- pred.label
seurat.unlabeled@meta.data$platform <- type.platform

# PCA
seurat.unlabeled <- RunPCA(seurat.unlabeled, npcs = 100, verbose = F)

# UMAP
seurat.unlabeled <- RunUMAP(seurat.unlabeled, dims = 1:20, n.neighbors = 30)

# figure1: ture label
library(ggplot2)
library(scales)
plot.umap <- 
    DimPlot(seurat.unlabeled, reduction = "umap", label = T, repel = T, group.by = 'original.label') + 
    # xlim(-16, 14) +
    theme_bw() + 
    theme(axis.text = element_text(size = 9),
          panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 11))
ggsave(filename = paste0('cluster_', ref.dataset, '_', dataset, '.png'), 
       path = path.res, plot = plot.umap,
       units = 'cm', height = 18, width = 24)

# figure2: scRef plus label
plot.umap.scRef <- 
    DimPlot(seurat.unlabeled, reduction = "umap", label = T, repel = T, group.by = 'scRef.tag') + 
    scale_color_manual(values = c(hue_pal()(9), 'gray'),
                       breaks = c(names(table(pred.label)))) + 
    theme_bw() + 
    theme(axis.text = element_text(size = 9),
          panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 11))
ggsave(filename = paste0('cluster_scRef_', ref.dataset, '_', dataset, '.png'), 
       path = path.res, plot = plot.umap.scRef,
       units = 'cm', height = 18, width = 24)

# figure3: platform
plot.umap.scRef <- 
    DimPlot(seurat.unlabeled, reduction = "umap", label = T, repel = T, group.by = 'platform') + 
    # scale_color_manual(values = c(hue_pal()(10), 'gray'),
    #                    breaks = c(names(table(pred.label)))) + 
    theme_bw() + 
    theme(axis.text = element_text(size = 9),
          panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 11))
ggsave(filename = paste0('cluster_scRef_', ref.dataset, '_', dataset, '.png'), 
       path = path.res, plot = plot.umap.scRef,
       units = 'cm', height = 18, width = 24)

# combat
library(sva)
mtx.log <- as.matrix(seurat.unlabeled@assays$RNA@data)[VariableFeatures(seurat.unlabeled),]
mtx.combat <- ComBat(mtx.log, type.platform, mod = model.matrix(~as.factor(true.label)))
seurat.combat <- CreateSeuratObject(counts = mtx.combat)
seurat.combat <- NormalizeData(seurat.combat, normalization.method = "LogNormalize", 
                                  scale.factor = 10000)
seurat.combat@assays$RNA@data <- mtx.combat
seurat.combat <- FindVariableFeatures(seurat.combat, selection.method = "vst", nfeatures = 2000)
seurat.combat <- ScaleData(seurat.combat)

# add label
# use.cells <- dimnames(seurat.unlabeled@assays$RNA@counts)[[2]]
seurat.combat@meta.data$original.label <- true.label
seurat.combat@meta.data$scRef.tag <- pred.label
seurat.combat@meta.data$platform <- type.platform

# PCA
seurat.combat <- RunPCA(seurat.combat, npcs = 100, verbose = F)

# UMAP
seurat.combat <- RunUMAP(seurat.combat, dims = 1:20, n.neighbors = 30)

# figure4: ture label
library(ggplot2)
library(scales)
plot.umap <- 
    DimPlot(seurat.combat, reduction = "umap", label = T, repel = T, group.by = 'original.label') + 
    # xlim(-16, 14) +
    theme_bw() + 
    theme(axis.text = element_text(size = 9),
          panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 11))
ggsave(filename = paste0('cluster_', ref.dataset, '_', dataset, '.png'), 
       path = path.res, plot = plot.umap,
       units = 'cm', height = 18, width = 24)

# figure5: platform
plot.umap.scRef <- 
    DimPlot(seurat.combat, reduction = "umap", label = T, repel = T, group.by = 'platform') + 
    # scale_color_manual(values = c(hue_pal()(10), 'gray'),
    #                    breaks = c(names(table(pred.label)))) + 
    theme_bw() + 
    theme(axis.text = element_text(size = 9),
          panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 11))
ggsave(filename = paste0('cluster_scRef_', ref.dataset, '_', dataset, '.png'), 
       path = path.res, plot = plot.umap.scRef,
       units = 'cm', height = 18, width = 24)


