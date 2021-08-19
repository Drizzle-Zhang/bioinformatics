# import python package: sklearn.metrics
library(reticulate)
use_python('/home/zy/tools/anaconda3/bin/python3', required = T)
# py_config()
py_module_available('sklearn')
metrics <- import('sklearn.metrics')


############# regard sc-counts data as reference
path.output <- '/mdshare/node9/zy/scRef/Benchmark/mouse_brain/'
dataset <- 'Tasic'
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
ref.mtx <- OUT$mat_exp
ref.labels <- OUT$label[,1]
ref.dataset <- 'Tasic'

############### import unlabeled data
############### Habib
path.output <- '/mdshare/node9/zy/scRef/Benchmark/mouse_brain/'
dataset <- 'Campbell'
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_Habib <- OUT$mat_exp
label_Habib <- OUT$label
exp_sc_mat <- exp_Habib
label_sc <- label_Habib

# run methods
#############################################
### scMAGIC
library(scMAGIC)
output.scMAGIC <- scMAGIC(exp_sc_mat, ref.mtx, ref.labels, num_threads = 10)
pred.scMAGIC <- output.scMAGIC$scMAGIC.tag
table(label_sc[,1], pred.scMAGIC)

path.res <- '/mdshare/node9/zy/scRef/unassign/'
file.scMAGIC <- paste0(path.res, 'scMAGIC.Rdata')
saveRDS(output.scMAGIC, file = file.scMAGIC)
output.scMAGIC <- readRDS(file.scMAGIC)


cell_id.unassigned <- row.names(output.scMAGIC[output.scMAGIC$scMAGIC.tag == 'Unassigned',])
exp.unassigned <- exp_sc_mat[, cell_id.unassigned]
label.unassigned <- label_sc[cell_id.unassigned,]
data("MCA_ref")
output.unassigned <- scMAGIC(exp.unassigned, MCA_ref,
                             type_ref = 'sum-counts', use_RUVseq = F,
                             corr_use_HVGene1 = 2000, corr_use_HVGene2 = NULL,
                             num_threads = 8)
table(label.unassigned, output.unassigned$scMAGIC.tag)
output.unassigned[output.unassigned$scMAGIC.tag %in% c('Astrocyte_Atp1b2 high(Brain)',
                                                       'Astroglial cell(Bergman glia)(Brain)',
                                                       'Smooth muscle cell_Acta2 high(Pancreas)',
                                                       'Ependymal cell(Fetal_Brain)'),'scMAGIC.tag'] <- 'Unassigned'
output.unassigned[output.unassigned$scMAGIC.tag %in%
                      c('Hypothalamic ependymal cell(Brain)'),'scMAGIC.tag'] <-
    'Hypothalamic ependymal cell'
output.unassigned[output.unassigned$scMAGIC.tag %in%
                      c('Atrial cardiomyocyte_Acta2 high(Neonatal-Heart)'),'scMAGIC.tag'] <-
    'Atrial cardiomyocyte_Acta2 high'

path.res <- '/mdshare/node9/zy/scRef/unassign/'
file.unassign <- paste0(path.res, 'unassigned.Rdata')
saveRDS(output.unassigned, file = file.unassign)
output.unassigned <- readRDS(file.unassign)

### original plot
library(Seurat)
# data preparing
seurat.unlabeled <- CreateSeuratObject(counts = exp_sc_mat)
seurat.unlabeled <- NormalizeData(seurat.unlabeled, normalization.method = "LogNormalize",
                                  scale.factor = 10000)
seurat.unlabeled <- FindVariableFeatures(seurat.unlabeled, selection.method = "vst", nfeatures = 2000)
# VariableFeatures(seurat.unlabeled)
# VariableFeaturePlot(seurat.unlabeled)
# seurat.unlabeled[["percent.mt"]] <- PercentageFeatureSet(seurat.unlabeled, pattern = "^mt-")
# seurat.unlabeled <- ScaleData(seurat.unlabeled, vars.to.regress = "percent.mt")
# seurat.unlabeled <- ScaleData(seurat.unlabeled, features = all.genes)
seurat.unlabeled <- ScaleData(seurat.unlabeled)

# add label
use.cells <- dimnames(seurat.unlabeled@assays$RNA@counts)[[2]]
names(label_sc) <- 'cell_type'
# label_sc <- cbind(label_sc, seurat.unlabeled@reductions$umap@cell.embeddings)
seurat.unlabeled@meta.data$original.label <- label_sc[use.cells, ]
# pred.scRef <- cbind(pred.scRef, seurat.unlabeled@reductions$umap@cell.embeddings)
seurat.unlabeled@meta.data$scMAGIC.tag <- pred.scMAGIC
pred.new <- rbind(output.scMAGIC[output.scMAGIC$scMAGIC.tag != 'Unassigned',], output.unassigned)
pred.unassign <- pred.new[use.cells, 'scMAGIC.tag']
# pred.unassign <- cbind(pred.unassign, seurat.unlabeled@reductions$umap@cell.embeddings)
seurat.unlabeled@meta.data$new.tag <- pred.unassign

# PCA
seurat.unlabeled <- RunPCA(seurat.unlabeled)

# cluster
# seurat.unlabeled <- FindNeighbors(seurat.unlabeled, reduction = "pca", dims = 1:75, nn.eps = 0.5)
# seurat.unlabeled <- FindClusters(seurat.unlabeled, resolution = 3, n.start = 20)

# UMAP
seurat.unlabeled <- RunUMAP(seurat.unlabeled, dims = 1:50)

library(ggplot2)
library(ggpubr)
library(scales)
library(RColorBrewer)
# figure1: ture label
plot.umap <-
    DimPlot(seurat.unlabeled, reduction = "umap",
            label = T, repel = T, label.size = 2.5,
            group.by = 'original.label') +
    scale_color_manual(values = c('#FB9A99', '#E18A00', '#BE9C00', '#00BE70', '#1F78B4', '#8CAB00',
                                  '#00C1AB', '#B2DF8A', '#00ACFC', '#33A02C', '#D575FE'),
                       breaks = c('Astrocytes', 'Endothelial cells', 'Ependymocytes', 'Mural cells',
                                  'Neurons', 'Oligodendrocytes', 'OPC', 'Pars tuberalis',
                                  'PVMs & Microglia', 'Tanycytes', 'VLMCs')) +
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_text(size = 8, face = 'bold'),
          legend.text = element_text(size = 6),
          legend.position = 'bottom') +
    guides(color = guide_legend(ncol = 3,  byrow = TRUE, reverse = F,
                                override.aes = list(size=3),
                                keywidth = 0.1, keyheight = 0.1, default.unit = 'cm'))
ggsave(filename = paste0('cluster_', ref.dataset, '_', dataset, '.png'),
       path = path.res, plot = plot.umap,
       units = 'cm', height = 10.5, width = 8.5)

# figure2: cluster label
# DimPlot(seurat.unlabeled, reduction = "umap", label = T, group.by = 'seurat_clusters')
# figure3: scRef plus label
plot.umap.scRef <-
    DimPlot(seurat.unlabeled, reduction = "umap",
            label = T, repel = T, label.size = 2.5,
            group.by = 'scMAGIC.tag') +
    scale_color_manual(values = c('#FB9A99', '#E18A00', '#00ACFC', '#1F78B4', '#8CAB00', '#00C1AB', 'gray'),
                       breaks = c('Astrocyte', 'Endothelial Cell', 'Microglia', 'Neuron',
                                  'Oligodendrocyte', 'Oligodendrocyte Precursor Cell', 'Unassigned')) +
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_text(size = 8, face = 'bold'),
          legend.text = element_text(size = 6),
          legend.position = 'bottom') +
    guides(color = guide_legend(ncol = 3,  byrow = TRUE, reverse = F,
                                override.aes = list(size=3),
                                keywidth = 0.1, keyheight = 0.1, default.unit = 'cm'))
ggsave(filename = paste0('cluster_scRef_', ref.dataset, '_', dataset, '.png'),
       path = path.res, plot = plot.umap.scRef,
       units = 'cm', height = 10, width = 8.5)

# figure3: scRef and annotate unassigned
plot.umap.scRef.unassign <-
    DimPlot(seurat.unlabeled, reduction = "umap",
            label = T, repel = T, label.size = 2.5,
            group.by = 'new.tag') +
    scale_color_manual(values = c('#FB9A99', '#E18A00', '#00ACFC', '#1F78B4', '#8CAB00', '#00C1AB',
                                  '#BE9C00', '#800080', '#A52A2A', 'gray'),
                       breaks = c('Astrocyte', 'Endothelial Cell', 'Microglia', 'Neuron',
                                  'Oligodendrocyte', 'Oligodendrocyte Precursor Cell',
                                  'Hypothalamic ependymal cell', 'Atrial cardiomyocyte_Acta2 high',
                                  'Stromal cell(Fetal_Brain)', 'Unassigned')) +
    # xlim(-18, 15) +
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_text(size = 8, face = 'bold'),
          legend.text = element_text(size = 6),
          legend.position = 'bottom') +
    guides(color = guide_legend(ncol = 2,  byrow = TRUE, reverse = F,
                                override.aes = list(size=3),
                                keywidth = 0.1, keyheight = 0.1, default.unit = 'cm'))
ggsave(filename = paste0('cluster_scRef_unassign_', ref.dataset, '_', dataset, '.png'),
       path = path.res, plot = plot.umap.scRef.unassign,
       units = 'cm', height = 11, width = 8.5)

# ggarrange(plot.umap, plot.umap.scRef, plot.umap.scRef.unassign,
#           ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
