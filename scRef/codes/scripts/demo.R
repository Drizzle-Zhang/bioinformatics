library(reticulate)
use_python('/home/zy/tools/anaconda3/bin/python3', required = T)
py_config()

# reference
path.output <- '/mdshare/node9/zy/scRef/Benchmark/mouse_brain/'
list.Ref <- readRDS(paste0(path.output, ref.dataset, '.Rdata'))
exp_Tasic <- list.Ref$mat_exp
label_Tasic <- list.Ref$label
ref.labels <-label_Tasic[, 1]
ref.mtx <- exp_Tasic
ref.dataset <- 'Tasic'

# sample test dataset
dataset <- 'Campbell'
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
# exp_sc_mat <- OUT$mat_exp
# label_sc <- OUT$label
exp_Habib <- OUT$mat_exp
label_Habib <- OUT$label
set.seed(1234)
cells.sample <- sample(rownames(label_Habib),2000)
exp_sc_mat <- exp_Habib[,cells.sample]
label_sc <- label_Habib[cells.sample,]
list.demo <- list()
list.demo$mat_exp <- exp_sc_mat
list.demo$label <- label_sc
saveRDS(list.demo, file = '/home/zy/my_git/scMAGIC_scripts/data/Campbell_2k.Rdata')

# load data
path.output <- '/mdshare/zy/scRef/Benchmark/mouse_brain/'
ref.dataset <- 'Tasic'
list.Ref <- readRDS(paste0(path.output, ref.dataset, '.Rdata'))
ref.mtx <- list.Ref$mat_exp
ref.labels <-list.Ref$label[, 1]

list.demo <- readRDS('/local/zy/my_git/scMAGIC_scripts/data/Campbell_2k.Rdata')
exp_sc_mat <- list.demo$mat_exp
label_sc <-list.demo$label

source('/local/zy/my_git/scRef/main/scRef.v21.R')
setwd('~/my_git/scRef')
output.scMAGIC <- SCREF(exp_sc_mat, ref.mtx, ref.labels, CPU = 4)
print(output.scMAGIC$run.time)
pred.scMAGIC <- output.scMAGIC$final.out$scRef.tag
table(label_sc, pred.scMAGIC)

library(Seurat)
seurat.unlabeled <- CreateSeuratObject(counts = exp_sc_mat)
seurat.unlabeled <- NormalizeData(seurat.unlabeled)
seurat.unlabeled <- FindVariableFeatures(seurat.unlabeled, nfeatures = 2000)
seurat.unlabeled <- ScaleData(seurat.unlabeled)
seurat.unlabeled <- RunPCA(seurat.unlabeled)
seurat.unlabeled <- RunUMAP(seurat.unlabeled, dims = 1:50)
seurat.unlabeled@meta.data$original.label <- label_sc
seurat.unlabeled@meta.data$pred.tag <- pred.scMAGIC


library(ggplot2)
path.res <- "/mdshare/node9/zy/scRef/demo"
# origin label
plot.umap <-
    DimPlot(seurat.unlabeled, reduction = "umap",
            label = T, repel = T, group.by = 'original.label') +
    labs(title = 'True labels') + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = 'TrueLabels.png',
       path = path.res, plot = plot.umap,
       units = 'cm', height = 15, width = 20)
# pred label
plot.umap.pred <-
    DimPlot(seurat.unlabeled, reduction = "umap",
            label = T, repel = T, group.by = 'pred.tag') +
    labs(title = 'Prediction labels') + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = 'PredLabels.png',
       path = path.res, plot = plot.umap.pred,
       units = 'cm', height = 15, width = 22)


library(ggplot2)
plot.umap <-
    DimPlot(seurat.unlabeled, reduction = "umap",
            label = T, repel = T, label.size = 2.5,
            group.by = 'original.label') +
    scale_color_manual(values = c('#24B700', '#E18A00', '#BE9C00', '#00BE70', '#24B700', '#8CAB00',
                                  '#00C1AB', '#00BBDA', '#00ACFC', '#8B93FF', '#D575FE'),
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


plot.umap.pred <-
    DimPlot(seurat.unlabeled, reduction = "umap",
            label = T, repel = T, label.size = 2.5,
            group.by = 'pred.tag') +
    scale_color_manual(values = c('#24B700', '#E18A00', '#00ACFC', '#24B700', '#8CAB00', '#00C1AB', 'gray'),
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


# altas
path.output <- '/mdshare/node9/zy/scRef/atlas_anno/'
dataset <- 'Tasic2018'
OUT <- readRDS(paste0('/mdshare/node9/zy/scRef/Benchmark/mouse_brain/', dataset, '.Rdata'))
exp_sc_mat <- OUT$mat_exp
label_sc <- OUT$label
set.seed(123)
sel.sample <- rownames(label_sc)[label_sc[,1] %in%
                                     c('Endothelial_Endo', 'Endothelial_Peri', 'Endothelial_SMC',
                                       'Non-Neuronal_Astro', 'Non-Neuronal_Macrophage',
                                       'Non-Neuronal_Oligo', 'Non-Neuronal_VLMC')]
col.neu <- setdiff(colnames(exp_sc_mat), sel.sample)
sel.sample <- c(sel.sample, sample(col.neu,3500-length(sel.sample)))
exp_input <- exp_sc_mat[, sel.sample]
sel_label <- label_sc[sel.sample, ]

library(scMAGIC)
data("MCA_ref")
output.scMAGIC <- scMAGIC(exp_input, MCA_ref,
                          type_ref = 'sum-counts', use_RUVseq = F,
                          min_cell = 5, num_threads = 10)
pred.scMAGIC <- output.scMAGIC$scMAGIC.tag
true.tags <- sel_label
true.tags[true.tags %in% c('Endothelial_Endo')] <- 'Endothelial cell'
true.tags[true.tags %in% c('Endothelial_Peri')] <- 'Pericyte'
true.tags[true.tags %in% c('Endothelial_SMC')] <- 'Smooth muscle cell'
true.tags[true.tags %in% c('GABAergic_', 'GABAergic_Lamp5', 'GABAergic_Meis2',
                           'GABAergic_Pvalb', 'GABAergic_Serpinf1', 'GABAergic_Sncg',
                           'GABAergic_Sst', 'GABAergic_Vip')] <- 'GABAergic Neuron'
true.tags[true.tags %in% c('Glutamatergic_', 'Glutamatergic_CR', 'Glutamatergic_L2/3 IT',
                           'Glutamatergic_L4', 'Glutamatergic_L5 IT', "Glutamatergic_L5 PT",
                           'Glutamatergic_L6 CT',  'Glutamatergic_L6 IT',
                           'Glutamatergic_L6b', 'Glutamatergic_NP')] <- 'Glutamatergic Neuron'
true.tags[true.tags %in% c('Non-Neuronal_Astro')] <- 'Astrocyte'
true.tags[true.tags %in% c('Non-Neuronal_Macrophage')] <- 'PVM & Microglia'
true.tags[true.tags %in% c('Non-Neuronal_Oligo')] <- 'Oligodendrocyte & OPC'
true.tags[true.tags %in% c('Non-Neuronal_VLMC')] <- 'VLMC'
table(true.tags, pred.scMAGIC)


list.demo.atlas <- list()
list.demo.atlas$mat_exp <- exp_input
list.demo.atlas$label <- true.tags
saveRDS(list.demo.atlas, file = '/local/zy/my_git/scMAGIC_scripts/data/MouseNeocortex.Rdata')

library(scMAGIC)
# load target dataset
list.target <- readRDS('/local/zy/my_git/scMAGIC_scripts/data/MouseNeocortex.Rdata')
exp_sc_mat <- list.target$mat_exp
label_sc <-list.target$label
# load MCA
data("MCA_ref")
# run scMAGIC
output.scMAGIC <- scMAGIC(exp_sc_mat, MCA_ref,
                          type_ref = 'sum-counts', use_RUVseq = F,
                          min_cell = 5, num_threads = 10)
pred.scMAGIC <- output.scMAGIC$scMAGIC.tag
# classification results
table(label_sc, pred.scMAGIC)



# 5-fold CV
ref.dataset <- 'Tasic'
path.output <- '/mdshare/node9/zy/scRef/Benchmark/mouse_brain/'
list.Ref <- readRDS(paste0(path.output, ref.dataset, '.Rdata'))
exp_Tasic <- list.Ref$mat_exp
label_Tasic <- list.Ref$label

# ref
ref.col <- sample(colnames(exp_Tasic), floor(ncol(exp_Tasic)*0.8))
test.col <- setdiff(colnames(exp_Tasic), ref.col)
ref.mtx <- exp_Tasic[,ref.col]
ref.labels <- label_Tasic[ref.col, 1]
test.mat <- exp_Tasic[,test.col]
test.label <- label_Tasic[test.col, 1]

output.scMAGIC <- scMAGIC(test.mat, ref.mtx, ref.labels, num_threads = 4)
pred.scMAGIC <- output.scMAGIC$scMAGIC.tag
table(test.label, pred.scMAGIC)


