library(scMAGIC)
setwd("/mdshare/node9/zy/scRef/demo")
# load reference dataset
list.Ref <- readRDS('Tasic.Rdata')
ref.mtx <- list.Ref$mat_exp
ref.labels <-list.Ref$label[, 1]
# load target dataset
list.demo <- readRDS('Campbell_2k.Rdata')
exp_sc_mat <- list.demo$mat_exp
label_sc <-list.demo$label

output.scMAGIC <- scMAGIC(exp_sc_mat, ref.mtx, ref.labels, num_threads = 4)
pred.scMAGIC <- output.scMAGIC$scMAGIC.tag
table(label_sc, pred.scMAGIC)

# Visualization
library(Seurat)
Obj.seurat <- CreateSeuratObject(counts = exp_sc_mat)
Obj.seurat <- NormalizeData(Obj.seurat)
Obj.seurat <- FindVariableFeatures(Obj.seurat, nfeatures = 2000)
Obj.seurat <- ScaleData(Obj.seurat)
Obj.seurat <- RunPCA(Obj.seurat)
Obj.seurat <- RunUMAP(Obj.seurat, dims = 1:50)
Obj.seurat@meta.data$original.label <- label_sc
Obj.seurat@meta.data$pred.tag <- pred.scMAGIC

library(ggplot2)
DimPlot(Obj.seurat, reduction = "umap",
        label = T, repel = T, group.by = 'original.label') +
    labs(title = 'True labels') + theme(plot.title = element_text(hjust = 0.5))
DimPlot(Obj.seurat, reduction = "umap",
        label = T, repel = T, group.by = 'pred.tag') +
    labs(title = 'Prediction labels') + theme(plot.title = element_text(hjust = 0.5))

# Campbell
# load target dataset
list.target <- readRDS('Campbell.Rdata')
exp_sc_mat <- list.target$mat_exp
label_sc <-list.target$label[,1]

output.scMAGIC <- scMAGIC(exp_sc_mat, ref.mtx, ref.labels, num_threads = 10)
pred.scMAGIC <- output.scMAGIC$scMAGIC.tag
table(label_sc, pred.scMAGIC)


# unassigned cells
cell_id.unassigned <- row.names(output.scMAGIC[output.scMAGIC$scMAGIC.tag == 'Unassigned',])
exp.unassigned <- exp_sc_mat[, cell_id.unassigned]
label.unassigned <- list.target$label[cell_id.unassigned,]
data("MCA_ref")
output.unassigned <- scMAGIC(exp.unassigned, MCA_ref,
                             type_ref = 'sum-counts', use_RUVseq = F,
                             corr_use_HVGene1 = 2000, corr_use_HVGene2 = NULL,
                             num_threads = 8)
table(label.unassigned, output.unassigned$scMAGIC.tag)

# reference-free annotation
library(scMAGIC)
# load target dataset
list.target <- readRDS('Campbell.Rdata')
exp_sc_mat <- list.target$mat_exp
label_sc <-list.target$label[,1]
# load MCA
data("MCA_ref")
# run scMAGIC
output.scMAGIC <- scMAGIC(exp_sc_mat, MCA_ref,
                          type_ref = 'sum-counts', use_RUVseq = F,
                          num_threads = 10)
pred.scMAGIC <- output.scMAGIC$scMAGIC.tag
# classification results
table(true.tags, pred.scMAGIC)
