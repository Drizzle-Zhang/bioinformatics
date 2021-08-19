# import python package: sklearn.metrics
library(reticulate)
use_python('/home/zy/tools/anaconda3/bin/python3', required = T)
# py_config()
py_module_available('sklearn')
metrics <- import('sklearn.metrics')

# function of data preparation
prepare.data <- function(file.data.unlabeled, file.label.unlabeled, 
                         del.label = c('miss')) {
    library(stringr)
    data.unlabeled <- read.delim(file.data.unlabeled, row.names=1)
    data.unlabeled <- floor(data.unlabeled)
    names(data.unlabeled) <- str_replace_all(names(data.unlabeled), '_', '.')
    names(data.unlabeled) <- str_replace_all(names(data.unlabeled), '-', '.')
    # read label file
    file.label.unlabeled <- file.label.unlabeled
    label.unlabeled <- read.delim(file.label.unlabeled, row.names=1)
    row.names(label.unlabeled) <- str_replace_all(row.names(label.unlabeled), '_', '.')
    row.names(label.unlabeled) <- str_replace_all(row.names(label.unlabeled), '-', '.')
    col.name1 <- names(data.unlabeled)[1]
    if (substring(col.name1, 1, 1) == 'X') {
        row.names(label.unlabeled) <- paste0('X', row.names(label.unlabeled))
    }
    # filter data
    use.cols <- row.names(label.unlabeled)[!label.unlabeled[,1] %in% del.label]
    data.filter <- data.unlabeled[,use.cols]
    label.filter <- data.frame(label.unlabeled[use.cols,], row.names = use.cols)
    
    OUT <- list()
    OUT$data.filter <- data.filter
    OUT$label.filter <- label.filter
    return(OUT)
    
}

# evaluation
simple.evaluation <- function(true.tag, scRef.tag, df.cell.names) {
    # uniform tags
    for (j in 1:dim(df.cell.names)[1]) {
        scRef.tag[scRef.tag == df.cell.names[j, 'ref.name']] <- 
            df.cell.names[j, 'sc.name']
    }
    
    true.tag[true.tag %in% unknow.cell] <- 'Unassigned'
    true.labels <- unique(true.tag)
    our.tag <- scRef.tag
    weighted_macro_f1 <- metrics$f1_score(true.tag, our.tag, average = 'weighted', labels = true.labels)
    macro_f1 <- metrics$f1_score(true.tag, our.tag, average = 'macro', labels = true.labels)
    accuracy <- metrics$accuracy_score(true.tag, our.tag)
    
    f1 <- c()
    for (label in true.labels) {
        tmp.true.tag <- true.tag
        tmp.our.tag <- our.tag
        tmp.true.tag[tmp.true.tag != label] <- '0'
        tmp.our.tag[tmp.our.tag != label] <- '0'
        sub.f1 <- metrics$f1_score(tmp.true.tag, tmp.our.tag, average = 'binary', pos_label = label)
        f1 <- c(f1, sub.f1)
    }
    names(f1) <- true.labels
    
    out <- list()
    out$weighted_macro_f1 <- weighted_macro_f1
    out$macro_f1 <- macro_f1
    out$accuracy <- accuracy
    out$f1 <- f1
    out$conf <- table(true.tag, our.tag)
    
    return(out)
    
}

# reference
file.ref <- '/home/zy/scRef/try_data/scRef/Reference/MouseBrain_Bulk_Zhang2014/Reference_expression.txt'
file.ref <- file.ref
exp_ref_mat.origin <- read.table(file.ref, header=T, row.names=1, sep='\t', check.name=F)
exp_ref_mat <- exp_ref_mat.origin


############# regard sc-counts data as reference
path.input <- '/home/zy/scRef/summary/'
path.output <- '/home/zy/scRef/Benchmark/mouse_brain/'
dataset <- 'Tasic'
file.data.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat_cluster_merged.txt')
# OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, del.label = c('Unclassified'))
# saveRDS(OUT, file = paste0(path.output, dataset, '.Rdata'))
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_Tasic <- OUT$data.filter
label_Tasic <- OUT$label.filter
ref.labels <-label_Tasic$label.unlabeled.use.cols...
ref.mtx <- exp_Tasic
ref.dataset <- 'Tasic'

############### import unlabeled data
############### HochgernerA
file.in <- '/home/zy/scRef/sc_data/Hochgerner/GSE95315_10X_expression_data_v2.tab'
path.output <- '/home/zy/scRef/Benchmark/mouse_brain/'
dataset <- 'HochgernerA'
# df.in <- read.delim(file.in, row.names = 1)
# label.filter <- t(df.in[1:2, ])
# data.filter <- df.in[-(1:2),]
# genes <- rownames(data.filter)
# data.filter <- apply(data.filter, 2, as.numeric)
# rownames(data.filter) <- genes
# OUT <- list()
# OUT$label.filter <- label.filter
# OUT$data.filter <- data.filter
# saveRDS(OUT, file = paste0(path.output, dataset, '.Rdata'))
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_Hochgerner <- OUT$data.filter
label_Hochgerner <- OUT$label.filter
exp_sc_mat <- exp_Hochgerner
label_sc <- label_Hochgerner

ref.names <- unique(ref.labels)
# list of cell names
all.cell <- unique(label_sc[,2])
uniform.names <- c("Neuron", "Endothelial Cell", "Astrocyte", "Microglia", 
                   "Oligodendrocyte", "OPC")
df.ref.names <- data.frame(ref.name = ref.names, name = uniform.names)
uniform.names <- c("Endothelial Cell", "Unassigned", "Unassigned", "Microglia", "Unassigned", 
                   "Oligodendrocyte", "Unassigned", "OPC", "Astrocyte", "Unassigned", 
                   "Unassigned", "Unassigned", "Unassigned", "Neuron", "Neuron", 
                   "Neuron", "Neuron", "Neuron", "Neuron", "Neuron",
                   "Neuron", "Unassigned")
df.sc.names <- data.frame(sc.name = all.cell, name = uniform.names)


### scRef
source('/home/zy/my_git/scRef/main/scRef.v17.R')
setwd('~/my_git/scRef')
# result.scref <- SCREF(exp_sc_mat, exp_ref_mat,
#                       type_ref = 'sum-counts', use.RUVseq = T, 
#                       cluster.speed = T, cluster.cell = 10,
#                       min_cell = 10, CPU = 8)
result.scref <- SCREF(exp_sc_mat, ref.mtx, ref.labels,
                      type_ref = 'sc-counts', use.RUVseq = T, 
                      # method1 = 'kendall',
                      cluster.speed = T, cluster.cell = 10,
                      min_cell = 10, CPU = 8, GMM.ceiling_cutoff = 30)
pred.scRef <- result.scref$final.out$scRef.tag
table(label_Habib$label.unlabeled.use.cols..., pred.scRef)

### UMAP with label
library(ggplot2)
plot.umap <- supervised.UMAP(exp_sc_mat, pred.scRef)
path.plot <- '/home/zy/scRef/UMAP_label'
ggsave(plot.umap, path = path.plot, filename = 'UMAP_label.png')

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
seurat.unlabeled@meta.data$original.label <- label_sc[use.cells, 'cluster_name']
seurat.unlabeled@meta.data$scRef.tag <- pred.scRef
# seurat.unlabeled@meta.data$new.tag <- new.tag

# PCA
seurat.unlabeled <- RunPCA(seurat.unlabeled, npcs = 150, verbose = F)

# cluster
# seurat.unlabeled <- FindNeighbors(seurat.unlabeled, reduction = "pca", dims = 1:75, nn.eps = 0.5)
# seurat.unlabeled <- FindClusters(seurat.unlabeled, resolution = 3, n.start = 20)

# UMAP
seurat.unlabeled <- RunUMAP(seurat.unlabeled, dims = 1:150, n.neighbors = 100)
# figure1: ture label
plot.umap <- DimPlot(seurat.unlabeled, reduction = "umap", label = T, group.by = 'original.label')
# figure2: cluster label
DimPlot(seurat.unlabeled, reduction = "umap", label = T, group.by = 'seurat_clusters')
# figure3: scRef plus label
plot.umap.scRef <- DimPlot(seurat.unlabeled, reduction = "umap", label = T, group.by = 'scRef.tag')

# umap with label
# transform character label to index
labels <- pred.scRef
vec.labels <- c(setdiff(unique(labels), 'Unassigned'), 'Unassigned')
df.label.idx <- data.frame(label = vec.labels, idx = c(1:(length(vec.labels) - 1), -1))
for (j in 1:dim(df.label.idx)[1]) {
    labels[labels == df.label.idx[j, 'label']] <- df.label.idx[j, 'idx']
}

labels <- as.factor(label_sc[use.cells,])

# supervised UMAP
umap <- import('umap')
class.umap <- umap$UMAP(target_weight = 0)
mat.pca <- seurat.unlabeled@reductions$pca@cell.embeddings
embedding <- class.umap$fit_transform(X = mat.pca, y = as.numeric(labels))
dimnames(embedding)[[1]] <- dimnames(mat.pca)[[1]]
umap.label <- CreateDimReducObject(embeddings = embedding, key = 'UMAP_')
seurat.unlabeled@reductions$umap.label <- umap.label

plot.umap.label.scRef <- DimPlot(seurat.unlabeled, reduction = "umap.label", label = T, group.by = 'scRef.tag')

plot.umap.label <- DimPlot(seurat.unlabeled, reduction = "umap.label", label = T, group.by = 'original.label')

library(ggplot2)
pathout <- '/home/zy/scRef/figure'
ggsave(plot = plot.umap, path = pathout, filename = 'umap.png', 
       units = 'cm', height = 23, width = 30)
ggsave(plot = plot.umap.scRef, path = pathout, filename = 'umap_scRef.png', 
       units = 'cm', height = 24, width = 33)
ggsave(plot = plot.umap.label, path = pathout, filename = 'umap_label.png', 
       units = 'cm', height = 23, width = 30)
ggsave(plot = plot.umap.label.scRef, path = pathout, filename = 'umap_label_scRef.png', 
       units = 'cm', height = 24, width = 33)
