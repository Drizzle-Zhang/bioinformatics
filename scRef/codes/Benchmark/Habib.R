library(reticulate)
library(foreach)
library(doParallel)
# library(metap)

# import python package: sklearn.metrics
use_python('/home/zy/tools/anaconda3/bin/python3', required = T)
# py_config()
py_module_available('sklearn')
metrics <- import('sklearn.metrics')

setwd('/home/zy/scRef/try_data')
# input file
file.ref <- './scRef/Reference/MouseBrain_Bulk_Zhang2014/Reference_expression.txt'
# parameters
num.cpu <- 8

# reference file
file.ref <- file.ref
exp_ref_mat.origin <- read.table(file.ref, header=T, row.names=1, sep='\t', check.name=F)
# MCA cell name
ref.names <- names(exp_ref_mat.origin)
# folder saving results
path.out <- '/home/zy/scRef/try_data/Habib/'
if (!file.exists(path.out)) {
    dir.create(path.out)
}

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

# rate of error correction
REC <- function(true.tag, scRef.tag, new.tag, method = 'f1') {
    error.tag <- c()
    for (idx in 1:length(true.tag)) {
        if (true.tag[idx] == scRef.tag[idx]) {
            error.tag <- c(error.tag, 'known')
        } else {
            error.tag <- c(error.tag, 'unknown')
        }
    }
    new.error.tag <- new.tag
    new.error.tag[new.error.tag != 'unknown'] <- 'known'
    if (method == 'f1') {
        score <- metrics$f1_score(error.tag, new.error.tag, 
                                  average = 'binary', pos_label = 'unknown')
    }
    if (method == 'precision') {
        score <- metrics$precision_score(error.tag, new.error.tag, pos_label = 'unknown')
    }
    if (method == 'recall') {
        score <- metrics$recall_score(error.tag, new.error.tag, pos_label = 'unknown')
    }
    if (method == 'accuracy') {
        score <- metrics$accuracy_score(error.tag, new.error.tag)
    }
    
    return(score)
    
}

######################### Habib
dataset <- 'Habib'
file.data.unlabeled <- paste0('./summary/', dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0('./summary/', dataset, '_exp_sc_mat_cluster_original.txt')
# OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, del.label = c('miss'))
# 
# saveRDS(OUT, file = './Habib/Habib.Rdata')
OUT <- readRDS('./Habib/Habib.Rdata')
data.filter <- OUT$data.filter
label.filter <- OUT$label.filter

# list of cell names
all.cell <- unique(label.filter[,1])
sc.name <- c("Astrocyte", "Neurons", "OPC",
             "Oligodend", "Oligodend",
             "microglia", "EndothelialCells")
unknow.cell <- setdiff(all.cell, sc.name)
df.cell.names <- data.frame(ref.name = ref.names, sc.name = sc.name, idx = 1:length(sc.name))

exp_ref_mat <- exp_ref_mat.origin

# match MCA names with ref names
ref_MCA_names <- data.frame(MCA.name = c("Astrocyte", "Neuron", "Oligodendrocyte precursor cell",
                                          "Oligodendrocyte", "Myelinating oligodendrocyte",
                                          "Microglia", "Endothelial cell"),
                            ref.name = ref.names)

# library(Seurat)
# # data preparing
# seurat.unlabeled <- CreateSeuratObject(counts = data.filter)
# seurat.unlabeled <- NormalizeData(seurat.unlabeled, normalization.method = "LogNormalize",
#                                   scale.factor = 10000)
# seurat.unlabeled <- FindVariableFeatures(seurat.unlabeled, selection.method = "vst", nfeatures = 8000)
# use.genes <- VariableFeatures(seurat.unlabeled)
# # VariableFeaturePlot(seurat.unlabeled)
# exp_sc_mat <- data.filter[use.genes,]
exp_sc_mat <- data.filter
# df.Habib <- data.filter

# run scRef
source('/home/zy/my_git/scRef/main/scRef.v8.R')
setwd('~/my_git/scRef')
result.scref <- SCREF(exp_sc_mat, exp_ref_mat, type_ref = 'fpkm', 
                      cluster.speed = T, cluster.cell = 10,
                      min_cell = 30, CPU = num.cpu)
meta.tag <- merge(result.scref$final.out, label.filter, by = 'row.names')
row.names(meta.tag) <- meta.tag$Row.names
meta.tag$Row.names <- NULL
names(meta.tag) <- c('scRef.tag', 'log10Pval', 'ori.tag')

# save txt file
meta.tag <- meta.tag[order(meta.tag$log10Pval),]
write.table(meta.tag, 
            paste0(path.out, 'tags_', dataset, '_scRef', '.txt'),
            sep = '\t', quote = F)

# evaluation
ori.tag = meta.tag$ori.tag
scRef.tag = meta.tag$scRef.tag
# uniform tags
for (j in 1:dim(df.cell.names)[1]) {
    scRef.tag[scRef.tag == df.cell.names[j, 'ref.name']] <- 
        df.cell.names[j, 'sc.name']
}
meta.tag$scRef.tag <- scRef.tag

# default cutoff
true.tag <- meta.tag$ori.tag
true.tag[true.tag %in% unknow.cell] <- 'Unassigned'
our.tag <- meta.tag$scRef.tag
metrics$f1_score(true.tag, our.tag, average = 'weighted')
metrics$f1_score(true.tag, our.tag, average = 'macro')
metrics$accuracy_score(true.tag, our.tag)

use.cells <- dimnames(data.filter)[[2]]
mtx.tag <- as.matrix(result.scref$final.out)
scRef.tags <- mtx.tag[use.cells, 'scRef.tag'] 
plot.Habib <- supervised.UMAP(data.filter, scRef.tags)

vec.cutoff <- 1:90

one.eval <- function(cutoff, meta.tag) {
    df.sub <- data.frame(stringsAsFactors = F)
    # library(reticulate)
    # use_python('/home/zy/tools/anaconda3/bin/python3', required = T)
    # metrics <- import('sklearn.metrics')
    true.tag <- meta.tag$ori.tag
    true.tag[true.tag %in% unknow.cell] <- 'Unassigned'
    our.tag <- meta.tag$scRef.tag
    our.tag[meta.tag$log10Pval <= cutoff] <- 'Unassigned'
    df.sub[1, 'cutoff'] <- cutoff
    df.sub[1, 'weighted.f1'] <-
        metrics$f1_score(true.tag, our.tag, average = 'weighted')
    df.sub[1, 'macro.f1'] <-
        metrics$f1_score(true.tag, our.tag, average = 'macro')
    df.sub[1, 'micro.f1'] <-
        metrics$f1_score(true.tag, our.tag, average = 'micro')
    df.sub[1, 'accuracy'] <-
        metrics$accuracy_score(true.tag, our.tag)
    df.sub[1, 'REC.f1'] <-
        REC(true.tag, meta.tag$scRef.tag, our.tag, method = 'f1')
    df.sub[1, 'REC.precision'] <-
        REC(true.tag, meta.tag$scRef.tag, our.tag, method = 'precision')
    df.sub[1, 'REC.recall'] <-
        REC(true.tag, meta.tag$scRef.tag, our.tag, method = 'recall')
    df.sub[1, 'REC.accuracy'] <-
        REC(true.tag, meta.tag$scRef.tag, our.tag, method = 'accuracy')
    gc()
    return(df.sub)

}
registerDoParallel(5)
df.metrics <- foreach(cutoff = vec.cutoff, .combine = rbind) %dopar% one.eval(cutoff, meta.tag)
stopImplicitCluster()
# cl= makeCluster(num.cpu, outfile='')
# list.sub <- parLapply(cl = cl, vec.cutoff, one.eval, meta.tag = meta.tag,
#                       unknow.cell = unknow.cell)
# stopCluster(cl)
# df.metrics = data.frame()
# for(df.sub in RUN){
#     df.metrics=rbind(df.metrics, df.sub)}


# best cutoff
best.cutoff <- df.metrics$cutoff[df.metrics$weighted.f1 == max(df.metrics$weighted.f1)][1]
new.tag <- meta.tag$scRef.tag
new.tag[meta.tag$log10Pval <= best.cutoff] <- 'unknown'
meta.tag$new.tag <- new.tag
print(best.cutoff)

# cluster
library(Seurat)
# data preparing
seurat.unlabeled <- CreateSeuratObject(counts = data.filter)
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
seurat.unlabeled@meta.data$original.label <- label.filter[use.cells,]
mtx.tag <- as.matrix(meta.tag)
seurat.unlabeled@meta.data$scRef.tag <- mtx.tag[use.cells, 'scRef.tag'] 
# seurat.unlabeled@meta.data$new.tag <- new.tag

# PCA
seurat.unlabeled <- RunPCA(seurat.unlabeled, npcs = 75)

# cluster
seurat.unlabeled <- FindNeighbors(seurat.unlabeled, reduction = "pca", dims = 1:75, nn.eps = 0.5)
seurat.unlabeled <- FindClusters(seurat.unlabeled, resolution = 3, n.start = 20)

# UMAP
seurat.unlabeled <- RunUMAP(seurat.unlabeled, dims = 1:20, n.neighbors = 30)
# figure1: ture label
DimPlot(seurat.unlabeled, reduction = "umap", label = T, group.by = 'original.label')
# figure2: cluster label
DimPlot(seurat.unlabeled, reduction = "umap", label = T, group.by = 'seurat_clusters')
# figure3: scRef plus label
DimPlot(seurat.unlabeled, reduction = "umap", label = T, group.by = 'scRef.tag')

# Supervised UMAP
# use python in R
library(reticulate)
use_python('/home/zy/tools/anaconda3/bin/python3', required = T)
# py_config()
# import python package: umap
py_module_available('umap')
umap <- import('umap')
class.umap <- umap$UMAP()

mat.pca <- seurat.unlabeled@reductions$pca@cell.embeddings
label.seurat <- seurat.unlabeled@meta.data$scRef.tag
for (j in 1:dim(df.cell.names)[1]) {
    label.seurat[label.seurat == df.cell.names[j, 'sc.name']] <- df.cell.names[j, 'idx']
}
label.seurat[label.seurat == 'Unassigned'] <- '-1'
label.seurat <- as.numeric(label.seurat)
embedding <- class.umap$fit_transform(X = mat.pca, y = label.seurat)
dimnames(embedding)[[1]] <- dimnames(mat.pca)[[1]]

umap.label <- CreateDimReducObject(embeddings = embedding, key = 'UMAP_')
seurat.unlabeled@reductions$umap.label <- umap.label
DimPlot(seurat.unlabeled, reduction = "umap.label", label = T, group.by = 'scRef.tag')
DimPlot(seurat.unlabeled, reduction = "umap.label", label = T, group.by = 'original.label')

mtx.tag <- as.matrix(meta.tag)
label.in <- mtx.tag[use.cells, 'scRef.tag']
plot.label <- supervised.UMAP(data.filter, label.in)
