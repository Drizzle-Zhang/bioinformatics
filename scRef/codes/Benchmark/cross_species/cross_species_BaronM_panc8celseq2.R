# import python package: sklearn.metrics
library(reticulate)
use_python('/local/zy/tools/anaconda3/bin/python3', required = T)
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
    OUT$mat_exp <- data.filter
    OUT$label <- label.filter
    return(OUT)

}

# evaluation
simple.evaluation <- function(true.tag, scRef.tag, df.ref.names, df.sc.names) {
    # uniform tags
    for (j in 1:dim(df.ref.names)[1]) {
        scRef.tag[scRef.tag == df.ref.names[j, 'ref.name']] <- df.ref.names[j, 'name']
    }
    scRef.tag[!(scRef.tag %in% df.ref.names$name)] <- 'Unassigned'
    for (j in 1:dim(df.sc.names)[1]) {
        true.tag[true.tag == df.sc.names[j, 'sc.name']] <- df.sc.names[j, 'name']
    }

    percent.unassigned <- sum(scRef.tag == 'Unassigned')/sum(true.tag == 'Unassigned')


    # true.labels <- setdiff(unique(true.tag), 'Unassigned')
    true.labels <- unique(true.tag)
    our.tag <- scRef.tag
    weighted_macro_f1 <- metrics$f1_score(true.tag, our.tag, average = 'weighted', labels = true.labels)
    macro_f1 <- metrics$f1_score(true.tag, our.tag, average = 'macro', labels = true.labels)
    accuracy <- metrics$accuracy_score(true.tag, our.tag)
    balanced_acc <- metrics$balanced_accuracy_score(true.tag, our.tag)
    # rm unassigned in tags
    true.tag.rm <- true.tag[our.tag != 'Unassigned']
    our.tag.rm <- our.tag[our.tag != 'Unassigned']
    accuracy.rm.unassigned <- metrics$accuracy_score(our.tag.rm, true.tag.rm)
    macro_f1.rm.unassigned <- metrics$f1_score(true.tag.rm, our.tag.rm, average = 'macro', labels = unique(true.tag.rm))
    balanced.accuracy.rm.unassigned <-
        metrics$balanced_accuracy_score(our.tag.rm, true.tag.rm)

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

    our.labels <- setdiff(unique(our.tag), 'Unassigned')
    precision <- c()
    for (label in our.labels) {
        tmp.true.tag <- true.tag.rm
        tmp.our.tag <- our.tag.rm
        tmp.true.tag[tmp.true.tag != label] <- '0'
        tmp.our.tag[tmp.our.tag != label] <- '0'
        sub.precision <- metrics$precision_score(tmp.true.tag, tmp.our.tag, average = 'binary', pos_label = label)
        precision <- c(precision, sub.precision)

    }
    names(precision) <- our.labels
    mean.precision <- mean(precision)
    
    # remove unassigend in true tags
    true.tag.rm <- true.tag[true.tag != 'Unassigned']
    our.tag.rm <- our.tag[true.tag != 'Unassigned']
    accuracy.rm.newcell <- metrics$accuracy_score(true.tag.rm, our.tag.rm)
    balanced.accuracy.rm.newcell <- metrics$balanced_accuracy_score(true.tag.rm, our.tag.rm)
    
    out <- list()
    out$percent.unassigned <- percent.unassigned
    out$weighted_macro_f1 <- weighted_macro_f1
    out$macro_f1 <- macro_f1
    out$accuracy <- accuracy
    out$balanced.accuracy <- balanced_acc
    out$f1 <- f1
    out$accuracy.rm.unassigned <- accuracy.rm.unassigned
    out$macro_f1.rm.unassigned <- macro_f1.rm.unassigned
    out$precision.rm.unassigned <- precision
    out$balanced.accuracy.rm.unassigned <- balanced.accuracy.rm.unassigned
    out$mean.precision.rm.unassigned <- mean.precision
    out$accuracy.rm.newcell <- accuracy.rm.newcell
    out$balanced.accuracy.rm.newcell <- balanced.accuracy.rm.newcell
    out$conf <- table(true.tag, our.tag)

    return(out)

}

source('/local/zy/my_git/scRef/main/scRef.v20.R')

############# regard sc-counts data as reference
path.input <- '/mdshare/node9/zy/scRef/sc_data/'
path.output <- '/mdshare/node9/zy/scRef/Benchmark/cross_species/'
dataset <- 'BaronM'
file.data.unlabeled <- paste0(path.input, dataset, '/cell_exp.txt')
file.label.unlabeled <- paste0(path.input, dataset, '/cell_meta.txt')
# OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, del.label = c('Unclassified'))
# saveRDS(OUT, file = paste0(path.output, dataset, '.Rdata'))
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
ref.mtx <- OUT$mat_exp
ref.labels <- OUT$label[,1]
ref.dataset <- 'BaronM'

############### import unlabeled data
############### BaronH
library(Seurat)
library(SeuratData)
data("panc8")
dataset <- 'panc8_celseq2'
file.save <- paste0(path.output, dataset, '.Rdata')
# OUT <- list()
# OUT$mat_exp <- as.matrix(panc8@assays$RNA@counts[, panc8$dataset %in% c('celseq2')])
# OUT$label <- data.frame(
#     annotations = as.character(panc8$celltype)[panc8$dataset %in% c('celseq2')],
#     row.names = colnames(OUT$mat_exp))
# saveRDS(OUT, file = file.save)
OUT <- readRDS(file.save)
exp_sc_mat <- OUT$mat_exp
label_sc <- OUT$label

ref.names <- names(table(ref.labels))
# list of cell names
all.cell <- names(table(label_sc[,1]))
df.ref.names <- data.frame(ref.name = ref.names, name = ref.names)
uniform.names <- c("Unassigned", "activated_stellate", "alpha", "beta",
                   "delta", "ductal", "endothelial", "Unassigned",
                   "gamma", "macrophage", "Unassigned", "quiescent_stellate",
                   "schwann")
df.sc.names <- data.frame(sc.name = all.cell, name = uniform.names)

# run methods
#############################################
### scRef
source('/local/zy/my_git/scRef/main/scRef.v20.R')
exp_sc_mat <- transform.HomoloGene(exp_sc_mat)
setwd('/local/zy/my_git/scRef')
result.scref <- SCREF(exp_sc_mat, ref.mtx, ref.labels,
                      type_ref = 'sc-counts', use.RUVseq = T,
                      cluster.speed = F, cluster.resolution = 1,
                      GMM.num_cluster = 3, GMM.floor_cutoff = 2, GMM.ceiling_cutoff = 10,
                      min_cell = 1, CPU = 10)
pred.scMAGIC <- result.scref$final.out$scRef.tag
saveRDS(pred.scMAGIC, file = paste0(path.output, ref.dataset, '_', dataset, '_scMAGIC.Rdata'))

### sciBet
suppressMessages(library(tidyverse))
suppressMessages(library(scibet))
suppressMessages(library(viridis))
suppressMessages(library(ggsci))
train_set <- as.data.frame(t(ref.mtx))
train_set$label <- ref.labels
test_set <- as.data.frame(t(exp_sc_mat))
sciBet <- SciBet(train_set, test_set)
saveRDS(sciBet, file = paste0(path.output, ref.dataset, '_', dataset, '_sciBet.Rdata'))

### singleCellNet
library(singleCellNet)
library(dplyr)
out <- .get_overlap_genes(exp_sc_mat, ref.mtx)
train_set <- as.matrix(out$exp_ref_mat)
LabelsTrain <- data.frame(Annotation = ref.labels, row.names = colnames(ref.mtx))
test_set <- as.matrix(out$exp_sc_mat)
class_info <- scn_train(stTrain = LabelsTrain, expTrain = train_set, dLevel = "Annotation")
classRes <- scn_predict(cnProc=class_info[['cnProc']], expDat=test_set, nrand = 50)
classRes <- classRes[, colnames(test_set)]
tags.singleCellNet <- rownames(classRes)[apply(classRes,2,which.max)]
tags.singleCellNet[tags.singleCellNet == 'rand'] <- 'Unassigned'
saveRDS(tags.singleCellNet, file = paste0(path.output, ref.dataset, '_', dataset, '_singleCellNet.Rdata'))

### singleR
library(SingleR)
train_set <- as.matrix(ref.mtx)
test_set <- as.matrix(exp_sc_mat)
singler = SingleR(method = "single", sc_data = test_set,
                  ref_data = train_set,
                  types = ref.labels, numCores = 4)
pred.singleR <- singler$labels
saveRDS(pred.singleR, file = paste0(path.output, ref.dataset, '_', dataset, '_singleR.Rdata'))

###  scmap
library(scmap)
library(SingleCellExperiment)
out <- .get_overlap_genes(exp_sc_mat, ref.mtx)
train_set <- as.matrix(out$exp_ref_mat)
test_set <- as.matrix(out$exp_sc_mat)
sce <- SingleCellExperiment(list(normcounts = train_set),
                            colData = data.frame(cell_type1 = ref.labels))
logcounts(sce) <- log2(normcounts(sce) + 1)
# use gene names as feature symbols
rowData(sce)$feature_symbol <- rownames(sce)
sce <- selectFeatures(sce, suppress_plot = TRUE)

sce_test <- SingleCellExperiment(list(normcounts = test_set))
logcounts(sce_test) <- log2(normcounts(sce_test) + 1)
rowData(sce_test)$feature_symbol <- rownames(sce_test)
sce_test@rowRanges@elementMetadata@listData = sce@rowRanges@elementMetadata@listData

# scmap-cluster
sce <- indexCluster(sce)
scmapCluster_results <- scmapCluster(projection = sce_test,index_list = list(metadata(sce)$scmap_cluster_index))
pred.scmap.cluster <- scmapCluster_results$combined_labs
saveRDS(pred.scmap.cluster,
        file = paste0(path.output, ref.dataset, '_', dataset, '_scmap-cluster.Rdata'))

# scmap-cell
set.seed(1)
sce <- indexCell(sce)
scmapCell_results <- scmapCell(sce_test,list(metadata(sce)$scmap_cell_index))
scmapCell_clusters <- scmapCell2Cluster(scmapCell_results,list(as.character(colData(sce)$cell_type1)))
pred.scmap.cell <- scmapCell_clusters$combined_labs
saveRDS(pred.scmap.cell,
        file = paste0(path.output, ref.dataset, '_', dataset, '_scmap-cell.Rdata'))

### CHETAH
library(CHETAH)
library(SingleCellExperiment)
sce <- SingleCellExperiment(assays = list(counts = ref.mtx),
                            colData = data.frame(celltypes = ref.labels))
sce_test <- SingleCellExperiment(list(counts = exp_sc_mat))
# sce_test <- CHETAHclassifier(input = sce_test, ref_cells = sce)
sce_test <- CHETAHclassifier(input = sce_test, ref_cells = sce,
                             n_genes = median(colSums(ref.mtx != 0))/2)
pred.CHETAH <- sce_test$celltype_CHETAH
saveRDS(pred.CHETAH,
        file = paste0(path.output, ref.dataset, '_', dataset, '_CHETAH.Rdata'))

### scPred
library("scPred")
library("Seurat")
library("magrittr")
# scPred Training
reference <- CreateSeuratObject(counts = ref.mtx)
reference@meta.data$cell_type <- ref.labels
reference <- reference %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30)
reference <- getFeatureSpace(reference, "cell_type")
reference <- trainModel(reference)
# scPred Prediction
query <- CreateSeuratObject(counts = exp_sc_mat)
query <- NormalizeData(query)
query <- scPredict(query, reference)
pred.scPred <- query@meta.data$scpred_prediction
saveRDS(pred.scPred,
        file = paste0(path.output, ref.dataset, '_', dataset, '_scPred.Rdata'))

### scID
library(scID)
library(Seurat)
Train_Labels <- list(ref.labels)
names(Train_Labels[[1]]) <- colnames(ref.mtx)
scID_output <- scid_multiclass(as.data.frame(exp_sc_mat), as.data.frame(ref.mtx), Train_Labels[[1]])
pred.scID <- scID_output$labels
saveRDS(pred.scID,
        file = paste0(path.output, ref.dataset, '_', dataset, '_scID.Rdata'))

### scClassify
library("scClassify")
library(Matrix)
exprsMat_train <- as(as.matrix(log1p(ref.mtx)), "dgCMatrix")
exp_sc_mat <- as(as.matrix(log1p(exp_sc_mat)), "dgCMatrix")
scClassify_res <- scClassify(exprsMat_train = exprsMat_train,
                             cellTypes_train = ref.labels,
                             exprsMat_test = list(one = exp_sc_mat),
                             tree = "HOPACH",
                             algorithm = "WKNN",
                             selectFeatures = c("limma"),
                             similarity = c("pearson"),
                             returnList = FALSE,
                             verbose = FALSE)
pred.scClassify <- scClassify_res$testRes$one$pearson_WKNN_limma$predRes
saveRDS(pred.scClassify,
        file = paste0(path.output, ref.dataset, '_', dataset, '_scClassify.Rdata'))

#############################################


# evaluation
#############################################
true.tags <- label_sc[,1]
df.plot <- data.frame(stringsAsFactors = F)

rda.scMAGIC <- paste0(path.output, ref.dataset, '_', dataset, '_scMAGIC.Rdata')
pred.scMAGIC <- readRDS(rda.scMAGIC)
res.scMAGIC <- simple.evaluation(true.tags, pred.scMAGIC, df.ref.names, df.sc.names)
df.sub <- data.frame(term = 'Weighted macro F1', method = 'scMAGIC',
                     value = res.scMAGIC$weighted_macro_f1, stringsAsFactors = F)
df.sub <- rbind(df.sub,
                data.frame(term = 'Macro F1', method = 'scMAGIC',
                           value = res.scMAGIC$macro_f1, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Accuracy', method = 'scMAGIC',
                           value = res.scMAGIC$accuracy, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Balanced accuracy', method = 'scMAGIC',
                           value = res.scMAGIC$balanced.accuracy, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'labeled-Accuracy', method = 'scMAGIC',
                           value = res.scMAGIC$accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'labeled-Balanced accuracy', method = 'scMAGIC',
                           value = res.scMAGIC$balanced.accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Percent of unassigned', method = 'scMAGIC',
                           value = res.scMAGIC$percent.unassigned, stringsAsFactors = F))
df.plot <- rbind(df.plot, df.sub)

rda.sciBet <- paste0(path.output, ref.dataset, '_', dataset, '_sciBet.Rdata')
pred.sciBet <- readRDS(rda.sciBet)
res.sciBet <- simple.evaluation(true.tags, pred.sciBet, df.ref.names, df.sc.names)
df.sub <- data.frame(term = 'Weighted macro F1', method = 'sciBet',
                     value = res.sciBet$weighted_macro_f1, stringsAsFactors = F)
df.sub <- rbind(df.sub,
                data.frame(term = 'Macro F1', method = 'sciBet',
                           value = res.sciBet$macro_f1, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Accuracy', method = 'sciBet',
                           value = res.sciBet$accuracy, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Balanced accuracy', method = 'sciBet',
                           value = res.sciBet$balanced.accuracy, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'labeled-Accuracy', method = 'sciBet',
                           value = res.sciBet$accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'labeled-Balanced accuracy', method = 'sciBet',
                           value = res.sciBet$balanced.accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Percent of unassigned', method = 'sciBet',
                           value = res.sciBet$percent.unassigned, stringsAsFactors = F))
df.plot <- rbind(df.plot, df.sub)

rda.singleCellNet <- paste0(path.output, ref.dataset, '_', dataset, '_singleCellNet.Rdata')
pred.singleCellNet <- readRDS(rda.singleCellNet)
res.singleCellNet <- simple.evaluation(true.tags, pred.singleCellNet, df.ref.names, df.sc.names)
df.sub <- data.frame(term = 'Weighted macro F1', method = 'singleCellNet',
                     value = res.singleCellNet$weighted_macro_f1, stringsAsFactors = F)
df.sub <- rbind(df.sub,
                data.frame(term = 'Macro F1', method = 'singleCellNet',
                           value = res.singleCellNet$macro_f1, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Accuracy', method = 'singleCellNet',
                           value = res.singleCellNet$accuracy, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Balanced accuracy', method = 'singleCellNet',
                           value = res.singleCellNet$balanced.accuracy, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'labeled-Accuracy', method = 'singleCellNet',
                           value = res.singleCellNet$accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'labeled-Balanced accuracy', method = 'singleCellNet',
                           value = res.singleCellNet$balanced.accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Percent of unassigned', method = 'singleCellNet',
                           value = res.singleCellNet$percent.unassigned, stringsAsFactors = F))
df.plot <- rbind(df.plot, df.sub)

rda.singleR <- paste0(path.output, ref.dataset, '_', dataset, '_singleR.Rdata')
pred.singleR <- readRDS(rda.singleR)
res.singleR <- simple.evaluation(true.tags, pred.singleR, df.ref.names, df.sc.names)
df.sub <- data.frame(term = 'Weighted macro F1', method = 'singleR',
                     value = res.singleR$weighted_macro_f1, stringsAsFactors = F)
df.sub <- rbind(df.sub,
                data.frame(term = 'Macro F1', method = 'singleR',
                           value = res.singleR$macro_f1, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Accuracy', method = 'singleR',
                           value = res.singleR$accuracy, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Balanced accuracy', method = 'singleR',
                           value = res.singleR$balanced.accuracy, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'labeled-Accuracy', method = 'singleR',
                           value = res.singleR$accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'labeled-Balanced accuracy', method = 'singleR',
                           value = res.singleR$balanced.accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Percent of unassigned', method = 'singleR',
                           value = res.singleR$percent.unassigned, stringsAsFactors = F))
df.plot <- rbind(df.plot, df.sub)

rda.scmap.cluster <- paste0(path.output, ref.dataset, '_', dataset, '_scmap-cluster.Rdata')
pred.scmap.cluster <- readRDS(rda.scmap.cluster)
res.scmap.cluster <- simple.evaluation(true.tags, pred.scmap.cluster, df.ref.names, df.sc.names)
df.sub <- data.frame(term = 'Weighted macro F1', method = 'scmap-cluster',
                     value = res.scmap.cluster$weighted_macro_f1, stringsAsFactors = F)
df.sub <- rbind(df.sub,
                data.frame(term = 'Macro F1', method = 'scmap-cluster',
                           value = res.scmap.cluster$macro_f1, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Accuracy', method = 'scmap-cluster',
                           value = res.scmap.cluster$accuracy, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Balanced accuracy', method = 'scmap-cluster',
                           value = res.scmap.cluster$balanced.accuracy, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'labeled-Accuracy', method = 'scmap-cluster',
                           value = res.scmap.cluster$accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'labeled-Balanced accuracy', method = 'scmap-cluster',
                           value = res.scmap.cluster$balanced.accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Percent of unassigned', method = 'scmap-cluster',
                           value = res.scmap.cluster$percent.unassigned, stringsAsFactors = F))
df.plot <- rbind(df.plot, df.sub)

rda.scmap.cell <- paste0(path.output, ref.dataset, '_', dataset, '_scmap-cell.Rdata')
pred.scmap.cell <- readRDS(rda.scmap.cell)
res.scmap.cell <- simple.evaluation(true.tags, pred.scmap.cell, df.ref.names, df.sc.names)
df.sub <- data.frame(term = 'Weighted macro F1', method = 'scmap-cell',
                     value = res.scmap.cell$weighted_macro_f1, stringsAsFactors = F)
df.sub <- rbind(df.sub,
                data.frame(term = 'Macro F1', method = 'scmap-cell',
                           value = res.scmap.cell$macro_f1, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Accuracy', method = 'scmap-cell',
                           value = res.scmap.cell$accuracy, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Balanced accuracy', method = 'scmap-cell',
                           value = res.scmap.cell$balanced.accuracy, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'labeled-Accuracy', method = 'scmap-cell',
                           value = res.scmap.cell$accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'labeled-Balanced accuracy', method = 'scmap-cell',
                           value = res.scmap.cell$balanced.accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Percent of unassigned', method = 'scmap-cell',
                           value = res.scmap.cell$percent.unassigned, stringsAsFactors = F))
df.plot <- rbind(df.plot, df.sub)

rda.CHETAH <- paste0(path.output, ref.dataset, '_', dataset, '_CHETAH.Rdata')
pred.CHETAH <- readRDS(rda.CHETAH)
res.CHETAH <- simple.evaluation(true.tags, pred.CHETAH, df.ref.names, df.sc.names)
df.sub <- data.frame(term = 'Weighted macro F1', method = 'CHETAH',
                     value = res.CHETAH$weighted_macro_f1, stringsAsFactors = F)
df.sub <- rbind(df.sub,
                data.frame(term = 'Macro F1', method = 'CHETAH',
                           value = res.CHETAH$macro_f1, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Accuracy', method = 'CHETAH',
                           value = res.CHETAH$accuracy, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Balanced accuracy', method = 'CHETAH',
                           value = res.CHETAH$balanced.accuracy, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'labeled-Accuracy', method = 'CHETAH',
                           value = res.CHETAH$accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'labeled-Balanced accuracy', method = 'CHETAH',
                           value = res.CHETAH$balanced.accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Percent of unassigned', method = 'CHETAH',
                           value = res.CHETAH$percent.unassigned, stringsAsFactors = F))
df.plot <- rbind(df.plot, df.sub)

rda.scPred <- paste0(path.output, ref.dataset, '_', dataset, '_scPred.Rdata')
pred.scPred <- readRDS(rda.scPred)
res.scPred <- simple.evaluation(true.tags, pred.scPred, df.ref.names, df.sc.names)
df.sub <- data.frame(term = 'Weighted macro F1', method = 'scPred',
                     value = res.scPred$weighted_macro_f1, stringsAsFactors = F)
df.sub <- rbind(df.sub,
                data.frame(term = 'Macro F1', method = 'scPred',
                           value = res.scPred$macro_f1, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Accuracy', method = 'scPred',
                           value = res.scPred$accuracy, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Balanced accuracy', method = 'scPred',
                           value = res.scPred$balanced.accuracy, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'labeled-Accuracy', method = 'scPred',
                           value = res.scPred$accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'labeled-Balanced accuracy', method = 'scPred',
                           value = res.scPred$balanced.accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Percent of unassigned', method = 'scPred',
                           value = res.scPred$percent.unassigned, stringsAsFactors = F))
df.plot <- rbind(df.plot, df.sub)

rda.scID <- paste0(path.output, ref.dataset, '_', dataset, '_scID.Rdata')
pred.scID <- readRDS(rda.scID)
res.scID <- simple.evaluation(true.tags, pred.scID, df.ref.names, df.sc.names)
df.sub <- data.frame(term = 'Weighted macro F1', method = 'scID',
                     value = res.scID$weighted_macro_f1, stringsAsFactors = F)
df.sub <- rbind(df.sub,
                data.frame(term = 'Macro F1', method = 'scID',
                           value = res.scID$macro_f1, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Accuracy', method = 'scID',
                           value = res.scID$accuracy, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Balanced accuracy', method = 'scID',
                           value = res.scID$balanced.accuracy, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'labeled-Accuracy', method = 'scID',
                           value = res.scID$accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'labeled-Balanced accuracy', method = 'scID',
                           value = res.scID$balanced.accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Percent of unassigned', method = 'scID',
                           value = res.scID$percent.unassigned, stringsAsFactors = F))
df.plot <- rbind(df.plot, df.sub)

rda.scClassify <- paste0(path.output, ref.dataset, '_', dataset, '_scClassify.Rdata')
pred.scClassify <- readRDS(rda.scClassify)
res.scClassify <- simple.evaluation(true.tags, pred.scClassify, df.ref.names, df.sc.names)
df.sub <- data.frame(term = 'Weighted macro F1', method = 'scClassify',
                     value = res.scClassify$weighted_macro_f1, stringsAsFactors = F)
df.sub <- rbind(df.sub,
                data.frame(term = 'Macro F1', method = 'scClassify',
                           value = res.scClassify$macro_f1, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Accuracy', method = 'scClassify',
                           value = res.scClassify$accuracy, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Balanced accuracy', method = 'scClassify',
                           value = res.scClassify$balanced.accuracy, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'labeled-Accuracy', method = 'scClassify',
                           value = res.scClassify$accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'labeled-Balanced accuracy', method = 'scClassify',
                           value = res.scClassify$balanced.accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub,
                data.frame(term = 'Percent of unassigned', method = 'scClassify',
                           value = res.scClassify$percent.unassigned, stringsAsFactors = F))
df.plot <- rbind(df.plot, df.sub)

# sort
df.acc <- df.plot[df.plot$term == 'Accuracy', ]
df.plot$method <- factor(df.plot$method,
                         levels = df.acc$method[order(df.acc$value, decreasing = T)])

# save results
file.res <- paste0(path.output, 'results_', ref.dataset, '_', dataset, '.txt')
write.table(df.plot, file = file.res, sep = '\t')

# #########################################
# rda.scMAGIC <- paste0(path.output, ref.dataset, '_', dataset, '_scMAGIC.Rdata')
# pred.scMAGIC <- readRDS(rda.scMAGIC)
#
# true.tags <- label_sc$annotations
# table(true.tags, pred.scMAGIC)
#
# library(ggplot2)
# path.res <- '/home/zy/scRef/figure/cross_species'
#
# # heatmap
# method <- 'scMAGIC'
# mytable <- table(true.tags, pred.scMAGIC)
# mydata <- data.frame(stringsAsFactors = F)
# table.true <- table(true.tags)
# for (label1 in rownames(mytable)) {
#     row.sum <- table.true[label1]
#     for (label2 in colnames(mytable)) {
#         mydata <- rbind(mydata, data.frame(origin = label1, annotation = label2,
#                                            count = mytable[label1, label2],
#                                            prop = mytable[label1, label2]/row.sum))
#     }
# }
# mydata$origin <- factor(mydata$origin,
#                         levels = c(colnames(mytable), "acinar", "epsilon", "mast"))
# mydata$annotation <- factor(mydata$annotation,
#                             levels = colnames(mytable))
#
# plot.heatmap <-
#     ggplot(data = mydata, aes(x = origin, y = annotation)) +
#     geom_tile(aes(fill = prop)) +
#     # scale_fill_continuous(low = "#FFFAFA", high = "#A52A2A") +
#     scale_fill_gradient2(low = "#4169E1", high = "#FF4500", mid = '#FFCC00', midpoint = 0.5) +
#     labs(fill = 'Proportion') +
#     theme_bw() +
#     theme(
#         axis.ticks = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_blank(),
#         axis.title = element_blank(),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
#     ) +
#     geom_text(aes(label = round(prop, 2)), family = "Arial", size = 2.5)
# ggsave(filename = paste0('heatmap_', ref.dataset, '_', dataset, '_', method, '.png'),
#        path = path.res, plot = plot.heatmap,
#        units = 'cm', height = 16, width = 22)
#
#
# ### original plot
# library(Seurat)
# # data preparing
# exp_sc_mat <- OUT$mat_exp
# seurat.unlabeled <- CreateSeuratObject(counts = exp_sc_mat)
# seurat.unlabeled <- NormalizeData(seurat.unlabeled, normalization.method = "LogNormalize",
#                                   scale.factor = 10000)
# seurat.unlabeled <- FindVariableFeatures(seurat.unlabeled, selection.method = "vst", nfeatures = 2000)
# seurat.unlabeled <- ScaleData(seurat.unlabeled)
#
# # add label
# use.cells <- dimnames(seurat.unlabeled@assays$RNA@counts)[[2]]
# names(label_sc) <- 'cell_type'
# seurat.unlabeled@meta.data$original.label <- label_sc[use.cells, 'cell_type']
# seurat.unlabeled@meta.data$scRef.tag <- pred.scMAGIC
#
# # PCA
# seurat.unlabeled <- RunPCA(seurat.unlabeled, npcs = 100, verbose = F)
#
# # UMAP
# seurat.unlabeled <- RunUMAP(seurat.unlabeled, dims = 1:20, n.neighbors = 30)
#
# # figure1: ture label
# plot.umap <-
#     DimPlot(seurat.unlabeled, reduction = "umap", label = T, repel = T, group.by = 'original.label') +
#     # xlim(-16, 14) +
#     theme_bw() +
#     theme(axis.text = element_text(size = 9),
#           panel.grid = element_blank(),
#           axis.title = element_text(size = 12),
#           legend.text = element_text(size = 11))
# ggsave(filename = paste0('cluster_', ref.dataset, '_', dataset, '.png'),
#        path = path.res, plot = plot.umap,
#        units = 'cm', height = 18, width = 24)
#
# # figure2: cluster label
# # DimPlot(seurat.unlabeled, reduction = "umap", label = T, group.by = 'seurat_clusters')
# # figure3: scRef plus label
# library(scales)
# plot.umap.scRef <-
#     DimPlot(seurat.unlabeled, reduction = "umap", label = T, repel = T, group.by = 'scRef.tag') +
#     scale_color_manual(values = c(hue_pal()(10), 'gray'),
#                        breaks = c(names(table(pred.scMAGIC)))) +
#     theme_bw() +
#     theme(axis.text = element_text(size = 9),
#           panel.grid = element_blank(),
#           axis.title = element_text(size = 12),
#           legend.text = element_text(size = 11))
# ggsave(filename = paste0('cluster_scRef_', ref.dataset, '_', dataset, '.png'),
#        path = path.res, plot = plot.umap.scRef,
#        units = 'cm', height = 18, width = 24)
