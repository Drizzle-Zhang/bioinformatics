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
    
    # percent.unassigned <- sum(our.tag == 'Unassigned')/length(our.tag)
    # f1.unassigned <- metrics$f1_score(true.tag, our.tag, 
    #                                   average = 'binary', pos_label = 'Unassigned')
    scRef.tag[!(scRef.tag %in% df.cell.names$sc.name)] <- 'Unassigned'
    
    true.labels <- unique(true.tag)
    our.tag <- scRef.tag
    weighted_macro_f1 <- metrics$f1_score(true.tag, our.tag, average = 'weighted', labels = true.labels)
    macro_f1 <- metrics$f1_score(true.tag, our.tag, average = 'macro', labels = true.labels)
    accuracy <- metrics$accuracy_score(true.tag, our.tag)
    balanced_acc <- metrics$balanced_accuracy_score(true.tag, our.tag)
    # rm unassigned in tags
    
    true.tag.rm <- true.tag[our.tag != 'Unassigned']
    our.tag.rm <- our.tag[our.tag != 'Unassigned']
    accuracy.rm.unassigned <- metrics$accuracy_score(true.tag.rm, our.tag.rm)
    balanced.accuracy.rm.unassigned <- 
        metrics$balanced_accuracy_score(true.tag.rm, our.tag.rm)
    
    f1 <- c()
    precision <- c()
    for (label in true.labels) {
        tmp.true.tag <- true.tag
        tmp.our.tag <- our.tag
        tmp.true.tag[tmp.true.tag != label] <- '0'
        tmp.our.tag[tmp.our.tag != label] <- '0'
        sub.f1 <- metrics$f1_score(tmp.true.tag, tmp.our.tag, average = 'binary', pos_label = label)
        f1 <- c(f1, sub.f1)
        sub.precision <- metrics$precision_score(tmp.true.tag, tmp.our.tag, average = 'binary', pos_label = label)
        precision <- c(precision, sub.precision)
    }
    names(f1) <- true.labels
    names(precision) <- true.labels
    names.rm.unassigned <- setdiff(true.labels, 'Unassigned')
    mean.precision <- mean(precision[names.rm.unassigned])
    
    out <- list()
    out$weighted_macro_f1 <- weighted_macro_f1
    out$macro_f1 <- macro_f1
    out$accuracy <- accuracy
    out$balanced.accuracy <- balanced_acc
    out$f1 <- f1
    out$precision.rm.unassigned <- precision
    out$accuracy.rm.unassigned <- accuracy.rm.unassigned
    out$balanced.accuracy.rm.unassigned <- balanced.accuracy.rm.unassigned
    out$mean.precision.rm.unassigned <- mean.precision
    out$conf <- table(true.tag, our.tag)
    
    return(out)
    
}


############# regard sc-counts data as reference
path.output <- '/mdshare/node9/zy/scRef/Benchmark/mouse_brain/'
ref.dataset <- 'Mizrak'
OUT <- readRDS(paste0(path.output, ref.dataset, '.Rdata'))
ref.mtx <- OUT$mat_exp
ref.labels <- OUT$label$Cluster

############### import unlabeled data
############### Habib
dataset <- 'Tasic'
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_sc_mat <- OUT$mat_exp
label_sc <- OUT$label$CellType

ref.names <- unique(ref.labels)
# list of cell names
all.cell <- unique(label_sc[,1])
sc.name <- c("Neurons", "EndothelialCells",
             "Astrocyte", "microglia", "Oligodend", "OPC")
unknow.cell <- setdiff(all.cell, sc.name)
df.cell.names <- data.frame(ref.name = ref.names, sc.name = sc.name, idx = 1:length(sc.name))

# run methods
#############################################
### scRef
library(scMAGIC)
output.scMAGIC <- scMAGIC(exp_sc_mat, ref.mtx, ref.labels,
                          type_ref = 'sc-counts', num_threads = 10)
pred.scMAGIC <- output.scMAGIC$scMAGIC.tag

# source('/local/zy/my_git/scRef/main/scRef.v20.R')
# setwd('~/my_git/scRef')
# result.scref <- SCREF(exp_sc_mat, ref.mtx, ref.labels,
#                       type_ref = 'sc-counts', use.RUVseq = T, 
#                       cluster.speed = T, cluster.cell = 5,
#                       min_cell = 10, CPU = 8)
# pred.scRef <- result.scref$final.out$scRef.tag
table(label_sc, pred.scMAGIC)
saveRDS(pred.scMAGIC, file = paste0(path.output, ref.dataset, '_', dataset, '_scRef.Rdata'))

### sciBet
suppressMessages(library(tidyverse))
suppressMessages(library(scibet))
suppressMessages(library(viridis))
suppressMessages(library(ggsci))
train_set <- apply(ref.mtx, 2, as.numeric)
rownames(train_set) <- rownames(ref.mtx)
train_set <- as.data.frame(t(train_set))
train_set$label <- ref.labels
test_set <- as.data.frame(t(exp_sc_mat))
sciBet <- SciBet(train_set, test_set)
table(label_sc, sciBet)
saveRDS(sciBet, file = paste0(path.output, ref.dataset, '_', dataset, '_sciBet.Rdata'))

### singleCellNet
library(singleCellNet)
library(dplyr)
out <- get_overlap_genes(exp_sc_mat, ref.mtx)
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
                  types = ref.labels, numCores = 1)
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
saveRDS(pred.scmap.cluster, 
        file = paste0(path.output, ref.dataset, '_', dataset, '_scmap-cell.Rdata'))

### CHETAH
library(CHETAH)
library(SingleCellExperiment)
sce <- SingleCellExperiment(assays = list(counts = ref.mtx), 
                            colData = data.frame(celltypes = ref.labels))
sce_test <- SingleCellExperiment(list(counts = exp_sc_mat))
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
scID_output <- scid_multiclass(exp_sc_mat, ref.mtx, Train_Labels[[1]])
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
true.tags <- label_sc$label.unlabeled.use.cols...
df.plot <- data.frame(stringsAsFactors = F)

rda.scRef <- paste0(path.output, ref.dataset, '_', dataset, '_scRef.Rdata')
pred.scRef <- readRDS(rda.scRef)
res.scRef <- simple.evaluation(true.tags, pred.scRef, df.cell.names)
df.sub <- data.frame(term = 'Weighted macro F1', method = 'scMAGIC',
                     value = res.scRef$weighted_macro_f1, stringsAsFactors = F)
df.sub <- rbind(df.sub, 
                data.frame(term = 'Macro F1', method = 'scMAGIC',
                           value = res.scRef$macro_f1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Accuracy', method = 'scMAGIC',
                           value = res.scRef$accuracy, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Balanced accuracy', method = 'scMAGIC',
                           value = res.scRef$balanced.accuracy, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'labeled-Accuracy', method = 'scMAGIC',
                           value = res.scRef$accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'labeled-Balanced accuracy', method = 'scMAGIC',
                           value = res.scRef$balanced.accuracy.rm.unassigned, stringsAsFactors = F))
df.plot <- rbind(df.plot, df.sub)

rda.sciBet <- paste0(path.output, ref.dataset, '_', dataset, '_sciBet.Rdata')
pred.sciBet <- readRDS(rda.sciBet)
res.sciBet <- simple.evaluation(true.tags, pred.sciBet, df.cell.names)
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
df.plot <- rbind(df.plot, df.sub)

rda.singleCellNet <- paste0(path.output, ref.dataset, '_', dataset, '_singleCellNet.Rdata')
pred.singleCellNet <- readRDS(rda.singleCellNet)
res.singleCellNet <- simple.evaluation(true.tags, pred.singleCellNet, df.cell.names)
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
df.plot <- rbind(df.plot, df.sub)

rda.singleR <- paste0(path.output, ref.dataset, '_', dataset, '_singleR.Rdata')
pred.singleR <- readRDS(rda.singleR)
res.singleR <- simple.evaluation(true.tags, pred.singleR, df.cell.names)
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
df.plot <- rbind(df.plot, df.sub)

rda.scmap.cluster <- paste0(path.output, ref.dataset, '_', dataset, '_scmap-cluster.Rdata')
pred.scmap.cluster <- readRDS(rda.scmap.cluster)
res.scmap.cluster <- simple.evaluation(true.tags, pred.scmap.cluster, df.cell.names)
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
df.plot <- rbind(df.plot, df.sub)

rda.scmap.cell <- paste0(path.output, ref.dataset, '_', dataset, '_scmap-cell.Rdata')
pred.scmap.cell <- readRDS(rda.scmap.cell)
res.scmap.cell <- simple.evaluation(true.tags, pred.scmap.cell, df.cell.names)
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
df.plot <- rbind(df.plot, df.sub)

rda.CHETAH <- paste0(path.output, ref.dataset, '_', dataset, '_CHETAH.Rdata')
pred.CHETAH <- readRDS(rda.CHETAH)
res.CHETAH <- simple.evaluation(true.tags, pred.CHETAH, df.cell.names)
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
df.plot <- rbind(df.plot, df.sub)

rda.scPred <- paste0(path.output, ref.dataset, '_', dataset, '_scPred.Rdata')
pred.scPred <- readRDS(rda.scPred)
res.scPred <- simple.evaluation(true.tags, pred.scPred, df.cell.names)
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
df.plot <- rbind(df.plot, df.sub)

rda.scID <- paste0(path.output, ref.dataset, '_', dataset, '_scID.Rdata')
pred.scID <- readRDS(rda.scID)
res.scID <- simple.evaluation(true.tags, pred.scID, df.cell.names)
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
df.plot <- rbind(df.plot, df.sub)

rda.scClassify <- paste0(path.output, ref.dataset, '_', dataset, '_scClassify.Rdata')
pred.scClassify <- readRDS(rda.scClassify)
res.scClassify <- simple.evaluation(true.tags, pred.scClassify, df.cell.names)
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
df.plot <- rbind(df.plot, df.sub)

# sort
df.acc <- df.plot[df.plot$term == 'Accuracy', ]
df.plot$method <- factor(df.plot$method, 
                            levels = df.acc$method[order(df.acc$value, decreasing = T)])

# save results
file.res <- paste0(path.output, 'results_', ref.dataset, '_', dataset, '.txt')
write.table(df.plot, file = file.res, sep = '\t')


library(ggplot2)
plot.bar <- ggplot(df.plot,
       aes(x = method, y = value, fill = term)) +
    geom_bar(position = 'dodge', stat = 'identity') +
    theme_bw() +
    # facet_wrap(~ Evaluator, scales = 'free_x', ncol = 2) +
    labs(title = "", y = 'Score', x = 'Method', fill = 'Metrics methods') + 
    theme(axis.text = element_text(size = 9),
          panel.grid = element_blank(),
          panel.grid.major.y = element_line(color = 'grey', size = 0.2),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(size = 12))
ggsave(filename = paste0(ref.dataset, '_', dataset, '.png'), 
       path = path.output, plot = plot.bar,
       units = 'cm', height = 12, width = 24)

