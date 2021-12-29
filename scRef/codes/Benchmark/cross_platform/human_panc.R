# import python package: sklearn.metrics
library(reticulate)
use_python('/home/zy/tools/anaconda3/bin/python3', required = T)
# py_config()
py_module_available('sklearn')
metrics <- import('sklearn.metrics')

# import data
library(Seurat)
library(SeuratData)
data("panc8")

# evaluation
simple.evaluation <- function(true.tag, scRef.tag, df.cell.names) {
    # uniform tags
    for (j in 1:dim(df.cell.names)[1]) {
        scRef.tag[scRef.tag == df.cell.names[j, 'ref.name']] <- 
            df.cell.names[j, 'sc.name']
    }
    scRef.tag[!(scRef.tag %in% df.cell.names$sc.name)] <- 'Unassigned'
    
    true.tag[true.tag %in% unknow.cell] <- 'Unassigned'
    true.labels <- unique(true.tag)
    our.tag <- scRef.tag
    weighted_macro_f1 <- metrics$f1_score(true.tag, our.tag, average = 'weighted', labels = true.labels)
    macro_f1 <- metrics$f1_score(true.tag, our.tag, average = 'macro', labels = true.labels)
    accuracy <- metrics$accuracy_score(true.tag, our.tag)
    # rm unassigned in tags
    true.tag.rm <- true.tag[our.tag != 'Unassigned']
    our.tag.rm <- our.tag[our.tag != 'Unassigned']
    accuracy.rm.unassigned <- metrics$accuracy_score(true.tag.rm, our.tag.rm)
    
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
    out$f1 <- f1
    out$precision.rm.unassigned <- precision
    out$accuracy.rm.unassigned <- accuracy.rm.unassigned
    out$mean.precision.rm.unassigned <- mean.precision
    out$conf <- table(true.tag, our.tag)
    
    return(out)
    
}

source('/home/zy/my_git/scRef/main/scRef.v20.R')

############# regard sc-counts data as reference
ref.dataset <- 'fluidigmc1'
ref.mtx <- as.matrix(panc8@assays$RNA@counts[, panc8$dataset %in% c('celseq2')])
ref.labels <- as.character(panc8$celltype)[panc8$dataset %in% c('celseq2')]

############### import unlabeled data
############### Habib
dataset <- 'smartseq2'
exp_sc_mat <- as.matrix(panc8@assays$RNA@counts[, panc8$dataset %in% c('smartseq2')])
label_sc <- as.character(panc8$celltype)[panc8$dataset %in% c('smartseq2')]

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
source('/home/zy/my_git/scRef/main/scRef.v20.R')
setwd('~/my_git/scRef')
result.scref <- SCREF(exp_sc_mat, ref.mtx, ref.labels,
                      type_ref = 'sc-counts', use.RUVseq = T, 
                      cluster.speed = F,  out.group = 'HCA',
                      GMM.floor_cutoff = 3, GMM.ceiling_cutoff = 20,
                      min_cell = 1, CPU = 10)
pred.scRef <- result.scref$final.out$scRef.tag
table(label_sc, pred.scRef)
saveRDS(pred.scRef, file = paste0(path.output, ref.dataset, '_', dataset, '_scRef.Rdata'))

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
saveRDS(pred.scmap.cluster, 
        file = paste0(path.output, ref.dataset, '_', dataset, '_scmap-cell.Rdata'))

### CHETAH
library(CHETAH)
library(SingleCellExperiment)
sce <- SingleCellExperiment(assays = list(counts = ref.mtx), 
                            colData = data.frame(celltypes = ref.labels))
sce_test <- SingleCellExperiment(list(counts = exp_sc_mat))
sce_test <- CHETAHclassifier(input = sce_test, ref_cells = sce)
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
df.sub <- data.frame(term = 'Weighted macro F1', method = 'scRef',
                     value = res.scRef$weighted_macro_f1, stringsAsFactors = F)
df.sub <- rbind(df.sub, 
                data.frame(term = 'Macro F1', method = 'scRef',
                           value = res.scRef$macro_f1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Accuracy', method = 'scRef',
                           value = res.scRef$accuracy, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Accuracy (remove unassigned)', method = 'scRef',
                           value = res.scRef$accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Mean precision (remove unassigned)', method = 'scRef',
                           value = res.scRef$mean.precision.rm.unassigned, stringsAsFactors = F))
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
                data.frame(term = 'Accuracy (remove unassigned)', method = 'sciBet',
                           value = res.sciBet$accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Mean precision (remove unassigned)', method = 'sciBet',
                           value = res.sciBet$mean.precision.rm.unassigned, stringsAsFactors = F))
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
                data.frame(term = 'Accuracy (remove unassigned)', method = 'singleCellNet',
                           value = res.singleCellNet$accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Mean precision (remove unassigned)', method = 'singleCellNet',
                           value = res.singleCellNet$mean.precision.rm.unassigned, stringsAsFactors = F))
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
                data.frame(term = 'Accuracy (remove unassigned)', method = 'singleR',
                           value = res.singleR$accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Mean precision (remove unassigned)', method = 'singleR',
                           value = res.singleR$mean.precision.rm.unassigned, stringsAsFactors = F))
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
                data.frame(term = 'Accuracy (remove unassigned)', method = 'scmap-cluster',
                           value = res.scmap.cluster$accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Mean precision (remove unassigned)', method = 'scmap-cluster',
                           value = res.scmap.cluster$mean.precision.rm.unassigned, stringsAsFactors = F))
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
                data.frame(term = 'Accuracy (remove unassigned)', method = 'scmap-cell',
                           value = res.scmap.cell$accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Mean precision (remove unassigned)', method = 'scmap-cell',
                           value = res.scmap.cell$mean.precision.rm.unassigned, stringsAsFactors = F))
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
                data.frame(term = 'Accuracy (remove unassigned)', method = 'CHETAH',
                           value = res.CHETAH$accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Mean precision (remove unassigned)', method = 'CHETAH',
                           value = res.CHETAH$mean.precision.rm.unassigned, stringsAsFactors = F))
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
                data.frame(term = 'Accuracy (remove unassigned)', method = 'scPred',
                           value = res.scPred$accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Mean precision (remove unassigned)', method = 'scPred',
                           value = res.scPred$mean.precision.rm.unassigned, stringsAsFactors = F))
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
                data.frame(term = 'Accuracy (remove unassigned)', method = 'scID',
                           value = res.scID$accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Mean precision (remove unassigned)', method = 'scID',
                           value = res.scID$mean.precision.rm.unassigned, stringsAsFactors = F))
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
                data.frame(term = 'Accuracy (remove unassigned)', method = 'scClassify',
                           value = res.scClassify$accuracy.rm.unassigned, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Mean precision (remove unassigned)', method = 'scClassify',
                           value = res.scClassify$mean.precision.rm.unassigned, stringsAsFactors = F))
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

