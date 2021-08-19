set.seed(1234)
############### Tasic2018
path.output <- '/mdshare/node9/zy/scRef/Benchmark/mouse_brain/'
dataset <- 'Tasic2018'
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
ref.mtx <- OUT$mat_exp
ref.labels <- OUT$label
dim(ref.mtx)
sel.sample <- rownames(ref.labels)[ref.labels[,1] %in% 
                                     c('Endothelial_Endo', 'Endothelial_Peri', 'Endothelial_SMC',
                                       'Non-Neuronal_Astro', 'Non-Neuronal_Macrophage',
                                       'Non-Neuronal_Oligo', 'Non-Neuronal_VLMC')]
col.neu <- setdiff(colnames(ref.mtx), sel.sample)
cells.ref <- c(sel.sample, sample(col.neu,5000-length(sel.sample)))
mat.ref <- ref.mtx[,cells.ref]
tags.ref <- ref.labels[colnames(ref.mtx) %in% cells.ref,]
dim(mat.ref)
true.tags <- tags.ref
true.tags[true.tags %in% c('Endothelial_Endo')] <- 'Endothelial cell'
true.tags[true.tags %in% c('Endothelial_Peri')] <- 'Pericyte'
true.tags[true.tags %in% c('Endothelial_SMC')] <- 'Smooth muscle cell'
true.tags[true.tags %in% c('GABAergic_', 'GABAergic_Lamp5', 'GABAergic_Meis2', 
                           'GABAergic_Pvalb', 'GABAergic_Serpinf1', 'GABAergic_Sncg', 
                           'GABAergic_Sst', 'GABAergic_Vip')] <- 'Neuron'
true.tags[true.tags %in% c('Glutamatergic_', 'Glutamatergic_CR', 'Glutamatergic_L2/3 IT', 
                           'Glutamatergic_L4', 'Glutamatergic_L5 IT', "Glutamatergic_L5 PT", 
                           'Glutamatergic_L6 CT',  'Glutamatergic_L6 IT', 
                           'Glutamatergic_L6b', 'Glutamatergic_NP')] <- 'Neuron'
true.tags[true.tags %in% c('Non-Neuronal_Astro')] <- 'Astrocyte'
true.tags[true.tags %in% c('Non-Neuronal_Macrophage')] <- 'PVM & Microglia'
true.tags[true.tags %in% c('Non-Neuronal_Oligo')] <- 'Oligodendrocyte & OPC'
true.tags[true.tags %in% c('Non-Neuronal_VLMC')] <- 'VLMC'
tags.ref <- true.tags

############### Campbell
path.output <- '/mdshare/node9/zy/scRef/Benchmark/mouse_brain/'
dataset <- 'Campbell'
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
target_mat <- OUT$mat_exp
label_sc <- OUT$label


#### ref 5000 | target 2000/6000/10000/15000/20000
# vec.num <- c(2000, 6000, 9000, 12000, 15000, 20000)
vec.num <- c(15000, 20000)
df.time <- data.frame(stringsAsFactors = F)
for (n.target in vec.num) {
    cells.target <- sample(colnames(target_mat),n.target)
    exp_sc_mat <- target_mat[,cells.target]
    ### scID
    library(scID)
    library(Seurat)
    time1 <- Sys.time()
    Train_Labels <- list(tags.ref)
    names(Train_Labels[[1]]) <- colnames(mat.ref)
    scID_output <- scid_multiclass(exp_sc_mat, mat.ref, Train_Labels[[1]], logFC = 0.3)
    pred.scID <- scID_output$labels
    time2 <- Sys.time()
    time.diff <- difftime(time2, time1, units = 'mins')
    df.time <- rbind(df.time, data.frame(num_cell = n.target, method = 'scID',
                                         run_time = time.diff))

    
    print(df.time)
}

file.runtime <- '/mdshare/node9/zy/scRef/run_time/ref_5000_query_2000_20000_scID.txt'
write.table(df.time, file = file.runtime, sep = '\t', quote = F)

# num_cell method      run_time
# 1     2000   scID 2.044026 mins
# 2     6000   scID 2.331806 mins
# 3     9000   scID 2.615747 mins
# 4    12000   scID 2.888483 mins
# 1    15000   scID 3.773262 mins
# 2    20000   scID 3.987875 mins

# 1     2000   scID 2.044026 
# 2     6000   scID 2.331806 
# 3     9000   scID 2.615747 
# 4    12000   scID 2.888483 
# 1    15000   scID 3.773262 
# 2    20000   scID 3.987875 

# ## scMAGIC
# library(scMAGIC)
# time1 <- Sys.time()
# output.scMAGIC <- scMAGIC(exp_sc_mat, mat.ref, tags.ref, num_threads = 4)
# time2 <- Sys.time()
# time.diff <- difftime(time2, time1, units = 'mins')
# df.time <- rbind(df.time, data.frame(num_cell = n.target, method = 'scMAGIC',
#                                      run_time = time.diff))
# 
# ## sciBet
# time1 <- Sys.time()
# suppressMessages(library(tidyverse))
# suppressMessages(library(scibet))
# suppressMessages(library(viridis))
# suppressMessages(library(ggsci))
# train_set <- as.data.frame(t(mat.ref))
# train_set$label <- tags.ref
# test_set <- as.data.frame(t(exp_sc_mat))
# sciBet <- SciBet(train_set, test_set)
# time2 <- Sys.time()
# time.diff <- difftime(time2, time1, units = 'mins')
# df.time <- rbind(df.time, data.frame(num_cell = n.target, method = 'sciBet',
#                                      run_time = time.diff))
# 
# ### singleCellNet
# time1 <- Sys.time()
# library(singleCellNet)
# library(dplyr)
# out <- get_overlap_genes(exp_sc_mat, mat.ref)
# train_set <- as.matrix(out$exp_ref_mat)
# LabelsTrain <- data.frame(Annotation = tags.ref, row.names = colnames(mat.ref))
# test_set <- as.matrix(out$exp_sc_mat)
# class_info <- scn_train(stTrain = LabelsTrain, expTrain = train_set, dLevel = "Annotation")
# classRes <- scn_predict(cnProc=class_info[['cnProc']], expDat=test_set, nrand = 50)
# classRes <- classRes[, colnames(test_set)]
# tags.singleCellNet <- rownames(classRes)[apply(classRes,2,which.max)]
# tags.singleCellNet[tags.singleCellNet == 'rand'] <- 'Unassigned'
# time2 <- Sys.time()
# time.diff <- difftime(time2, time1, units = 'mins')
# df.time <- rbind(df.time, data.frame(num_cell = n.target, method = 'singleCellNet',
#                                      run_time = time.diff))
# 
# ### singleR
# time1 <- Sys.time()
# library(SingleR)
# train_set <- as.matrix(mat.ref)
# test_set <- as.matrix(exp_sc_mat)
# singler = SingleR(method = "single", sc_data = test_set,
#                   ref_data = train_set,
#                   types = tags.ref, numCores = 4)
# pred.singleR <- singler$labels
# time2 <- Sys.time()
# time.diff <- difftime(time2, time1, units = 'mins')
# df.time <- rbind(df.time, data.frame(num_cell = n.target, method = 'singleR',
#                                      run_time = time.diff))
# 
# ###  scmap
# time1 <- Sys.time()
# library(scmap)
# library(SingleCellExperiment)
# out <- get_overlap_genes(exp_sc_mat, mat.ref)
# train_set <- as.matrix(out$exp_ref_mat)
# test_set <- as.matrix(out$exp_sc_mat)
# sce <- SingleCellExperiment(list(normcounts = train_set),
#                             colData = data.frame(cell_type1 = tags.ref))
# logcounts(sce) <- log2(normcounts(sce) + 1)
# # use gene names as feature symbols
# rowData(sce)$feature_symbol <- rownames(sce)
# sce <- selectFeatures(sce, suppress_plot = TRUE)
# 
# sce_test <- SingleCellExperiment(list(normcounts = test_set))
# logcounts(sce_test) <- log2(normcounts(sce_test) + 1)
# rowData(sce_test)$feature_symbol <- rownames(sce_test)
# sce_test@rowRanges@elementMetadata@listData = sce@rowRanges@elementMetadata@listData
# time2 <- Sys.time()
# time.diff.1 <- difftime(time2, time1, units = 'mins')
# 
# # scmap-cluster
# time1 <- Sys.time()
# sce <- indexCluster(sce)
# scmapCluster_results <- scmapCluster(projection = sce_test,index_list = list(metadata(sce)$scmap_cluster_index))
# pred.scmap.cluster <- scmapCluster_results$combined_labs
# time2 <- Sys.time()
# time.diff.2 <- difftime(time2, time1, units = 'mins')
# df.time <- rbind(df.time, data.frame(num_cell = n.target, method = 'scmap-cluster',
#                                      run_time = time.diff.1 + time.diff.2))
# 
# # scmap-cell
# time1 <- Sys.time()
# set.seed(1)
# sce <- indexCell(sce)
# scmapCell_results <- scmapCell(sce_test,list(metadata(sce)$scmap_cell_index))
# scmapCell_clusters <- scmapCell2Cluster(scmapCell_results,list(as.character(colData(sce)$cell_type1)))
# pred.scmap.cell <- scmapCell_clusters$combined_labs
# time2 <- Sys.time()
# time.diff.3 <- difftime(time2, time1, units = 'mins')
# df.time <- rbind(df.time, data.frame(num_cell = n.target, method = 'scmap-cell',
#                                      run_time = time.diff.1 + time.diff.3))
# 
# ### CHETAH
# time1 <- Sys.time()
# library(CHETAH)
# library(SingleCellExperiment)
# sce <- SingleCellExperiment(assays = list(counts = mat.ref),
#                             colData = data.frame(celltypes = tags.ref))
# sce_test <- SingleCellExperiment(list(counts = exp_sc_mat))
# # sce_test <- CHETAHclassifier(input = sce_test, ref_cells = sce)
# sce_test <- CHETAHclassifier(input = sce_test, ref_cells = sce,
#                              n_genes = median(colSums(mat.ref != 0))/2)
# pred.CHETAH <- sce_test$celltype_CHETAH
# time2 <- Sys.time()
# time.diff <- difftime(time2, time1, units = 'mins')
# df.time <- rbind(df.time, data.frame(num_cell = n.target, method = 'CHETAH',
#                                      run_time = time.diff))
# 
# ### scPred
# time1 <- Sys.time()
# library("scPred")
# library("Seurat")
# library("magrittr")
# # scPred Training
# reference <- CreateSeuratObject(counts = mat.ref)
# reference@meta.data$cell_type <- tags.ref
# reference <- reference %>%
#     NormalizeData() %>%
#     FindVariableFeatures() %>%
#     ScaleData() %>%
#     RunPCA() %>%
#     RunUMAP(dims = 1:30)
# reference <- getFeatureSpace(reference, "cell_type")
# reference <- trainModel(reference, allowParallel = T)
# # scPred Prediction
# query <- CreateSeuratObject(counts = exp_sc_mat)
# query <- NormalizeData(query)
# query <- scPredict(query, reference)
# pred.scPred <- query@meta.data$scpred_prediction
# time2 <- Sys.time()
# time.diff <- difftime(time2, time1, units = 'mins')
# df.time <- rbind(df.time, data.frame(num_cell = n.target, method = 'scPred',
#                                      run_time = time.diff))
# 
# ### scClassify
# time1 <- Sys.time()
# library("scClassify")
# library(Matrix)
# exprsMat_train <- as(as.matrix(log1p(mat.ref)), "dgCMatrix")
# exp_sc_mat <- as(as.matrix(log1p(exp_sc_mat)), "dgCMatrix")
# scClassify_res <- scClassify(exprsMat_train = exprsMat_train,
#                              cellTypes_train = tags.ref,
#                              exprsMat_test = list(one = exp_sc_mat),
#                              tree = "HOPACH",
#                              algorithm = "WKNN",
#                              selectFeatures = c("limma"),
#                              similarity = c("pearson"),
#                              returnList = FALSE,parallel = T,
#                              verbose = FALSE)
# pred.scClassify <- scClassify_res$testRes$one$pearson_WKNN_limma$predRes
# time2 <- Sys.time()
# time.diff <- difftime(time2, time1, units = 'mins')
# df.time <- rbind(df.time, data.frame(num_cell = n.target, method = 'scClassify',
#                                      run_time = time.diff))
