run_SingleR<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,GeneOrderPath = NULL,NumGenes = NULL){
    "
  run SingleR
  Wrapper script to run SingleR on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.tsv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.tsv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
    
    Data <- read.delim(DataPath,row.names = 1)
    Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
    load(CV_RDataPath)
    Labels <- as.vector(Labels[,col_Index])
    Data <- Data[Cells_to_Keep,]
    Labels <- Labels[Cells_to_Keep]
    if(!is.null(GeneOrderPath) & !is.null(NumGenes)){
        GenesOrder = read.csv(GeneOrderPath)
    }
    
    #############################################################################
    #                               SingleR                                     #
    #############################################################################
    library(SingleR)
    library(Seurat)
    True_Labels_SingleR <- list()
    Pred_Labels_SingleR <- list()
    Total_Time_SingleR <- list()

    for (i in c(1:n_folds)){
      train_data <- as.matrix(Data[,Train_Idx[[i]]])
      train_label <- Labels[Train_Idx[[i]]]
      test_set <- as.matrix(Data[,Test_Idx[[i]]])
      if(!is.null(GeneOrderPath) & !is.null(NumGenes)){
            start_time <- Sys.time()
            singler = SingleR(method = "single", Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]], 
                              Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]], 
                              Labels[Train_Idx[[i]]], numCores = 1)
            end_time <- Sys.time()
        }
        else{
            start_time <- Sys.time()
            singler = SingleR(method = "single", sc_data = test_set, 
                              ref_data = train_data, 
                              types = train_label, numCores = 1)
            end_time <- Sys.time()
        }
        Total_Time_SingleR[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
        
        True_Labels_SingleR[i] <- list(Labels[Test_Idx[[i]]])
        Pred_Labels_SingleR[i] <- list(as.vector(singler$labels))
    }
    True_Labels_SingleR <- as.vector(unlist(True_Labels_SingleR))
    Pred_Labels_SingleR <- as.vector(unlist(Pred_Labels_SingleR))
    Total_Time_SingleR <- as.vector(unlist(Total_Time_SingleR))
    
    setwd(OutputDir)
    
    if(!is.null(GeneOrderPath) & !is.null(NumGenes)){
        write.csv(True_Labels_SingleR,paste('SingleR_',NumGenes,'_True_Labels.csv', sep = ''),row.names = FALSE)
        write.csv(Pred_Labels_SingleR,paste('SingleR_',NumGenes,'_Pred_Labels.csv', sep = ''),row.names = FALSE)
        write.csv(Total_Time_SingleR,paste('SingleR_',NumGenes,'_Total_Time.csv', sep = ''),row.names = FALSE)
    }
    else{
        write.csv(True_Labels_SingleR,'SingleR_True_Labels.csv',row.names = FALSE)
        write.csv(Pred_Labels_SingleR,'SingleR_Pred_Labels.csv',row.names = FALSE)
        write.csv(Total_Time_SingleR,'SingleR_Total_Time.csv',row.names = FALSE)
    }
}


run_scmap <- function(DataPath,LabelsPath,CV_RDataPath,OutputDir,
                      GeneOrderPath = NULL,NumGenes = NULL){
  "
  run scmap
  Wrapper script to run scmap on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  Data <- read.delim(DataPath,row.names = 1)
  Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                                 scmap                                     #
  #############################################################################
  library(scmap)
  library(SingleCellExperiment)
  True_Labels_scmapcluster <- list()
  Pred_Labels_scmapcluster <- list()
  True_Labels_scmapcell <- list()
  Pred_Labels_scmapcell <- list()
  Training_Time_scmapcluster <- list()
  Testing_Time_scmapcluster <- list()
  Training_Time_scmapcell <- list()
  Testing_Time_scmapcell <- list()
  Data = as.matrix(Data)
  
  for (i in c(1:n_folds)){
    if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
      sce <- SingleCellExperiment(list(normcounts = Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]]), 
                                  colData = data.frame(cell_type1 = Labels[Train_Idx[[i]]]))
      logcounts(sce) <- log2(normcounts(sce) + 1)
      # use gene names as feature symbols
      rowData(sce)$feature_symbol <- rownames(sce)
      sce <- selectFeatures(sce, n_features = NumGenes, suppress_plot = TRUE)
      
      sce_test <- SingleCellExperiment(list(normcounts = Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]]), 
                                       colData = data.frame(cell_type1 = Labels[Test_Idx[[i]]]))
      logcounts(sce_test) <- log2(normcounts(sce_test) + 1)
      rowData(sce_test)$feature_symbol <- rownames(sce_test)
      sce_test@rowRanges@elementMetadata@listData = sce@rowRanges@elementMetadata@listData
    }
    else{
      sce <- SingleCellExperiment(list(normcounts = Data[,Train_Idx[[i]]]), 
                                  colData = data.frame(cell_type1 = Labels[Train_Idx[[i]]]))
      logcounts(sce) <- log2(normcounts(sce) + 1)
      # use gene names as feature symbols
      rowData(sce)$feature_symbol <- rownames(sce)
      sce <- selectFeatures(sce, suppress_plot = TRUE)
      
      sce_test <- SingleCellExperiment(list(normcounts = Data[,Test_Idx[[i]]]), 
                                       colData = data.frame(cell_type1 = Labels[Test_Idx[[i]]]))
      logcounts(sce_test) <- log2(normcounts(sce_test) + 1)
      rowData(sce_test)$feature_symbol <- rownames(sce_test)
      sce_test@rowRanges@elementMetadata@listData = sce@rowRanges@elementMetadata@listData
    }
    
    # scmap-cluster
    start_time <- Sys.time()
    sce <- indexCluster(sce)
    end_time <- Sys.time()
    Training_Time_scmapcluster[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    start_time <- Sys.time()
    scmapCluster_results <- scmapCluster(projection = sce_test,index_list = list(metadata(sce)$scmap_cluster_index))
    end_time <- Sys.time()
    Testing_Time_scmapcluster[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    True_Labels_scmapcluster[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_scmapcluster[i] <- list(scmapCluster_results$combined_labs)
    
    # scmap-cell
    start_time <- Sys.time()
    set.seed(1)
    sce <- indexCell(sce)
    end_time <- Sys.time()
    Training_Time_scmapcell[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    start_time <- Sys.time()
    scmapCell_results <- scmapCell(sce_test,list(metadata(sce)$scmap_cell_index))
    scmapCell_clusters <- scmapCell2Cluster(scmapCell_results,list(as.character(colData(sce)$cell_type1)))
    end_time <- Sys.time()
    Testing_Time_scmapcell[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    True_Labels_scmapcell[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_scmapcell[i] <- list(scmapCell_clusters$combined_labs)
  }
  
  True_Labels_scmapcluster <- as.vector(unlist(True_Labels_scmapcluster))
  Pred_Labels_scmapcluster <- as.vector(unlist(Pred_Labels_scmapcluster))
  True_Labels_scmapcell <- as.vector(unlist(True_Labels_scmapcell))
  Pred_Labels_scmapcell <- as.vector(unlist(Pred_Labels_scmapcell))
  Training_Time_scmapcluster <- as.vector(unlist(Training_Time_scmapcluster))
  Testing_Time_scmapcluster <- as.vector(unlist(Testing_Time_scmapcluster))
  Training_Time_scmapcell <- as.vector(unlist(Training_Time_scmapcell))
  Testing_Time_scmapcell <- as.vector(unlist(Testing_Time_scmapcell))
  
  setwd(OutputDir)
  
  if (!is.null(GeneOrderPath) & !is.null (NumGenes)){
    write.csv(True_Labels_scmapcluster,paste('scmapcluster_',NumGenes,'_True_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Pred_Labels_scmapcluster,paste('scmapcluster_',NumGenes,'_Pred_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(True_Labels_scmapcell,paste('scmapcell_',NumGenes,'_True_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Pred_Labels_scmapcell,paste('scmapcell_',NumGenes,'_Pred_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Training_Time_scmapcluster,paste('scmapcluster_',NumGenes,'_Training_Time.csv', sep = ''),row.names = FALSE)
    write.csv(Testing_Time_scmapcluster,paste('scmapcluster_',NumGenes,'_Testing_Time.csv', sep = ''),row.names = FALSE)
    write.csv(Training_Time_scmapcell,paste('scmapcell_',NumGenes,'_Training_Time.csv', sep = ''),row.names = FALSE)
    write.csv(Testing_Time_scmapcell,paste('scmapcell_',NumGenes,'_Testing_Time.csv', sep = ''),row.names = FALSE)
  }
  else{
    write.csv(True_Labels_scmapcluster,'scmapcluster_True_Labels.csv',row.names = FALSE)
    write.csv(Pred_Labels_scmapcluster,'scmapcluster_Pred_Labels.csv',row.names = FALSE)
    write.csv(True_Labels_scmapcell,'scmapcell_True_Labels.csv',row.names = FALSE)
    write.csv(Pred_Labels_scmapcell,'scmapcell_Pred_Labels.csv',row.names = FALSE)
    write.csv(Training_Time_scmapcluster,'scmapcluster_Training_Time.csv',row.names = FALSE)
    write.csv(Testing_Time_scmapcluster,'scmapcluster_Testing_Time.csv',row.names = FALSE)
    write.csv(Training_Time_scmapcell,'scmapcell_Training_Time.csv',row.names = FALSE)
    write.csv(Testing_Time_scmapcell,'scmapcell_Testing_Time.csv',row.names = FALSE)
  }
}


run_CHETAH<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,GeneOrderPath = NULL,NumGenes = NULL){
  "
  run CHETAH
  Wrapper script to run CHETAH on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  Data <- read.delim(DataPath,row.names = 1)
  Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                                CHETAH                                     #
  #############################################################################
  library(CHETAH)
  library(SingleCellExperiment)
  True_Labels_CHETAH <- list()
  Pred_Labels_CHETAH <- list()
  Total_Time_CHETAH <- list()
  Data = as.matrix(Data)
  
  for (i in c(1:n_folds)){
    if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
      sce <- SingleCellExperiment(assays = list(counts = Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]]), 
                                  colData = data.frame(celltypes = Labels[Train_Idx[[i]]]))
      
      sce_test <- SingleCellExperiment(assays = list(counts = Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]]), 
                                       colData = data.frame(celltypes = Labels[Test_Idx[[i]]]))
      start_time <- Sys.time()
      sce_test <- CHETAHclassifier(input = sce_test, ref_cells = sce, n_genes = NumGenes)
      end_time <- Sys.time()
    }
    else{
      sce <- SingleCellExperiment(assays = list(counts = Data[,Train_Idx[[i]]]), 
                                  colData = data.frame(celltypes = Labels[Train_Idx[[i]]]))
      
      sce_test <- SingleCellExperiment(assays = list(counts = Data[,Test_Idx[[i]]]), 
                                       colData = data.frame(celltypes = Labels[Test_Idx[[i]]]))
      start_time <- Sys.time()
      sce_test <- CHETAHclassifier(input = sce_test, ref_cells = sce)
      end_time <- Sys.time()
    }
    
    Total_Time_CHETAH[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    True_Labels_CHETAH[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_CHETAH[i] <- list(sce_test$celltype_CHETAH)
  }
  True_Labels_CHETAH <- as.vector(unlist(True_Labels_CHETAH))
  Pred_Labels_CHETAH <- as.vector(unlist(Pred_Labels_CHETAH))
  Total_Time_CHETAH <- as.vector(unlist(Total_Time_CHETAH))
  
  setwd(OutputDir)
  
  if (!is.null(GeneOrderPath) & !is.null (NumGenes)){
    write.csv(True_Labels_CHETAH,paste('CHETAH_',NumGenes,'_True_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Pred_Labels_CHETAH,paste('CHETAH_',NumGenes,'_Pred_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Total_Time_CHETAH,paste('CHETAH_',NumGenes,'_Total_Time.csv', sep = ''),row.names = FALSE)
  }
  else{
    write.csv(True_Labels_CHETAH,'CHETAH_True_Labels.csv',row.names = FALSE)
    write.csv(Pred_Labels_CHETAH,'CHETAH_Pred_Labels.csv',row.names = FALSE)
    write.csv(Total_Time_CHETAH,'CHETAH_Total_Time.csv',row.names = FALSE)
  }
}


run_scPred<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,
                     GeneOrderPath = NULL,NumGenes = NULL){
  "
  run scPred
  Wrapper script to run scPred on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  Data <- read.delim(DataPath,row.names = 1)
  Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  print(dim(Data))
  print(length(Labels))
  
  Labels <- Labels[Cells_to_Keep]
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                                scPred                                     #
  #############################################################################
  library(scPred)
  library(tidyverse)
  library(SingleCellExperiment)
  True_Labels_scPred <- list()
  Pred_Labels_scPred <- list()
  Training_Time_scPred <- list()
  Testing_Time_scPred <- list()
  Data = as.matrix(round(Data))
  
  for (i in c(1:n_folds)){
    if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
      sce <- SingleCellExperiment(list(normcounts = Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]]), 
                                  colData = data.frame(cell_type1 = Labels[Train_Idx[[i]]]))
      sce_counts <- normcounts(sce)
      sce_cpm <- apply(sce_counts, 2, function(x) (x/sum(x))*1000000)
      sce_metadata <- as.data.frame(colData(sce))
      
      sce_test <- SingleCellExperiment(list(normcounts = Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]]), 
                                       colData = data.frame(cell_type1 = Labels[Test_Idx[[i]]]))
      sce_counts_test <- normcounts(sce_test)
      sce_cpm_test <- apply(sce_counts_test, 2, function(x) (x/sum(x))*1000000)
      sce_metadata_test <- as.data.frame(colData(sce_test))
    }else{
      sce <- SingleCellExperiment(list(normcounts = Data[,Train_Idx[[i]]]), 
                                  colData = data.frame(cell_type1 = Labels[Train_Idx[[i]]]))
      sce_counts <- normcounts(sce)
      sce_cpm <- apply(sce_counts, 2, function(x) (x/sum(x))*1000000)
      sce_metadata <- as.data.frame(colData(sce))
      
      sce_test <- SingleCellExperiment(list(normcounts = Data[,Test_Idx[[i]]]), 
                                       colData = data.frame(cell_type1 = Labels[Test_Idx[[i]]]))
      sce_counts_test <- normcounts(sce_test)
      sce_cpm_test <- apply(sce_counts_test, 2, function(x) (x/sum(x))*1000000)
      sce_metadata_test <- as.data.frame(colData(sce_test))
    }
    
    
    # scPred Training    
    start_time <- Sys.time()
    set.seed(1234)
    scp <- eigenDecompose(sce_cpm)
    scPred::metadata(scp) <- sce_metadata
    scp <- getFeatureSpace(scp, pVar = 'cell_type1')
    # plotEigen(scp, group = 'cell_type1')
    scp <- trainModel(scp)
    # plotTrainProbs(scp)
    end_time <- Sys.time()
    Training_Time_scPred[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    # scPred Prediction
    start_time <- Sys.time()
    scp <- scPredict(scp,newData = sce_cpm_test)
    end_time <- Sys.time()
    Testing_Time_scPred[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    True_Labels_scPred[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_scPred[i] <- list(getPredictions(scp)$predClass)
  }
  True_Labels_scPred <- as.vector(unlist(True_Labels_scPred))
  Pred_Labels_scPred <- as.vector(unlist(Pred_Labels_scPred))
  Training_Time_scPred <- as.vector(unlist(Training_Time_scPred))
  Testing_Time_scPred <- as.vector(unlist(Testing_Time_scPred))
  Pred_Labels_scPred <- gsub('pyramidal.SS','pyramidal SS', Pred_Labels_scPred)
  Pred_Labels_scPred <- gsub('pyramidal.CA1','pyramidal CA1', Pred_Labels_scPred)
  Pred_Labels_scPred <- gsub('endothelial.mural','endothelial-mural', Pred_Labels_scPred)
  
  setwd(OutputDir)
  
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    write.csv(True_Labels_scPred,paste('scPred_',NumGenes,'_True_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Pred_Labels_scPred,paste('scPred_',NumGenes,'_Pred_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Training_Time_scPred,paste('scPred_',NumGenes,'_Training_Time.csv', sep = ''),row.names = FALSE)
    write.csv(Testing_Time_scPred,paste('scPred_',NumGenes,'_Testing_Time.csv', sep = ''),row.names = FALSE)
  }
  else{
    write.csv(True_Labels_scPred,'scPred_True_Labels.csv',row.names = FALSE)
    write.csv(Pred_Labels_scPred,'scPred_Pred_Labels.csv',row.names = FALSE)
    write.csv(Training_Time_scPred,'scPred_Training_Time.csv',row.names = FALSE)
    write.csv(Testing_Time_scPred,'scPred_Testing_Time.csv',row.names = FALSE)
  }
}


run_sciBet<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,
                     GeneOrderPath = NULL,NumGenes = NULL){
  "
  run sciBet
  Wrapper script to run sciBet on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.tsv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.tsv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  Data <- read.delim(DataPath,row.names = 1)
  Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
  load(CV_RDataPath)
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  if(!is.null(GeneOrderPath) & !is.null(NumGenes)){
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                               sciBet                                     #
  #############################################################################
  suppressMessages(library(tidyverse))
  suppressMessages(library(scibet))
  suppressMessages(library(viridis))
  suppressMessages(library(ggsci))
  True_Labels_sciBet <- list()
  Pred_Labels_sciBet <- list()
  Total_Time_sciBet <- list()
  # Data = t(as.matrix(Data))
  # names(Data) <- 1:dim(Data)[2]
  Data <- Data / 1.0
  
  for (i in c(1:n_folds)){
    train_set <- as.data.frame(t(Data[,Train_Idx[[i]]]))
    train_set$label <- Labels[Train_Idx[[i]]]
    test_set <- as.data.frame(t(Data[,Test_Idx[[i]]]))
    if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
      start_time <- Sys.time()
      # sciBet = sciBet(method = "single", Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]], 
      #                   Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]], 
      #                   Labels[Train_Idx[[i]]], numCores = 1)
      end_time <- Sys.time()
    } else {
      start_time <- Sys.time()
      sciBet = SciBet(train_set, test_set)
      end_time <- Sys.time()
    }
    Total_Time_sciBet[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    True_Labels_sciBet[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_sciBet[i] <- list(sciBet)
  }
  True_Labels_sciBet <- as.vector(unlist(True_Labels_sciBet))
  Pred_Labels_sciBet <- as.vector(unlist(Pred_Labels_sciBet))
  Total_Time_sciBet <- as.vector(unlist(Total_Time_sciBet))
  
  setwd(OutputDir)
  
  if(!is.null(GeneOrderPath) & !is.null(NumGenes)){
    write.csv(True_Labels_sciBet,paste('sciBet_',NumGenes,'_True_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Pred_Labels_sciBet,paste('sciBet_',NumGenes,'_Pred_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Total_Time_sciBet,paste('sciBet_',NumGenes,'_Total_Time.csv', sep = ''),row.names = FALSE)
  }
  else{
    write.csv(True_Labels_sciBet,'sciBet_True_Labels.csv',row.names = FALSE)
    write.csv(Pred_Labels_sciBet,'sciBet_Pred_Labels.csv',row.names = FALSE)
    write.csv(Total_Time_sciBet,'sciBet_Total_Time.csv',row.names = FALSE)
  }
}


run_scRef<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,
                    GeneOrderPath = NULL,NumGenes = NULL){
  "
  run scRef
  Wrapper script to run scPred on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.tsv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.tsv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  Data <- read.delim(DataPath,row.names = 1)
  Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
  load(CV_RDataPath)
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  if(!is.null(GeneOrderPath) & !is.null(NumGenes)){
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                               scRef                                     #
  #############################################################################
  # source('/home/zy/my_git/scRef/main/scRef.v8.R')
  source('/home/drizzle_zhang/my_git/scRef/main/scRef.v8.R')
  True_Labels_scRef <- list()
  Pred_Labels_scRef <- list()
  Total_Time_scRef <- list()

  for (i in c(1:n_folds)) {
    train_data <- Data[,Train_Idx[[i]]]
    train_label <- data.frame(cell_id = names(train_data), tag = Labels[Train_Idx[[i]]])
    train_set <- .generate_ref(train_data, train_label)
    test_set <- Data[,Test_Idx[[i]]]
    if (!is.null(GeneOrderPath) & !is.null(NumGenes)) {
      start_time <- Sys.time()
      # scRef = scRef(method = "single", Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]], 
      #                   Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]], 
      #                   Labels[Train_Idx[[i]]], numCores = 1)
      end_time <- Sys.time()
    } else {
      start_time <- Sys.time()
      setwd('~/my_git/scRef')
      # scRef = SCREF(test_set, train_set, type_ref = 'count',
      #               cluster.speed = T, cluster.cell = 3, min_cell = 3, CPU = 4)
      scRef = SCREF(test_set, train_set, type_ref = 'count',
                    cluster.num.pc = 50, cluster.resolution = 1, 
                    cluster.speed = F, min_cell = 1, CPU = 6)
      # scRef = SCREF(test_set, train_set, 
      #               identify_unassigned = F, CPU = 6)
      label.scRef <- as.character(scRef$final.out$scRef.tag)
      end_time <- Sys.time()
    }
    Total_Time_scRef[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    True_Labels_scRef[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_scRef[i] <- list(label.scRef)
  }
  True_Labels_scRef <- as.vector(unlist(True_Labels_scRef))
  Pred_Labels_scRef <- as.vector(unlist(Pred_Labels_scRef))
  Total_Time_scRef <- as.vector(unlist(Total_Time_scRef))
  
  setwd(OutputDir)
  
  if(!is.null(GeneOrderPath) & !is.null(NumGenes)){
    write.csv(True_Labels_scRef,paste('scRef_',NumGenes,'_True_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Pred_Labels_scRef,paste('scRef_',NumGenes,'_Pred_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Total_Time_scRef,paste('scRef_',NumGenes,'_Total_Time.csv', sep = ''),row.names = FALSE)
  } else{
    write.csv(True_Labels_scRef,'scRef_True_Labels.csv',row.names = FALSE)
    write.csv(Pred_Labels_scRef,'scRef_Pred_Labels.csv',row.names = FALSE)
    write.csv(Total_Time_scRef,'scRef_Total_Time.csv',row.names = FALSE)
  }
}


run_singleCellNet<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,
                            GeneOrderPath = NULL,NumGenes = NULL){
  "
  run singleCellNet
  Wrapper script to run singleCellNet on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  Data <- read.delim(DataPath,row.names = 1)
  genes <- rownames(Data)
  genes <- gsub('_', '.', genes)
  rownames(Data) <- genes
  Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                              singleCellNet                                #
  #############################################################################
  library(singleCellNet)
  library(dplyr)
  True_Labels_singleCellNet <- list()
  Pred_Labels_singleCellNet <- list()
  Training_Time_singleCellNet <- list()
  Testing_Time_singleCellNet <- list()
  # Data = as.matrix(Data)            # deals also with sparse matrix
  
  for(i in c(1:n_folds)){
    if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
      DataTrain <- Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]]
      DataTest <- Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]]
    } else{
      DataTrain <- as.matrix(Data[,Train_Idx[[i]]])
      LabelsTrain <- data.frame(Annotation = Labels[Train_Idx[[i]]], row.names = colnames(DataTrain))
      DataTest <- as.matrix(Data[,Test_Idx[[i]]])
    }
    
    start_time <- Sys.time()
    class_info <- scn_train(stTrain = LabelsTrain, expTrain = DataTrain, dLevel = "Annotation")
    end_time <- Sys.time()
    Training_Time_singleCellNet[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    start_time <- Sys.time()
    classRes <-scn_predict(cnProc=class_info[['cnProc']], expDat=DataTest, nrand = 50)
    end_time <- Sys.time()
    Testing_Time_singleCellNet[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    True_Labels_singleCellNet[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_singleCellNet[i] <- list((rownames(classRes)[apply(classRes,2,which.max)])[1:length(Test_Idx[[i]])])
  }
  True_Labels_singleCellNet <- as.vector(unlist(True_Labels_singleCellNet))
  Pred_Labels_singleCellNet <- as.vector(unlist(Pred_Labels_singleCellNet))
  Training_Time_singleCellNet <- as.vector(unlist(Training_Time_singleCellNet))
  Testing_Time_singleCellNet <- as.vector(unlist(Testing_Time_singleCellNet))
  
  setwd(OutputDir)
  
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    write.csv(True_Labels_singleCellNet,paste('singleCellNet_',NumGenes,'_True_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Pred_Labels_singleCellNet,paste('singleCellNet_',NumGenes,'_Pred_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Training_Time_singleCellNet,paste('singleCellNet_',NumGenes,'_Training_Time.csv', sep = ''),row.names = FALSE)
    write.csv(Testing_Time_singleCellNet,paste('singleCellNet_',NumGenes,'_Testing_Time.csv', sep = ''),row.names = FALSE)
  }
  else{
    write.csv(True_Labels_singleCellNet,'singleCellNet_True_Labels.csv',row.names = FALSE)
    write.csv(Pred_Labels_singleCellNet,'singleCellNet_Pred_Labels.csv',row.names = FALSE)
    write.csv(Training_Time_singleCellNet,'singleCellNet_Training_Time.csv',row.names = FALSE)
    write.csv(Testing_Time_singleCellNet,'singleCellNet_Testing_Time.csv',row.names = FALSE)
  }
}


run_CaSTLe<-function(DataPath,LabelsPath,CV_RDataPath, OutputDir, GeneOrderPath = NULL, NumGenes = NULL){
  "
  run CaSTLe
  Wrapper script to run CaSTLe on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  Data <- read.delim(DataPath,row.names = 1)
  Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                                CaSTLe                                     #
  #############################################################################
  library(igraph)
  library(xgboost)
  True_Labels_Castle <- list()
  Pred_Labels_Castle <- list()
  Training_Time_Castle <- list()
  Testing_Time_Castle <- list()
  Data <- t(as.matrix(Data))
  
  BREAKS=c(-1, 0, 1, 6, Inf)
  nFeatures = 100
  
  for(i in c(1:n_folds)){
    # 1. Load datasets
    if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
      ds1 = Data[Train_Idx[[i]],as.vector(GenesOrder[c(1:NumGenes),i])+1]
      ds2 = Data[Test_Idx[[i]],as.vector(GenesOrder[c(1:NumGenes),i])+1]
    }
    else{
      ds1 = Data[Train_Idx[[i]],]
      ds2 = Data[Test_Idx[[i]],]
    }
    
    sourceCellTypes = as.factor(Labels[Train_Idx[[i]]])
    targetCellTypes = as.factor(Labels[Test_Idx[[i]]])
    
    start_time <- Sys.time()
    # 2. Unify sets, excluding low expressed genes
    source_n_cells_counts = apply(ds1, 2, function(x) { sum(x > 0) } )
    target_n_cells_counts = apply(ds2, 2, function(x) { sum(x > 0) } )
    common_genes = intersect( colnames(ds1)[source_n_cells_counts>10], 
                              colnames(ds2)[target_n_cells_counts>10])
    remove(source_n_cells_counts, target_n_cells_counts)
    ds1 = ds1[, colnames(ds1) %in% common_genes]
    ds2 = ds2[, colnames(ds2) %in% common_genes]
    ds = rbind(ds1[,common_genes], ds2[,common_genes])
    isSource = c(rep(TRUE,nrow(ds1)), rep(FALSE,nrow(ds2)))
    remove(ds1, ds2)
    
    # 3. Highest mean in both source and target
    topFeaturesAvg = colnames(ds)[order(apply(ds, 2, mean), decreasing = T)]
    end_time <- Sys.time()
    Training_Time_Castle[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    start_time <- Sys.time()
    # for each cell - what is the most probable classification?
    L = length(levels(sourceCellTypes))
    targetClassification = as.data.frame(matrix(rep(0,L*sum(!isSource)), nrow=L), row.names = levels(sourceCellTypes))
    
    for (cellType in levels(sourceCellTypes)) {
      
      inSourceCellType = as.factor(ifelse(sourceCellTypes == cellType, cellType, paste0("NOT",cellType)))
      
      # 4. Highest mutual information in source
      topFeaturesMi = names(sort(apply(ds[isSource,],2,function(x) { compare(cut(x,breaks=BREAKS),inSourceCellType,method = "nmi") }), decreasing = T))
      
      # 5. Top n genes that appear in both mi and avg
      selectedFeatures = union(head(topFeaturesAvg, nFeatures) , head(topFeaturesMi, nFeatures) )
      
      # 6. remove correlated features
      tmp = cor(ds[,selectedFeatures], method = "pearson")
      tmp[!lower.tri(tmp)] = 0
      selectedFeatures = selectedFeatures[apply(tmp,2,function(x) any(x < 0.9))]
      remove(tmp)
      
      # 7,8. Convert data from continous to binned dummy vars
      # break datasets to bins
      dsBins = apply(ds[, selectedFeatures], 2, cut, breaks= BREAKS)
      # use only bins with more than one value
      nUniq = apply(dsBins, 2, function(x) { length(unique(x)) })
      # convert to dummy vars
      ds0 = model.matrix(~ . , as.data.frame(dsBins[,nUniq>1]))
      remove(dsBins, nUniq)
      
      cat(paste0("<h2>Classifier for ",cellType,"</h2>"))
      
      inTypeSource = sourceCellTypes == cellType
      # 9. Classify
      xg=xgboost(data=ds0[isSource,] , 
                 label=inTypeSource,
                 objective="binary:logistic", 
                 eta=0.7 , nthread=1, nround=20, verbose=0,
                 gamma=0.001, max_depth=5, min_child_weight=10)
      
      # 10. Predict
      inTypeProb = predict(xg, ds0[!isSource, ])
      
      targetClassification[cellType,] = inTypeProb
    }
    end_time <- Sys.time()
    Testing_Time_Castle[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    True_Labels_Castle[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_Castle[i] <- list(rownames(targetClassification)[apply(targetClassification,2,which.max)])
  }
  True_Labels_Castle <- as.vector(unlist(True_Labels_Castle))
  Pred_Labels_Castle <- as.vector(unlist(Pred_Labels_Castle))
  Training_Time_Castle <- as.vector(unlist(Training_Time_Castle))
  Testing_Time_Castle <- as.vector(unlist(Testing_Time_Castle))
  
  setwd(OutputDir)
  
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    write.csv(True_Labels_Castle,paste('True_Labels_Castle_',NumGenes,'.csv', sep = ''),row.names = FALSE)
    write.csv(Pred_Labels_Castle,paste('Pred_Labels_Castle_',NumGenes,'.csv', sep = ''),row.names = FALSE)
    write.csv(Training_Time_Castle,paste('Training_Time_Castle_',NumGenes,'.csv', sep = ''),row.names = FALSE)
    write.csv(Testing_Time_Castle,paste('Testing_Time_Castle_',NumGenes,'.csv', sep = ''),row.names = FALSE)
  }
  else{
    write.csv(True_Labels_Castle,'CaSTLe_True_Labels.csv',row.names = FALSE)
    write.csv(Pred_Labels_Castle,'CaSTLe_Pred_Labels.csv',row.names = FALSE)
    write.csv(Training_Time_Castle,'CaSTLe_Training_Time.csv',row.names = FALSE)
    write.csv(Testing_Time_Castle,'CaSTLe_Testing_Time.csv',row.names = FALSE)
  }
  
}


run_scID<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,GeneOrderPath = NULL,NumGenes = NULL){
  "
  run scID
  Wrapper script to run scID on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  Data <- read.delim(DataPath,row.names = 1)
  Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                                 scID                                      #
  #############################################################################
  library(scID)
  library(Seurat)
  True_Labels_scID <- list()
  Pred_Labels_scID <- list()
  Total_Time_scID <- list()
  Data = as.matrix(Data)
  
  for (i in c(1:n_folds)){
    if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
      Train_Labels <- list(Labels[Train_Idx[[i]]])
      names(Train_Labels[[1]]) <- colnames(Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]])
      start_time <- Sys.time()
      scID_output <- scid_multiclass(Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]], 
                                     Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]], 
                                     Train_Labels[[1]])
      end_time <- Sys.time()
    }
    else{
      Train_Labels <- list(Labels[Train_Idx[[i]]])
      names(Train_Labels[[1]]) <- colnames(Data[,Train_Idx[[i]]])
      start_time <- Sys.time()
      scID_output <- scid_multiclass(Data[,Test_Idx[[i]]], Data[,Train_Idx[[i]]], Train_Labels[[1]])
      end_time <- Sys.time()
    }
    Total_Time_scID[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    True_Labels_scID[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_scID[i] <- list(as.vector(scID_output$labels))
  }
  True_Labels_scID <- as.vector(unlist(True_Labels_scID))
  Pred_Labels_scID <- as.vector(unlist(Pred_Labels_scID))
  Total_Time_scID <- as.vector(unlist(Total_Time_scID))
  
  setwd(OutputDir)
  
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    write.csv(True_Labels_scID,paste('scID_',NumGenes,'_True_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Pred_Labels_scID,paste('scID_',NumGenes,'_Pred_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Total_Time_scID,paste('scID_',NumGenes,'_Total_Time.csv', sep = ''),row.names = FALSE)
  }
  else{
    write.csv(True_Labels_scID,'scID_True_Labels.csv',row.names = FALSE)
    write.csv(Pred_Labels_scID,'scID_Pred_Labels.csv',row.names = FALSE)
    write.csv(Total_Time_scID,'scID_Total_Time.csv',row.names = FALSE)
  }
}
