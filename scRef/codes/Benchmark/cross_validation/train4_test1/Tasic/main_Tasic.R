# import python package: sklearn.metrics
library(reticulate)
use_python('/local/zy/tools/anaconda3/bin/python3', required = T)
# py_config()
py_module_available('sklearn')

setwd('/local/zy/my_git/scRef/Benchmark/cross_validation/train4_test1/Tasic')
source('./Cross_Validation.R')
source('./method_functions.R')
source('./evaluate.R')

path.input <- '/mdshare/node9/zy/scRef/sc_data/'
path.output <- '/mdshare/node9/zy/scRef/cross_validation/train4_test1/Tasic/'

# setwd('/home/drizzle_zhang/my_git/scRef/Benchmark/cross_validation/train4_test1')
# source('./Cross_Validation.R')
# source('./method_functions.R')
# source('./evaluate.R')
# 
# path.input <- '/home/drizzle_zhang/scRef/'
# path.output <- '/home/drizzle_zhang/scRef/cross_validation/train4_test1/'

# filter data
dataset <- 'Tasic'
DataPath <- paste0(path.output, dataset, '.Rdata')
OutputDir <- path.output
# # function of data preparation
# prepare.data <- function(file.data.unlabeled, file.label.unlabeled,
#                          del.label = c('miss')) {
#     library(stringr)
#     data.unlabeled <- read.delim(file.data.unlabeled, row.names=1)
#     data.unlabeled <- floor(data.unlabeled)
#     names(data.unlabeled) <- str_replace_all(names(data.unlabeled), '_', '.')
#     names(data.unlabeled) <- str_replace_all(names(data.unlabeled), '-', '.')
#     # read label file
#     file.label.unlabeled <- file.label.unlabeled
#     label.unlabeled <- read.delim(file.label.unlabeled, row.names=1)
#     row.names(label.unlabeled) <- str_replace_all(row.names(label.unlabeled), '_', '.')
#     row.names(label.unlabeled) <- str_replace_all(row.names(label.unlabeled), '-', '.')
#     col.name1 <- names(data.unlabeled)[1]
#     if (substring(col.name1, 1, 1) == 'X') {
#         row.names(label.unlabeled) <- paste0('X', row.names(label.unlabeled))
#     }
#     # filter data
#     use.cols <- row.names(label.unlabeled)[!label.unlabeled[,1] %in% del.label]
#     data.filter <- data.unlabeled[,use.cols]
#     label.filter <- data.frame(label.unlabeled[use.cols,], row.names = use.cols)
# 
#     OUT <- list()
#     OUT$data.filter <- data.filter
#     OUT$label.filter <- label.filter
#     return(OUT)
# 
# }
# file.data.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat.txt')
# file.label.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat_cluster_original.txt')
# OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, del.label = c('Unclassified'))
# saveRDS(OUT, file = DataPath)
# 
# # generate cross validation dataset
# OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
# label.filter <- OUT$label.filter
# LabelsPath <- paste0(path.input, 'Habib_label.txt')
# write.table(label.filter, file = LabelsPath, sep = '\t', quote = F)
# Cross_Validation(LabelsPath, OutputDir)

CV_RDataPath <- paste0(path.output, 'CV_folds.RData')

# SingleR
run_SingleR(DataPath,LabelsPath,CV_RDataPath,OutputDir)
# scmap
run_scmap(DataPath,LabelsPath,CV_RDataPath,OutputDir)
# CHETAH
run_CHETAH(DataPath,LabelsPath,CV_RDataPath,OutputDir)
# scPred
run_scPred(DataPath,LabelsPath,CV_RDataPath,OutputDir)
# sciBet
run_sciBet(DataPath,LabelsPath,CV_RDataPath,OutputDir)
# scRef
run_scRef(DataPath,LabelsPath,CV_RDataPath,OutputDir)
# singleCellNet
run_singleCellNet(DataPath,LabelsPath,CV_RDataPath,OutputDir)
# scID
run_scID(DataPath,LabelsPath,CV_RDataPath,OutputDir)
# scClassify
run_scClassify(DataPath,LabelsPath,CV_RDataPath,OutputDir)


# heatmap
df.heatmap <- data.frame(stringsAsFactors = F)

# scRef
# run_scRef(DataPath,LabelsPath,CV_RDataPath,OutputDir)
TrueLabelsPath <- paste0(OutputDir, 'scRef_True_Labels.csv')
# PredLabelsPath <- paste0(OutputDir, 'scRef_Pred_Labels_cell.csv')
PredLabelsPath <- paste0(OutputDir, 'scRef_Pred_Labels.csv')
res.scRef <- evaluate(TrueLabelsPath, PredLabelsPath)
df.sub <- data.frame(term = names(res.scRef$F1), 
                     method = rep('scMAGIC', length(res.scRef$F1)),
                     value = res.scRef$F1, stringsAsFactors = F)
df.sub <- rbind(df.sub, 
                data.frame(term = 'macro F1', method = 'scMAGIC',
                           value = res.scRef$Mean_F1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Accuracy', method = 'scMAGIC',
                           value = res.scRef$Acc, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Weighted macro F1', method = 'scMAGIC',
                           value = res.scRef$WMean_F1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Balanced accuracy', method = 'scMAGIC',
                           value = res.scRef$balanced_acc, stringsAsFactors = F))
df.heatmap <- rbind(df.heatmap, df.sub)

# SingleR
# run_SingleR(DataPath,LabelsPath,CV_RDataPath,OutputDir)
TrueLabelsPath <- paste0(OutputDir, 'SingleR_True_Labels.csv')
PredLabelsPath <- paste0(OutputDir, 'SingleR_Pred_Labels.csv')
res.SingleR <- evaluate(TrueLabelsPath, PredLabelsPath)
df.sub <- data.frame(term = names(res.SingleR$F1), 
                     method = rep('SingleR', length(res.SingleR$F1)),
                     value = res.SingleR$F1, stringsAsFactors = F)
df.sub <- rbind(df.sub, 
                data.frame(term = 'macro F1', method = 'SingleR',
                           value = res.SingleR$Mean_F1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Accuracy', method = 'SingleR',
                           value = res.SingleR$Acc, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Weighted macro F1', method = 'SingleR',
                           value = res.SingleR$WMean_F1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Balanced accuracy', method = 'SingleR',
                           value = res.SingleR$balanced_acc, stringsAsFactors = F))
df.heatmap <- rbind(df.heatmap, df.sub)


# scmap
# run_scmap(DataPath,LabelsPath,CV_RDataPath,OutputDir)
TrueLabelsPath <- paste0(OutputDir, 'scmapcell_True_Labels.csv')
PredLabelsPath <- paste0(OutputDir, 'scmapcell_Pred_Labels.csv')
res.scmapcell <- evaluate(TrueLabelsPath, PredLabelsPath)

df.sub <- data.frame(term = names(res.scmapcell$F1), 
                     method = rep('scmap-cell', length(res.scmapcell$F1)),
                     value = res.scmapcell$F1, stringsAsFactors = F)
df.sub <- rbind(df.sub, 
                data.frame(term = 'macro F1', method = 'scmap-cell',
                           value = res.scmapcell$Mean_F1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Accuracy', method = 'scmap-cell',
                           value = res.scmapcell$Acc, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Weighted macro F1', method = 'scmap-cell',
                           value = res.scmapcell$WMean_F1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Balanced accuracy', method = 'scmap-cell',
                           value = res.scmapcell$balanced_acc, stringsAsFactors = F))
df.heatmap <- rbind(df.heatmap, df.sub)

PredLabelsPath <- paste0(OutputDir, 'scmapcluster_Pred_Labels.csv')
res.scmapcluster <- evaluate(TrueLabelsPath, PredLabelsPath)

df.sub <- data.frame(term = names(res.scmapcluster$F1), 
                     method = rep('scmap-cluster', length(res.scmapcluster$F1)),
                     value = res.scmapcluster$F1, stringsAsFactors = F)
df.sub <- rbind(df.sub, 
                data.frame(term = 'macro F1', method = 'scmap-cluster',
                           value = res.scmapcluster$Mean_F1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Accuracy', method = 'scmap-cluster',
                           value = res.scmapcluster$Acc, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Weighted macro F1', method = 'scmap-cluster',
                           value = res.scmapcluster$WMean_F1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Balanced accuracy', method = 'scmap-cluster',
                           value = res.scmapcluster$balanced_acc, stringsAsFactors = F))
df.heatmap <- rbind(df.heatmap, df.sub)

# CHETAH
# run_CHETAH(DataPath,LabelsPath,CV_RDataPath,OutputDir)
TrueLabelsPath <- paste0(OutputDir, 'CHETAH_True_Labels.csv')
PredLabelsPath <- paste0(OutputDir, 'CHETAH_Pred_Labels.csv')
res.CHETAH <- evaluate(TrueLabelsPath, PredLabelsPath)
df.sub <- data.frame(term = names(res.CHETAH$F1), 
                     method = rep('CHETAH', length(res.CHETAH$F1)),
                     value = res.CHETAH$F1, stringsAsFactors = F)
df.sub <- rbind(df.sub, 
                data.frame(term = 'macro F1', method = 'CHETAH',
                           value = res.CHETAH$Mean_F1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Accuracy', method = 'CHETAH',
                           value = res.CHETAH$Acc, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Weighted macro F1', method = 'CHETAH',
                           value = res.CHETAH$WMean_F1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Balanced accuracy', method = 'CHETAH',
                           value = res.CHETAH$balanced_acc, stringsAsFactors = F))
df.heatmap <- rbind(df.heatmap, df.sub)

# scPred
# run_scPred(DataPath,LabelsPath,CV_RDataPath,OutputDir)
TrueLabelsPath <- paste0(OutputDir, 'scPred_True_Labels.csv')
PredLabelsPath <- paste0(OutputDir, 'scPred_Pred_Labels.csv')
res.scPred <- evaluate(TrueLabelsPath, PredLabelsPath)
df.sub <- data.frame(term = names(res.scPred$F1), 
                     method = rep('scPred', length(res.scPred$F1)),
                     value = res.scPred$F1, stringsAsFactors = F)
df.sub <- rbind(df.sub, 
                data.frame(term = 'macro F1', method = 'scPred',
                           value = res.scPred$Mean_F1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Accuracy', method = 'scPred',
                           value = res.scPred$Acc, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Weighted macro F1', method = 'scPred',
                           value = res.scPred$WMean_F1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Balanced accuracy', method = 'scPred',
                           value = res.scPred$balanced_acc, stringsAsFactors = F))
df.heatmap <- rbind(df.heatmap, df.sub)

# sciBet
# run_sciBet(DataPath,LabelsPath,CV_RDataPath,OutputDir)
TrueLabelsPath <- paste0(OutputDir, 'sciBet_True_Labels.csv')
PredLabelsPath <- paste0(OutputDir, 'sciBet_Pred_Labels.csv')
res.sciBet <- evaluate(TrueLabelsPath, PredLabelsPath)
df.sub <- data.frame(term = names(res.sciBet$F1), 
                     method = rep('sciBet', length(res.sciBet$F1)),
                     value = res.sciBet$F1, stringsAsFactors = F)
df.sub <- rbind(df.sub, 
                data.frame(term = 'macro F1', method = 'sciBet',
                           value = res.sciBet$Mean_F1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Accuracy', method = 'sciBet',
                           value = res.sciBet$Acc, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Weighted macro F1', method = 'sciBet',
                           value = res.sciBet$WMean_F1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Balanced accuracy', method = 'sciBet',
                           value = res.sciBet$balanced_acc, stringsAsFactors = F))
df.heatmap <- rbind(df.heatmap, df.sub)

# singleCellNet
# run_singleCellNet(DataPath,LabelsPath,CV_RDataPath,OutputDir)
TrueLabelsPath <- paste0(OutputDir, 'singleCellNet_True_Labels.csv')
PredLabelsPath <- paste0(OutputDir, 'singleCellNet_Pred_Labels.csv')
res.singleCellNet <- evaluate(TrueLabelsPath, PredLabelsPath)
df.sub <- data.frame(term = names(res.singleCellNet$F1), 
                     method = rep('singleCellNet', length(res.singleCellNet$F1)),
                     value = res.singleCellNet$F1, stringsAsFactors = F)
df.sub <- rbind(df.sub, 
                data.frame(term = 'macro F1', method = 'singleCellNet',
                           value = res.singleCellNet$Mean_F1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Accuracy', method = 'singleCellNet',
                           value = res.singleCellNet$Acc, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Weighted macro F1', method = 'singleCellNet',
                           value = res.singleCellNet$WMean_F1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Balanced accuracy', method = 'singleCellNet',
                           value = res.singleCellNet$balanced_acc, stringsAsFactors = F))
df.heatmap <- rbind(df.heatmap, df.sub)

# CaSTLe
# run_CaSTLe(DataPath,LabelsPath,CV_RDataPath,OutputDir)
# TrueLabelsPath <- paste0(OutputDir, 'CaSTLe_True_Labels.csv')
# PredLabelsPath <- paste0(OutputDir, 'CaSTLe_Pred_Labels.csv')
# res.CaSTLe <- evaluate(TrueLabelsPath, PredLabelsPath)
# df.sub <- data.frame(term = names(res.CaSTLe$F1), 
#                      method = rep('CaSTLe', length(res.CaSTLe$F1)),
#                      value = res.CaSTLe$F1, stringsAsFactors = F)
# df.sub <- rbind(df.sub, 
#                 data.frame(term = 'macro F1', method = 'CaSTLe',
#                            value = res.CaSTLe$Mean_F1, stringsAsFactors = F))
# df.sub <- rbind(df.sub, 
#                 data.frame(term = 'Accuracy', method = 'CaSTLe',
#                            value = res.CaSTLe$Acc, stringsAsFactors = F))
# df.heatmap <- rbind(df.heatmap, df.sub)

# scID
# run_scID(DataPath,LabelsPath,CV_RDataPath,OutputDir)
TrueLabelsPath <- paste0(OutputDir, 'scID_True_Labels.csv')
PredLabelsPath <- paste0(OutputDir, 'scID_Pred_Labels.csv')
res.scID <- evaluate(TrueLabelsPath, PredLabelsPath)
df.sub <- data.frame(term = names(res.scID$F1), 
                     method = rep('scID', length(res.scID$F1)),
                     value = res.scID$F1, stringsAsFactors = F)
df.sub <- rbind(df.sub, 
                data.frame(term = 'macro F1', method = 'scID',
                           value = res.scID$Mean_F1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Accuracy', method = 'scID',
                           value = res.scID$Acc, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Weighted macro F1', method = 'scID',
                           value = res.scID$WMean_F1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Balanced accuracy', method = 'scID',
                           value = res.scID$balanced_acc, stringsAsFactors = F))
df.heatmap <- rbind(df.heatmap, df.sub)

# scClassify
# run_scClassify(DataPath,LabelsPath,CV_RDataPath,OutputDir)
TrueLabelsPath <- paste0(OutputDir, 'scClassify_True_Labels.csv')
PredLabelsPath <- paste0(OutputDir, 'scClassify_Pred_Labels.csv')
res.scClassify <- evaluate(TrueLabelsPath, PredLabelsPath)
df.sub <- data.frame(term = names(res.scClassify$F1), 
                     method = rep('scClassify', length(res.scClassify$F1)),
                     value = res.scClassify$F1, stringsAsFactors = F)
df.sub <- rbind(df.sub, 
                data.frame(term = 'macro F1', method = 'scClassify',
                           value = res.scClassify$Mean_F1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Accuracy', method = 'scClassify',
                           value = res.scClassify$Acc, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Weighted macro F1', method = 'scClassify',
                           value = res.scClassify$WMean_F1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Balanced accuracy', method = 'scClassify',
                           value = res.scClassify$balanced_acc, stringsAsFactors = F))
df.heatmap <- rbind(df.heatmap, df.sub)

unique.term <- unique(df.heatmap$term)
df.heatmap$term <- factor(df.heatmap$term, levels = unique.term)
df.acc <- df.heatmap[df.heatmap$term == 'Accuracy', ]
df.heatmap$method <- factor(df.heatmap$method, 
                            levels = df.acc$method[order(df.acc$value, decreasing = T)])

# save results
file.res <- paste0(path.output, 'results_', dataset, '.txt')
write.table(df.heatmap, file = file.res, sep = '\t', quote = F, row.names = F)


# plot heatmap
library(ggplot2)
plot.heatmap <- ggplot(data = df.heatmap, aes(method, term)) + 
    geom_tile(aes(fill = value)) + 
    scale_fill_continuous(low = "#FFFAFA", high = "#A52A2A") + 
    theme_bw() +
    theme(
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.6)
    ) + 
    geom_text(aes(label = round(value, 2)), family = "Arial", size = 2.5)
ggsave(filename = 'heatmap.png', path = path.output, plot = plot.heatmap,
       units = 'cm', height = 18, width = 18)

# plot barplot
df.barplot <- df.heatmap[as.character(df.heatmap$term) %in% c('Accuracy', 'macro F1', 'Weighted macro F1'),]
plot.bar <- ggplot(df.barplot,
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
ggsave(filename = 'barplot.png', 
       path = path.output, plot = plot.bar,
       units = 'cm', height = 12, width = 18)




