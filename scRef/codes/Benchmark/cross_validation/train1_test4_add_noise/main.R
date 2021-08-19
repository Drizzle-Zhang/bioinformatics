# setwd('/home/zy/my_git/scRef/Benchmark/cross_validation/train1_test4')
# source('./Cross_Validation.R')
# source('./method_functions.R')
# source('./evaluate.R')
# 
# path.input <- '/home/zy/scRef/'
# path.output <- '/home/zy/scRef/cross_validation/train1_test4/'

setwd('/home/drizzle_zhang/my_git/scRef/Benchmark/cross_validation/train1_test4_add_noise')
source('./Cross_Validation.R')
source('./method_functions.R')
source('./evaluate.R')

path.input <- '/home/drizzle_zhang/scRef/'
path.output <- '/home/drizzle_zhang/scRef/cross_validation/train1_test4_add_noise/'

# generate cross validation dataset
LabelsPath <- paste0(path.input, 'summary/Zeisel_exp_sc_mat_cluster_original.txt')
OutputDir <- path.output
if (!file.exists(OutputDir)) {
    dir.create(OutputDir)
}
# Cross_Validation(LabelsPath, OutputDir)

DataPath.origin <- paste0(path.input, 'summary/Zeisel_exp_sc_mat.txt')
############
# add noise
# Data <- read.delim(DataPath.origin,row.names = 1)
# set.seed(123)
# addNOI=function(x){
#     M=mean(x)
#     y=x+M/3*(runif(length(x))*2-1)
#     return(y)
# }
# nData=t(apply(Data,1,addNOI))
# nData[which(nData<0)]=0
# rownames(nData)=rownames(Data)
# colnames(nData)=colnames(Data)
# nData <- round(nData, digits = 4)
# file.noise <- paste0(OutputDir, 'Data_noise.txt')
# write.table(nData, file = file.noise, sep = '\t')
############

DataPath <- paste0(OutputDir, 'Data_noise.txt')
CV_RDataPath <- paste0(path.output, 'CV_folds.RData')

# run methods
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
# run_singleCellNet(DataPath,LabelsPath,CV_RDataPath,OutputDir)
# CaSTLe
run_CaSTLe(DataPath,LabelsPath,CV_RDataPath,OutputDir)
# scID
run_scID(DataPath,LabelsPath,CV_RDataPath,OutputDir)


# heatmap
df.heatmap <- data.frame(stringsAsFactors = F)

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
df.heatmap <- rbind(df.heatmap, df.sub)

# scRef
# run_scRef(DataPath,LabelsPath,CV_RDataPath,OutputDir)
TrueLabelsPath <- paste0(OutputDir, 'scRef_True_Labels.csv')
# PredLabelsPath <- paste0(OutputDir, 'scRef_Pred_Labels_cell.csv')
PredLabelsPath <- paste0(OutputDir, 'scRef_Pred_Labels.csv')
res.scRef <- evaluate(TrueLabelsPath, PredLabelsPath)
df.sub <- data.frame(term = names(res.scRef$F1), 
                     method = rep('scRef', length(res.scRef$F1)),
                     value = res.scRef$F1, stringsAsFactors = F)
df.sub <- rbind(df.sub, 
                data.frame(term = 'macro F1', method = 'scRef',
                           value = res.scRef$Mean_F1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Accuracy', method = 'scRef',
                           value = res.scRef$Acc, stringsAsFactors = F))
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
df.heatmap <- rbind(df.heatmap, df.sub)

# CaSTLe
# run_CaSTLe(DataPath,LabelsPath,CV_RDataPath,OutputDir)
TrueLabelsPath <- paste0(OutputDir, 'CaSTLe_True_Labels.csv')
PredLabelsPath <- paste0(OutputDir, 'CaSTLe_Pred_Labels.csv')
res.CaSTLe <- evaluate(TrueLabelsPath, PredLabelsPath)
df.sub <- data.frame(term = names(res.CaSTLe$F1), 
                     method = rep('CaSTLe', length(res.CaSTLe$F1)),
                     value = res.CaSTLe$F1, stringsAsFactors = F)
df.sub <- rbind(df.sub, 
                data.frame(term = 'macro F1', method = 'CaSTLe',
                           value = res.CaSTLe$Mean_F1, stringsAsFactors = F))
df.sub <- rbind(df.sub, 
                data.frame(term = 'Accuracy', method = 'CaSTLe',
                           value = res.CaSTLe$Acc, stringsAsFactors = F))
df.heatmap <- rbind(df.heatmap, df.sub)

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
df.heatmap <- rbind(df.heatmap, df.sub)
unique.term <- unique(df.heatmap$term)
df.heatmap$term <- factor(df.heatmap$term, levels = unique.term)
df.acc <- df.heatmap[df.heatmap$term == 'macro F1', ]
df.heatmap$method <- factor(df.heatmap$method, 
                            levels = df.acc$method[order(df.acc$value, decreasing = T)])

# plot heatmap
library(ggplot2)
plot.heatmap <- ggplot(data = df.heatmap, aes(method, term)) + 
    geom_tile(aes(fill = value)) + 
    scale_fill_continuous(low = "#FFFAFA", high = "#A52A2A") + 
    labs(fill = '') + 
    theme_bw() +
    theme(
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.6)
    ) + 
    geom_text(aes(label = round(value, 3)), family = "Arial", size = 2.5)
path <- '/home/drizzle_zhang/scRef/cross_validation/train1_test4_add_noise'
ggsave(filename = 'heatmap.png', path = path, plot = plot.heatmap,
       units = 'cm', height = 10, width = 18)

