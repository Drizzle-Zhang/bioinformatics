path.scHCL <- '/mdshare/node9/zy/scRef/scHCL/'

library(scMAGIC)
data("HCL_ref")
data("MCA_ref")

library(scHCL)
gettissue <- function(x,Num=3){
    top_markers <-order(x,decreasing = T)[1:Num]
    return(top_markers)
}

scHCL <- function(scdata,ref.expr,numbers_plot=3){
    library(tidyverse)
    tst.expr <- data.frame(matrix(nrow = dim(ref.expr)[1],ncol=dim(scdata)[2]))
    rownames(tst.expr)<-rownames(ref.expr)
    colnames(tst.expr)<-colnames(scdata)
    for (i in rownames(tst.expr)) {
        if(i%in%rownames(scdata))tst.expr[i,]<- scdata[i,]
    }
    tst.expr[is.na(tst.expr)]<-0
    tst.expr<-as.matrix(t(t(tst.expr)/colSums(tst.expr))*100000)
    tst.expr<-log(tst.expr+1)
    cors <- cor(log(ref.expr+1),tst.expr)
    
    cors_index <- apply(cors,2,gettissue,numbers_plot)
    cors_index <- sort(unique(as.integer(cors_index)))
    scblast.result <- apply(cors,2,function(x) rownames(cors)[which.max(x)])
    cors_in = cors[cors_index,]
    colnames(cors_in)<-colnames(scdata)
    cors_out = reshape2::melt(cors_in)[c(2,1,3)]
    colnames(cors_out)<- c("Cell","Cell type","Score")
    cors_out <- as.data.frame(cors_out %>% group_by(Cell) %>% top_n(n=numbers_plot,wt=Score))
    
    
    result <- list()
    cors[which(is.na(cors))]<-0
    result[["cors_matrix"]] <- cors
    result[['top_cors']]<-numbers_plot
    result[['scHCL']]<-scblast.result
    result[['scHCL_probility']]<-cors_out
    return(result)
}

#################### pbmcsca
library(Seurat)
library(SeuratData)
data("pbmcsca")
dataset <- 'hPBMC_10Xv2'
exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method %in%
                                                      c('10x Chromium (v2)')])
label_sc <- data.frame(
    annotation = as.character(pbmcsca$CellType)[pbmcsca$Method %in%
                                                    c('10x Chromium (v2)')],
    row.names = colnames(exp_sc_mat))

list.overlap <- get_overlap_genes(exp_sc_mat, HCL_ref)

exp_mat <- list.overlap$exp_sc_mat
time1 <- Sys.time()
hcl_result_pbmc <- scHCL(scdata = exp_mat, ref.expr = HCL_ref, numbers_plot = 3)
time2 <- Sys.time()
time.diff <- difftime(time2, time1, units = 'mins')
hcl_result_pbmc$run_time <- time.diff
table(label_sc$annotation, hcl_result_pbmc$scHCL)
file.pbmc <- paste0(path.scHCL, 'scHCL_', dataset, '.Rdata')
saveRDS(hcl_result_pbmc, file.pbmc)
# file.mat.pbmc <- paste0(path.scHCL, 'mat_', dataset, '.txt')
# write.table(exp_sc_mat[sel.genes,], file = file.mat.pbmc, sep = '\t', quote = F)

# assessment
file.pbmc <- paste0(path.scHCL, 'scHCL_', dataset, '.Rdata')
hcl_result_pbmc <- readRDS(file.pbmc)
true.tags <- label_sc$annotation
true.tags[true.tags %in% c('CD14+ monocyte', 'CD16+ monocyte')] <- 'Monocyte'
true.tags[true.tags %in% c('CD4+ T cell', 'Cytotoxic T cell')] <- 'T cell'
true.tags[true.tags %in% c('Plasmacytoid dendritic cell')] <- 'Dendritic cell'

df.simple <- data.frame(check.names = F)
for (col in colnames(HCL_ref)) {
    col.split <- strsplit(col, split = '_')[[1]]
    col.simple <- paste(col.split[1:(length(col.split)-1)], collapse = '-')
    df.simple <- rbind(df.simple, data.frame(col = col,
                                             col.simple = col.simple))
}
row.names(df.simple) <- df.simple$col
pred.scHCL <- hcl_result_pbmc$scHCL
for (cell in unique(pred.scHCL)) {
    pred.scHCL[pred.scHCL == cell] <- df.simple[cell, 'col.simple']
}
pred.scHCL[pred.scHCL %in% c('B cell (Centrocyte)', 'B cell(Centrocyte)',
                             'B cell(Plasmocyte)', 'Proliferating  B cell')] <- 'B cell'
pred.scHCL[pred.scHCL %in% c('Conventional dendritic cell', 'Dendritic cell-FCER1A high',
                             'Plasmacytoid dendritic cell')] <- 'Dendritic cell'
pred.scHCL[pred.scHCL %in% c('CD14+ monocyte', 'CD16+ monocyte')] <- 'Megakaryocyte'
pred.scHCL[pred.scHCL %in% c('Monocyte-FCGR3A high', 'Monocyte-IGHG4 high',
                             'Monocyte-TPPP3 high')] <- 'Monocyte'
pred.scHCL[pred.scHCL %in% c('NK cell', 'CD16+ monocyte')] <- 'Natural killer cell'
pred.scHCL[pred.scHCL %in% c('activative T cell', 'Actived T cell', 'Proliferating T cell',
                             'T cell-CCL5 high', 'T cell-GNLY high', 'T cell-IL7R high',
                             'T cell-TRAC high')] <- 'T cell'
pred.scHCL[pred.scHCL %in% c('Arterial endothelial cell-GJA5 high', 'Multipotential  progenitor cell',
                             'Proliferating cell')] <- 'Impurity'
pred.scHCL[pred.scHCL %in% c('Macrophage-CXCL2 high', 'Macrophage-FCGR3A high',
                             'Macrophage-HLA-DRA high', 'Motile liver macrophage')] <- 'Macrophage'
pred.scHCL[pred.scHCL %in% c('Neutrophil-ELANE high', 'Neutrophil-MMP high',
                             'Neutrophil-MMP9 high')] <- 'Neutrophil'
table(true.tags, pred.scHCL)

library(reticulate)
use_python('/local/zy/tools/anaconda3/bin/python3', required = T)
# py_config()
py_module_available('sklearn')
metrics <- import('sklearn.metrics')

accuracy <- metrics$accuracy_score(true.tags, pred.scHCL)
balanced_acc <- metrics$balanced_accuracy_score(true.tags, pred.scHCL)
balanced.accuracy.rm.unassigned <- metrics$balanced_accuracy_score(pred.scHCL, true.tags)
res.scHCL <- list()
res.scHCL$accuracy <- accuracy
res.scHCL$accuracy.rm.unassigned <- accuracy
res.scHCL$balanced.accuracy <- balanced_acc
res.scHCL$balanced.accuracy.rm.unassigned <- balanced.accuracy.rm.unassigned
path.res <- '/mdshare/node9/zy/scRef/figure/atlas_anno/'
file.res.scHCL <- paste0(path.res, 'RES_HCL_', dataset, '_scHCL.Rdata')
saveRDS(res.scHCL, file.res.scHCL)


######################### panc8
data("panc8")
dataset <- 'panc8_indrop'
exp_sc_mat <- as.matrix(panc8@assays$RNA@counts[, panc8$dataset %in%
                                                    c('indrop1', 'indrop2', 'indrop3', 'indrop4')])
label_sc <- data.frame(
    annotations = as.character(panc8$celltype)[panc8$dataset %in%
                                                   c('indrop1', 'indrop2', 'indrop3', 'indrop4')],
    row.names = colnames(exp_sc_mat))

list.overlap <- get_overlap_genes(exp_sc_mat, HCL_ref)

exp_mat <- list.overlap$exp_sc_mat
time1 <- Sys.time()
hcl_result_panc <- scHCL(scdata = exp_mat, ref.expr = HCL_ref, numbers_plot = 3)
time2 <- Sys.time()
time.diff <- difftime(time2, time1, units = 'mins')
hcl_result_panc$run_time <- time.diff
table(label_sc$annotation, hcl_result_panc$scHCL)
file.panc <- paste0(path.scHCL, 'scHCL_', dataset, '.Rdata')
saveRDS(hcl_result_panc, file.panc)

# file.mat.panc <- paste0(path.scHCL, 'mat_', dataset, '.txt')
# write.table(exp_sc_mat[sel.genes,], file = file.mat.panc, sep = '\t', quote = F)

# assessment
file.panc <- paste0(path.scHCL, 'scHCL_', dataset, '.Rdata')
hcl_result_panc <- readRDS(file.panc)
true.tags <- label_sc$annotation
true.tags[true.tags %in% c('acinar')] <- 'Acinar cell'
true.tags[true.tags %in% c('activated_stellate')] <- 'ASC'
true.tags[true.tags %in% c('alpha')] <- 'Alpha cell'
true.tags[true.tags %in% c('beta')] <- 'Beta cell'
true.tags[true.tags %in% c('delta')] <- 'Delta cell'
true.tags[true.tags %in% c('ductal')] <- 'Dutal cell'
true.tags[true.tags %in% c('endothelial')] <- 'Endothelial Cell'
true.tags[true.tags %in% c('epsilon')] <- 'Epsilon cell'
true.tags[true.tags %in% c('gamma')] <- 'Gamma cell'
true.tags[true.tags %in% c('macrophage')] <- 'Macrophage'
true.tags[true.tags %in% c('mast')] <- 'Mast cell'
true.tags[true.tags %in% c('quiescent_stellate')] <- 'QSC'
true.tags[true.tags %in% c('schwann')] <- 'Schwann cell'
df.simple <- data.frame(check.names = F)
for (col in colnames(HCL_ref)) {
    col.split <- strsplit(col, split = '_')[[1]]
    col.simple <- paste(col.split[1:(length(col.split)-1)], collapse = '-')
    df.simple <- rbind(df.simple, data.frame(col = col,
                                             col.simple = col.simple))
}
row.names(df.simple) <- df.simple$col
pred.scHCL <- hcl_result_pbmc$scHCL
for (cell in unique(pred.scHCL)) {
    pred.scHCL[pred.scHCL == cell] <- df.simple[cell, 'col.simple']
}
pred.scHCL[pred.scHCL %in% c('Acinar cell-CPA1 high', 'Acinar cell-REG1B high',
                             'Acniar cell-ANXA4 high')] <- 'Acinar cell'
pred.scHCL[pred.scHCL %in% c('Adipocyte-SPP1 high')] <- 'Adipocyte'
pred.scHCL[pred.scHCL %in% c('Arterial endothelial cell', 'Endothelial cell in EMT',
                             'Endothelial cell-ACKR1 high', 'Endothelial cell-APC',
                             'Endothelial cell-CCL21 high', 'Endothelial cell-ESM1 high',
                             'Endothelial cell-FABP4 high', 'Endothelial cell-IGFBP3 high',
                             'Endothelial cell-IGFBP5 high', 'Endothelial cell-SELE high',
                             'Endothelial cell-SOCS3 high', 'Endothelial cell-SPARCL1 high',
                             'Endothelial cell-TMEM100 high', 'Endothelial cell-VWF high',
                             'Vascular endothelial cell-IGFBP3 high')] <- 'Endothelial cell'
pred.scHCL[pred.scHCL %in% c('Basal cell-KRT6A high')] <- 'Basal  cell'
pred.scHCL[pred.scHCL %in% c('Chromaffin cell-VIP high')] <- 'Chromaffin cell'
pred.scHCL[pred.scHCL %in% c('Actived T cell', 'CD8 T cell',
                             'Proliferating T cell', 'T cell-GNLY high', 'T cell-IL7R high',
                             'T cell-TRAC high')] <- 'T cell'
pred.scHCL[pred.scHCL %in% c('Epithelial cell-KRT17 high', 'Luminal epithelium ', 
                             'Mucous Epithelial cell-REG1A high', 
                             'Mucous Epithelial cell-TFF1 high', 
                             'Secretory epithelial cell',
                             'Unknown Epithelial cell-FOS high',
                             'Enterocyte-SLC26A3 high')] <- 'Epithelial cell'
pred.scHCL[pred.scHCL %in% c('Exocrine cell-SAA1 high')] <- 'Exocrine cell'
pred.scHCL[pred.scHCL %in% c('Fibroblast-APOD high', 'Fibroblast-PTX3 high', 
                             'Myofibroblast-POSTN high')] <- 'Fibroblast'
pred.scHCL[pred.scHCL %in% c('M1 Macrophage', 'Macrophage-C1QB high',
                             'Macrophage-RNASE1 high')] <- 'Macrophage'
pred.scHCL[pred.scHCL %in% c('Neutrophil-IL1B high')] <- 'Neutrophil'
pred.scHCL[pred.scHCL %in% c('Neuroendocrine cell-SST high')] <- 'Neuroendocrine cell'
pred.scHCL[pred.scHCL %in% c('Pit cell-FOXQ1 high')] <- 'Pit cell'
pred.scHCL[pred.scHCL %in% c('Smooth muscel cell', 'Smooth muscle cell-CCL19 high',
                             'Smooth muscle cell-MYL9 high', 
                             'Smooth muscle cell-PDK4 high')] <- 'Smooth muscle cell'
pred.scHCL[pred.scHCL %in% c('Stromal cell-ERRFI1 high', 'Stromal cell-LUM high',
                             'Stromal cell-PLA2G2A high', 
                             'Smooth muscle cell-PDK4 high')] <- 'Stromal cell'
mytable <- table(true.tags, pred.scHCL)


list.query.cell <- list()
list.query.cell[['Acinar cell']] <- c('Acinar cell', 'Exocrine cell')
list.query.cell[['Alpha cell']] <- c('Alpha cell', 'Endocrine cell')
list.query.cell[['ASC']] <- c('Stromal cell')
list.query.cell[['Beta cell']] <- c('Beta cell', 'Endocrine cell')
list.query.cell[['Delta cell']] <- c('Endocrine cell')
list.query.cell[['Ductal cell']] <- c('Ductal cell')
list.query.cell[['Endothelial Cell']] <- c('Endothelial Cell')
list.query.cell[['Epsilon cell']] <- c('Endocrine cell')
list.query.cell[['Gamma cell']] <- c('Endocrine cell')
list.query.cell[['Macrophage']] <- c('Macrophage')
list.query.cell[['Mast cell']] <- c('Mast cell')
list.query.cell[['QSC']] <- c('Stromal cell')
list.query.cell[['Schwann cell']] <- c('Schwann cell')

total.cell.num <- length(true.tags)
sum.correct <- 0
vec.acc <- c()
for (cell in rownames(mytable)) {
    sub.pred <- pred.scHCL[true.tags == cell]
    sub.num <- length(sub.pred)
    sub.correct.num <- sum(sub.pred %in% list.query.cell[[cell]])
    vec.acc <- c(vec.acc, sub.correct.num/sub.num)
    sum.correct <- sum.correct + sub.correct.num
}
accuracy <- sum.correct/total.cell.num
balanced_acc <- mean(vec.acc)

list.ref.cell <- list()
list.ref.cell[['Acinar cell']] <- c('Acinar cell')
list.ref.cell[['Alpha cell']] <- c('Alpha cell')
list.ref.cell[['Beta cell']] <- c('Beta cell')
list.ref.cell[['Ductal cell']] <- c('Ductal cell')
list.ref.cell[['Endocrine cell']] <- c('Alpha cell', 'Beta cell', 'Delta cell', 
                                       'Epsilon cell','Gamma cell')
list.ref.cell[['Endothelial Cell']] <- c('Endothelial Cell')
list.ref.cell[['Exocrine cell']] <- c('Acinar cell')
list.ref.cell[['Macrophage']] <- c('Macrophage')
list.ref.cell[['Mast cell']] <- c('Mast cell')
list.ref.cell[['Stromal cell']] <- c('ASC', 'QSC')

ref.cells <- names(table(pred.scHCL))[table(pred.scHCL) > 15]
vec.acc <- c()
for (cell in ref.cells) {
    sub.true <- true.tags[pred.scHCL == cell]
    sub.num <- length(sub.true)
    sub.correct.num <- sum(sub.true %in% list.ref.cell[[cell]])
    vec.acc <- c(vec.acc, sub.correct.num/sub.num)
}
balanced.accuracy.rm.unassigned <- mean(vec.acc)

res.scHCL <- list()
res.scHCL$accuracy <- accuracy
res.scHCL$accuracy.rm.unassigned <- accuracy
res.scHCL$balanced.accuracy <- balanced_acc
res.scHCL$balanced.accuracy.rm.unassigned <- balanced.accuracy.rm.unassigned
path.res <- '/mdshare/node9/zy/scRef/figure/atlas_anno/'
file.res.scHCL <- paste0(path.res, 'RES_HCL_', dataset, '_scHCL.Rdata')
saveRDS(res.scHCL, file.res.scHCL)


###################### Neocortex
list.target <- readRDS('/local/zy/my_git/scMAGIC_scripts/data/MouseNeocortex.Rdata')
dataset <- 'Tasic2018'
exp_sc_mat <- list.target$mat_exp
label_sc <-list.target$label

list.overlap <- get_overlap_genes(exp_sc_mat, MCA_ref)

exp_mat <- list.overlap$exp_sc_mat
time1 <- Sys.time()
hcl_result_brain <- scHCL(scdata = exp_mat, ref.expr = MCA_ref, numbers_plot = 3)
time2 <- Sys.time()
time.diff <- difftime(time2, time1, units = 'mins')
hcl_result_brain$run_time <- time.diff
table(label_sc, hcl_result_brain$scHCL)
file.brain <- paste0(path.scHCL, 'scHCL_', dataset, '.Rdata')
saveRDS(hcl_result_brain, file.brain)

# assessment
file.brain <- paste0(path.scHCL, 'scHCL_', dataset, '.Rdata')
hcl_result_brain <- readRDS(file.brain)
true.tags <- label_sc
true.tags[true.tags %in% c('GABAergic Neuron', 'Glutamatergic Neuron')] <- 'Neuron'
true.tags[true.tags %in% c('VLMC', 'Pericyte')] <- 'Unassigned'

df.simple <- data.frame(check.names = F)
for (col in colnames(MCA_ref)) {
    col.split <- strsplit(col, split = "(", fixed = T)[[1]]
    col.simple <- col.split[1]
    # col.simple <- paste(col.split[1:(length(col.split)-1)], collapse = '-')
    df.simple <- rbind(df.simple, data.frame(col = col,
                                             col.simple = col.simple))
}
row.names(df.simple) <- df.simple$col
pred.scHCL <- hcl_result_brain$scHCL
for (cell in unique(pred.scHCL)) {
    pred.scHCL[pred.scHCL == cell] <- df.simple[cell, 'col.simple']
}
pred.scHCL[pred.scHCL %in% c('Astrocyte_Atp1b2 high')] <- 'Astrocyte'
pred.scHCL[pred.scHCL %in% c('Atrial cardiomyocyte_Acta2 high', 
                             'Ventricle cardiomyocyte_Kcnj8 high')] <- 'Cardiomyocyte'
pred.scHCL[pred.scHCL %in% c('Dopaminergic neurons', 'Granule neurons',
                             'Hippocampus neurons_Asic4 high', 'Neuron', 'Purkinje cell')] <- 'Neuron'
pred.scHCL[pred.scHCL %in% c('Pyramidal neuron cell', 'Schwann cell')] <- 'Other Neuron'
pred.scHCL[pred.scHCL %in% c('Endothelial cell_Fabp4 high', 'Endothelial cell_Ly6c1 high',
                             'Endothelial cell_Tm4sf1 high', 'Endothelial cells_Vwf high')] <- 'Endothelial cell'
pred.scHCL[pred.scHCL %in% c('Muscle cell_Lrrc15 high', 'Muscle cell_Mgp high',
                             'Muscle cell_Myl9 high')] <- 'Muscle cell'
pred.scHCL[pred.scHCL %in% c('Mesenchymal stromal cell', 'Stromal cell_Car3 high', 
                             'Stromal cell_Col3a1 high', 'Stromal cell_Cxcl14 high',
                             'Stromal cell_Dcn high', 'Stromal cell_Fmod high', 'Stromal cell_Gas6 high',
                             'Stromal cell_Mfap4 high', 'Stromal cell_Ptgds high',
                             'Stromal cell_Smoc2 high')] <- 'Stromal cell'
pred.scHCL[pred.scHCL %in% c('Oligodendrocyte precursor cell')] <- 'Oligodendrocyte & OPC'
pred.scHCL[pred.scHCL %in% c('Macrophage', 'Macrophage_C1qc high',
                             'Macrophage_Pf4 high', 'Microglia')] <- 'PVM & Microglia'
pred.scHCL[pred.scHCL %in% c('Smooth muscle cell_Mgp high', 
                             'Smooth muscle cell_Rgs5 high')] <- 'Smooth muscle cell'
pred.scHCL[pred.scHCL %in% c('Erythroblast_Hbb-y high', 'MSC',
                             'T cell_Pclaf high')] <- 'Other cells'
table(true.tags, pred.scHCL)

library(reticulate)
use_python('/local/zy/tools/anaconda3/bin/python3', required = T)
# py_config()
py_module_available('sklearn')
metrics <- import('sklearn.metrics')

accuracy <- metrics$accuracy_score(true.tags, pred.scHCL)
balanced_acc <- metrics$balanced_accuracy_score(true.tags, pred.scHCL)
balanced.accuracy.rm.unassigned <- metrics$balanced_accuracy_score(pred.scHCL, true.tags)
res.scHCL <- list()
res.scHCL$accuracy <- accuracy
res.scHCL$accuracy.rm.unassigned <- accuracy
res.scHCL$balanced.accuracy <- balanced_acc
res.scHCL$balanced.accuracy.rm.unassigned <- balanced.accuracy.rm.unassigned
path.res <- '/mdshare/node9/zy/scRef/figure/atlas_anno/'
file.res.scHCL <- paste0(path.res, 'RES_MCA_', dataset, '_scHCL.Rdata')
saveRDS(res.scHCL, file.res.scHCL)


###################### Duodenum
path.input <- '/mdshare/node9/zy/scRef/sc_data/MouseSmallIntestinalEpithelium_SingleCell_Haber2018'
dataset <- 'Haber_Duodenum'
file.data.unlabeled <- paste0(path.input, '/raw_count/cell_exp.txt')
file.label.unlabeled <- paste0(path.input, '/raw_count/cell_meta.txt')
data.unlabeled <- read.delim(file.data.unlabeled, row.names=1)
label.unlabeled <- read.delim(file.label.unlabeled, row.names=1)
label.unlabeled$location <- unlist(lapply(strsplit(row.names(label.unlabeled), '_'),
                                          function(x) {x[2]}))
label.unlabeled$location[!(label.unlabeled$location %in%
                               c('Ileum', 'Duodenum', 'Jejunum'))] <- 'Other'
label.Duodenum <- label.unlabeled[label.unlabeled$location == 'Duodenum',]
mat.Duodenum <- data.unlabeled[, rownames(label.Duodenum)]
mat.Duodenum <- mat.Duodenum[!(is.na(mat.Duodenum[,1])),]
exp_sc_mat <- mat.Duodenum
label_sc <- label.Duodenum

list.overlap <- get_overlap_genes(exp_sc_mat, MCA_ref)

exp_mat <- list.overlap$exp_sc_mat
time1 <- Sys.time()
hcl_result_Duodenum <- scHCL(scdata = exp_mat, ref.expr = MCA_ref, numbers_plot = 3)
time2 <- Sys.time()
time.diff <- difftime(time2, time1, units = 'mins')
hcl_result_Duodenum$run_time <- time.diff
table(label_sc$Cluster, hcl_result_Duodenum$scHCL)
file.Duodenum <- paste0(path.scHCL, 'scHCL_', dataset, '.Rdata')
saveRDS(hcl_result_Duodenum, file.Duodenum)

# assessment
file.Duodenum <- paste0(path.scHCL, 'scHCL_', dataset, '.Rdata')
hcl_result_Duodenum <- readRDS(file.Duodenum)
true.tags <- label_sc$Cluster
true.tags[true.tags %in% c('Goblet')] <- 'Goblet cell'
true.tags[true.tags %in% c('Enteroendocrine')] <- 'Enteroendocrine cell'
true.tags[true.tags %in% c('Paneth')] <- 'Paneth cell'
true.tags[true.tags %in% c('Tuft')] <- 'Tuft cell'
true.tags[true.tags %in% c('EP', 'TA', 'Stem')] <- 'Unassigned'

df.simple <- data.frame(check.names = F)
for (col in colnames(MCA_ref)) {
    col.split <- strsplit(col, split = "(", fixed = T)[[1]]
    col.simple <- col.split[1]
    # col.simple <- paste(col.split[1:(length(col.split)-1)], collapse = '-')
    df.simple <- rbind(df.simple, data.frame(col = col,
                                             col.simple = col.simple))
}
row.names(df.simple) <- df.simple$col
pred.scHCL <- hcl_result_Duodenum$scHCL
for (cell in unique(pred.scHCL)) {
    pred.scHCL[pred.scHCL == cell] <- df.simple[cell, 'col.simple']
}
pred.scHCL[pred.scHCL %in% c('Enteroendocrine', 'S cell_Chgb high', 
                             'S cell_Gip high')] <- 'Enteroendocrine cell'
pred.scHCL[pred.scHCL %in% c('Epithelial cell_Bex1 high', 'Epithelial cell_Lgals2 high',
                             'Epithelial cell_Gkn3 high')] <- 'Epithelial cell'
pred.scHCL[pred.scHCL %in% c('Epithelium of small intestinal villi_Fabp1 high', 
                             'Epithelium of small intestinal villi_mt-Nd1 high',
                             'Epithelium of small intestinal villi_S100g high')] <- 'Enterocyte'
pred.scHCL[pred.scHCL %in% c('Erythroblast_Mt2 high')] <- 'Erythroblast'
pred.scHCL[pred.scHCL %in% c('Glandular epithelium_Ltf high')] <- 'Glandular epithelium'
pred.scHCL[pred.scHCL %in% c('G cell', 'Dendritic cell')] <- 'Other cells'
table(true.tags, pred.scHCL)

library(reticulate)
use_python('/local/zy/tools/anaconda3/bin/python3', required = T)
# py_config()
py_module_available('sklearn')
metrics <- import('sklearn.metrics')

accuracy <- metrics$accuracy_score(true.tags, pred.scHCL)
balanced_acc <- metrics$balanced_accuracy_score(true.tags, pred.scHCL)
balanced.accuracy.rm.unassigned <- metrics$balanced_accuracy_score(pred.scHCL, true.tags)
res.scHCL <- list()
res.scHCL$accuracy <- accuracy
res.scHCL$accuracy.rm.unassigned <- accuracy
res.scHCL$balanced.accuracy <- balanced_acc
res.scHCL$balanced.accuracy.rm.unassigned <- balanced.accuracy.rm.unassigned
path.res <- '/mdshare/node9/zy/scRef/figure/atlas_anno/'
file.res.scHCL <- paste0(path.res, 'RES_MCA_', dataset, '_scHCL.Rdata')
saveRDS(res.scHCL, file.res.scHCL)
