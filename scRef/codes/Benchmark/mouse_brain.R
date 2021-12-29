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
    meta.tag$scRef.tag <- scRef.tag
    
    true.tag[true.tag %in% unknow.cell] <- 'Unassigned'
    true.labels <- unique(true.tag)
    our.tag <- meta.tag$scRef.tag
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

source('/home/zy/my_git/scRef/main/scRef.v10.R')

############# regard counts data as reference
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
label.in <- data.frame(cell_id = row.names(label_Tasic), tag = label_Tasic$label.unlabeled.use.cols...)
exp.Tasic.sum <- .generate_ref(exp_Tasic, label.in, M='SUM')
exp_ref_mat <- exp.Tasic.sum

############### import unlabeled data
############### Habib
dataset <- 'Habib'
file.data.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat_cluster_original.txt')
# OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, del.label = c('miss'))
# saveRDS(OUT, file = paste0(path.output, dataset, '.Rdata'))
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_Habib <- OUT$data.filter
label_Habib <- OUT$label.filter

ref.names <- colnames(exp_ref_mat)
# list of cell names
all.cell <- unique(label_Habib[,1])
sc.name <- c("Astrocyte", "EndothelialCells",
             "microglia", "Neurons", "Oligodend", "OPC")
unknow.cell <- setdiff(all.cell, sc.name)
df.cell.names <- data.frame(ref.name = ref.names, sc.name = sc.name, idx = 1:length(sc.name))


#################
### scRef
source('/home/zy/my_git/scRef/main/scRef.v11.R')
setwd('~/my_git/scRef')
exp_sc_mat <- exp_Habib
result.scref <- SCREF(exp_sc_mat, exp_ref_mat, type_ref = 'count', 
                      cluster.speed = T, cluster.cell = 10,
                      min_cell = 10, CPU = 8)

meta.tag <- merge(result.scref$final.out, label_Habib, by = 'row.names')
row.names(meta.tag) <- meta.tag$Row.names
meta.tag$Row.names <- NULL
names(meta.tag) <- c('scRef.tag', 'log10Pval', 'ori.tag')

### evaluation
true.tag = meta.tag$ori.tag
scRef.tag = meta.tag$scRef.tag

simple.evaluation(true.tag, scRef.tag, df.cell.names)

### sciBet
suppressMessages(library(tidyverse))
suppressMessages(library(scibet))
suppressMessages(library(viridis))
suppressMessages(library(ggsci))
out <- .get_overlap_genes(exp_Habib, exp_Tasic)
train_set <- as.data.frame(t(out$exp_ref_mat))
train_set$label <- label_Tasic$label.unlabeled.use.cols...
test_set <- as.data.frame(t(out$exp_sc_mat))
sciBet <- SciBet(train_set, test_set)
simple.evaluation(label_Habib$label.unlabeled.use.cols..., sciBet, df.cell.names)

### singleCellNet
library(singleCellNet)
library(dplyr)
out <- .get_overlap_genes(exp_Habib, exp_Tasic)
train_set <- as.matrix(out$exp_ref_mat)
LabelsTrain <- label_Tasic
names(LabelsTrain) <- 'Annotation'
test_set <- as.matrix(out$exp_sc_mat)
class_info <- scn_train(stTrain = LabelsTrain, expTrain = train_set, dLevel = "Annotation")
classRes <- scn_predict(cnProc=class_info[['cnProc']], expDat=test_set, nrand = 50)
classRes <- classRes[, colnames(test_set)]
tags.singleCellNet <- rownames(classRes)[apply(classRes,2,which.max)]
tags.singleCellNet[tags.singleCellNet == 'rand'] <- 'Unassigned'
simple.evaluation(label_Habib$label.unlabeled.use.cols..., tags.singleCellNet, df.cell.names)


### 


############### Zeisel
dataset <- 'Zeisel'
file.data.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat_cluster_merged.txt')
OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, del.label = c('miss'))
saveRDS(OUT, file = paste0(path.output, dataset, '.Rdata'))
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_Zeisel <- OUT$data.filter
label_Zeisel <- OUT$label.filter

ref.names <- colnames(exp_ref_mat)
# list of cell names
all.cell <- unique(label_Zeisel[,1])
sc.name <- c("astrocytes_ependymal", "endothelial-mural",
             "microglia", "neurons", "oligodendrocytes", "Oligodendrocyte Precursor Cell")
unknow.cell <- setdiff(all.cell, sc.name)
df.cell.names <- data.frame(ref.name = ref.names, sc.name = sc.name, idx = 1:length(sc.name))

# scRef
source('/home/zy/my_git/scRef/main/scRef.v11.R')
setwd('~/my_git/scRef')
exp_sc_mat <- exp_Zeisel
result.scref <- SCREF(exp_sc_mat, exp_ref_mat, type_ref = 'count', 
                      cluster.speed = F, cluster.cell = 5, CPU = 8)

meta.tag <- merge(result.scref$final.out, label_Zeisel, by = 'row.names')
row.names(meta.tag) <- meta.tag$Row.names
meta.tag$Row.names <- NULL
names(meta.tag) <- c('scRef.tag', 'log10Pval', 'ori.tag')

### evaluation
true.tag = meta.tag$ori.tag
scRef.tag = meta.tag$scRef.tag

simple.evaluation(true.tag, scRef.tag, df.cell.names)

### sciBet
suppressMessages(library(tidyverse))
suppressMessages(library(scibet))
suppressMessages(library(viridis))
suppressMessages(library(ggsci))
out <- .get_overlap_genes(exp_Zeisel, exp_Tasic)
train_set <- as.data.frame(t(out$exp_ref_mat))
train_set$label <- label_Tasic$label.unlabeled.use.cols...
test_set <- as.data.frame(t(out$exp_sc_mat))
sciBet <- SciBet(train_set, test_set)
simple.evaluation(label_Zeisel$label.unlabeled.use.cols..., sciBet, df.cell.names)

### singleCellNet
library(singleCellNet)
library(dplyr)
out <- .get_overlap_genes(exp_Zeisel, exp_Tasic)
train_set <- as.matrix(out$exp_ref_mat)
LabelsTrain <- label_Tasic
names(LabelsTrain) <- 'Annotation'
test_set <- as.matrix(out$exp_sc_mat)
class_info <- scn_train(stTrain = LabelsTrain, expTrain = train_set, dLevel = "Annotation")
classRes <- scn_predict(cnProc=class_info[['cnProc']], expDat=test_set, nrand = 50)
classRes <- classRes[, colnames(test_set)]
tags.singleCellNet <- rownames(classRes)[apply(classRes,2,which.max)]
tags.singleCellNet[tags.singleCellNet == 'rand'] <- 'Unassigned'
simple.evaluation(label_Zeisel$label.unlabeled.use.cols..., tags.singleCellNet, df.cell.names)


############### YeEmilyWu
dataset <- 'YeEmilyWu'
file.data.unlabeled <- paste0(path.input, dataset, '/cell_exp.txt')
file.label.unlabeled <- paste0(path.input, dataset, '/cell_meta.txt')
# OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, del.label = c('miss'))
# saveRDS(OUT, file = paste0(path.output, dataset, '.Rdata'))
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_YeEmilyWu <- OUT$data.filter
label_YeEmilyWu <- OUT$label.filter

ref.names <- colnames(exp_ref_mat)
# list of cell names
all.cell <- unique(label_YeEmilyWu[,1])
sc.name <- c("AS", "EN", "MG", "N", "OL", "OPC")
unknow.cell <- setdiff(all.cell, sc.name)
df.cell.names <- data.frame(ref.name = ref.names, sc.name = sc.name, idx = 1:length(sc.name))

# scRef
source('/home/zy/my_git/scRef/main/scRef.v11.R')
setwd('~/my_git/scRef')
exp_sc_mat <- exp_YeEmilyWu
result.scref <- SCREF(exp_sc_mat, exp_ref_mat, type_ref = 'count', 
                      cluster.speed = T, cluster.cell = 10, CPU = 8)

meta.tag <- merge(result.scref$final.out, label_YeEmilyWu, by = 'row.names')
row.names(meta.tag) <- meta.tag$Row.names
meta.tag$Row.names <- NULL
names(meta.tag) <- c('scRef.tag', 'log10Pval', 'ori.tag')

### evaluation
true.tag = meta.tag$ori.tag
scRef.tag = meta.tag$scRef.tag

simple.evaluation(true.tag, scRef.tag, df.cell.names)

### sciBet
suppressMessages(library(tidyverse))
suppressMessages(library(scibet))
suppressMessages(library(viridis))
suppressMessages(library(ggsci))
out <- .get_overlap_genes(exp_YeEmilyWu, exp_Tasic)
train_set <- as.data.frame(t(out$exp_ref_mat))
train_set$label <- label_Tasic$label.unlabeled.use.cols...
test_set <- as.data.frame(t(out$exp_sc_mat))
sciBet <- SciBet(train_set, test_set)
simple.evaluation(label_YeEmilyWu$label.unlabeled.use.cols..., sciBet, df.cell.names)

### singleCellNet
library(singleCellNet)
library(dplyr)
out <- .get_overlap_genes(exp_YeEmilyWu, exp_Tasic)
train_set <- as.matrix(out$exp_ref_mat)
LabelsTrain <- label_Tasic
names(LabelsTrain) <- 'Annotation'
test_set <- as.matrix(out$exp_sc_mat)
class_info <- scn_train(stTrain = LabelsTrain, expTrain = train_set, dLevel = "Annotation")
classRes <- scn_predict(cnProc=class_info[['cnProc']], expDat=test_set, nrand = 50)
classRes <- classRes[, colnames(test_set)]
tags.singleCellNet <- rownames(classRes)[apply(classRes,2,which.max)]
tags.singleCellNet[tags.singleCellNet == 'rand'] <- 'Unassigned'
simple.evaluation(label_YeEmilyWu$label.unlabeled.use.cols..., tags.singleCellNet, df.cell.names)


############### Mizrak
dataset <- 'Mizrak'
######################
library(stringr)
file.mtx.rep1 <- paste0(path.input, dataset, '/GSE109447_29319_cells.matrix.txt')
mtx.rep1 <- read.delim(file.mtx.rep1, header = F, stringsAsFactors = F)
file.cellid.rep1 <- paste0(path.input, dataset, '/GSE109447_29319_cells_id_repinfo.txt')
df.cellid <- read.table(file.cellid.rep1, stringsAsFactors = F, sep = '\t')
df.cellid <- str_replace_all(df.cellid, '_', '.')
df.cellid <- str_replace_all(df.cellid, '-', '.')
file.cells.rep1 <- paste0(path.input, dataset, '/GSE109447_Rep1_29319cells_Basic.txt')
df.cells <- read.table(file.cells.rep1, stringsAsFactors = F, sep = '\t')

data.unlabeled <- mtx.rep1
data.genes <- data.unlabeled[, 'V2']
data.unlabeled[, 'V1'] <- NULL
data.unlabeled[, 'V2'] <- NULL
data.unlabeled <- as.matrix(data.unlabeled)
rownames(data.unlabeled) <- data.genes
colnames(data.unlabeled) <- df.cellid$V1
# read label file
label.unlabeled <- data.frame(Cluster = df.cells$V1, row.names = df.cellid$V1)

OUT <- list()
OUT$data.filter <- data.unlabeled
OUT$label.filter <- label.unlabeled
saveRDS(OUT, file = paste0(path.output, dataset, '.Rdata'))
##########

OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_Mizrak <- OUT$data.filter
label_Mizrak <- OUT$label.filter

ref.names <- colnames(exp_ref_mat)
# list of cell names
all.cell <- unique(label_Mizrak[,1])
sc.name <- c("Astrocyte", "Endothelial", "Microglia", "Neuron", "Oligodendrocyte", "OPC")
unknow.cell <- setdiff(all.cell, sc.name)
df.cell.names <- data.frame(ref.name = ref.names, sc.name = sc.name, idx = 1:length(sc.name))

# scRef
source('/home/zy/my_git/scRef/main/scRef.v11.R')
setwd('~/my_git/scRef')
exp_sc_mat <- exp_Mizrak
result.scref <- SCREF(exp_sc_mat, exp_ref_mat, type_ref = 'count', 
                      cluster.speed = T, cluster.cell = 10, CPU = 8)

meta.tag <- merge(result.scref$final.out, label_Mizrak, by = 'row.names')
row.names(meta.tag) <- meta.tag$Row.names
meta.tag$Row.names <- NULL
names(meta.tag) <- c('scRef.tag', 'log10Pval', 'ori.tag')

### evaluation
true.tag = meta.tag$ori.tag
scRef.tag = meta.tag$scRef.tag

simple.evaluation(true.tag, scRef.tag, df.cell.names)

### sciBet
suppressMessages(library(tidyverse))
suppressMessages(library(scibet))
suppressMessages(library(viridis))
suppressMessages(library(ggsci))
out <- .get_overlap_genes(exp_Mizrak, exp_Tasic)
train_set <- as.data.frame(t(out$exp_ref_mat))
train_set$label <- label_Tasic$label.unlabeled.use.cols...
test_set <- as.data.frame(t(out$exp_sc_mat))
sciBet <- SciBet(train_set, test_set)
simple.evaluation(label_Mizrak$Cluster, sciBet, df.cell.names)

### singleCellNet
library(singleCellNet)
library(dplyr)
out <- .get_overlap_genes(exp_Mizrak, exp_Tasic)
train_set <- as.matrix(out$exp_ref_mat)
LabelsTrain <- label_Tasic
names(LabelsTrain) <- 'Annotation'
test_set <- as.matrix(out$exp_sc_mat)
class_info <- scn_train(stTrain = LabelsTrain, expTrain = train_set, dLevel = "Annotation")
classRes <- scn_predict(cnProc=class_info[['cnProc']], expDat=test_set, nrand = 50)
classRes <- classRes[, colnames(test_set)]
tags.singleCellNet <- rownames(classRes)[apply(classRes,2,which.max)]
tags.singleCellNet[tags.singleCellNet == 'rand'] <- 'Unassigned'
simple.evaluation(label_Mizrak$Cluster, tags.singleCellNet, df.cell.names)



############# regard counts data as reference
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
label.in <- data.frame(cell_id = row.names(label_Tasic), tag = label_Tasic$label.unlabeled.use.cols...)
exp.Tasic.sum <- .generate_ref(exp_Tasic, label.in, M='SUM')
exp_ref_mat <- exp.Tasic.sum

file.ref <- '/home/disk/scRef/MouseAtlas_SingleCell_Han2018/MCA_ref/Brain_ref_mouse.txt'
df.MCA.Brain <- read.delim(file.ref, header = F)
cell.names <- as.vector(t(df.MCA.Brain[1, 1:(dim(df.MCA.Brain)[2])-1]))
df.MCA.Brain <- read.delim(file.ref, row.names = 1)
colnames(df.MCA.Brain) <- cell.names

# sum counts
list.split <- strsplit(cell.names, '(', fixed=T)
vec.cellname <- c()
for (i in 1:length(list.split)) {
    vec.cellname <- c(vec.cellname, list.split[[i]][1])
}
list.split <- strsplit(vec.cellname, '_', fixed=T)
cell_names <- c()
for (i in 1:length(list.split)) {
    cell_names <- c(cell_names, list.split[[i]][1])
}

df.MCA.Brain.sum <- data.frame(stringsAsFactors = F)
gene.names <- rownames(df.MCA.Brain)
i = 1
for (cell in unique(cell_names)) {
    df.sub <- df.MCA.Brain[, cell_names == cell]
    if (class(df.sub) == 'data.frame') {
        df.sum <- rowSums(df.sub)
    } else {
        df.sum <- df.sub
    }
    df.sum <- data.frame(df.sum, row.names = gene.names)
    names(df.sum) <- cell
    if (i ==1) {
        df.MCA.Brain.sum <- df.sum
    } else {
        df.MCA.Brain.sum <- cbind(df.MCA.Brain.sum, df.sum)
    }
    i = i+ 1
}


exp_ref_mat <- df.MCA.Brain.sum



############### import unlabeled data
############### Habib
dataset <- 'Habib'
file.data.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat_cluster_original.txt')
# OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, del.label = c('miss'))
# saveRDS(OUT, file = paste0(path.output, dataset, '.Rdata'))
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_Habib <- OUT$data.filter
label_Habib <- OUT$label.filter

ref.names <- colnames(exp_ref_mat)
# list of cell names
all.cell <- unique(label_Habib[,1])
sc.name <- c("Unassigned", "Unassigned", "Unassigned", "microglia", "Unassigned",
             "Unassigned", "Astrocyte", "Ependymocytes", "Oligodend", "Neurons", "OPC")
unknow.cell <- setdiff(all.cell, sc.name)
df.cell.names <- data.frame(ref.name = ref.names, sc.name = sc.name, idx = 1:length(sc.name))


#################
### scRef
source('/home/zy/my_git/scRef/main/scRef.v11.R')
setwd('~/my_git/scRef')
exp_sc_mat <- exp_Habib
result.scref <- SCREF(exp_sc_mat, exp_ref_mat, type_ref = 'count', 
                      cluster.speed = T, cluster.cell = 10,
                      min_cell = 10, CPU = 8)

meta.tag <- merge(result.scref$final.out, label_Habib, by = 'row.names')
row.names(meta.tag) <- meta.tag$Row.names
meta.tag$Row.names <- NULL
names(meta.tag) <- c('scRef.tag', 'log10Pval', 'ori.tag')

### evaluation
true.tag = meta.tag$ori.tag
scRef.tag = meta.tag$scRef.tag

simple.evaluation(true.tag, scRef.tag, df.cell.names)


### sciBet
suppressMessages(library(tidyverse))
suppressMessages(library(scibet))
suppressMessages(library(viridis))
suppressMessages(library(ggsci))
out <- .get_overlap_genes(exp_Habib, exp_ref_mat)
train_set <- as.data.frame(t(out$exp_ref_mat))
train_set$label <- colnames(exp_ref_mat)
test_set <- as.data.frame(t(out$exp_sc_mat))
sciBet <- SciBet(train_set, test_set)
simple.evaluation(label_Habib$label.unlabeled.use.cols..., sciBet, df.cell.names)
