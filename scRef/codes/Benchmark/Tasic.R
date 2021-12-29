library(reticulate)
# library(metap)

source('/home/zy/my_git/scRef/main/scRef.v5.R')
# import python package: sklearn.metrics
use_python('/home/zy/tools/anaconda3/bin/python3', required = T)
# py_config()
py_module_available('sklearn')
metrics <- import('sklearn.metrics')

setwd('/home/zy/scRef/try_data')
# input file
file.ref <- './scRef/Reference/MouseBrain_Bulk_Zhang2014/Reference_expression.txt'
# parameters
num.cpu <- 20

# reference file
file.ref <- file.ref
exp_ref_mat.origin <- read.table(file.ref, header=T, row.names=1, sep='\t', check.name=F)
# MCA cell name
ref.names <- names(exp_ref_mat.origin)
# ref.names <- c("Astrocyte", "Neuron", "Oligodendrocyte precursor cell",
#                "Oligodendrocyte", "Myelinating oligodendrocyte",
#                "Microglia", "Endothelial cell")
# names(exp_ref_mat.origin) <- ref.names

# function of data preparation
prepare.data <- function(file.data.unlabeled, file.label.unlabeled, 
                         del.label = c('Unclassified', 'OTHER')) {
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


######################### Tasic
dataset <- 'Tasic'
# folder saving results
path.out <- paste0('/home/zy/scRef/try_data/', dataset)
if (!file.exists(path.out)) {
    dir.create(path.out)
}
file.data.unlabeled <- paste0('./summary/', dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0('./summary/', dataset, '_exp_sc_mat_cluster_merged.txt')
# OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled)
# saveRDS(OUT, file = './Tasic/Tasic.Rdata')
OUT <- readRDS('./Tasic/Tasic.Rdata')
data.filter <- OUT$data.filter
label.filter <- OUT$label.filter
# list of cell names
all.cell <- unique(label.filter[,1])
sc.name <- c("Astrocyte", "Neuron", "Oligodendrocyte Precursor Cell",
             "Oligodendrocyte", "Oligodendrocyte",
             "Microglia", "Endothelial Cell")
unknow.cell <- setdiff(all.cell, sc.name)
df.cell.names <- data.frame(ref.name = ref.names, sc.name = sc.name)

exp_ref_mat <- exp_ref_mat.origin
exp_sc_mat = data.filter

# run scRef
result.scref <- SCREF(exp_sc_mat, exp_ref_mat, type_ref = 'fpkm', min_cell=5, CPU = num.cpu)
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
vec.cutoff <- 1:90

one.eval <- function(cutoff, meta.tag) {
    df.sub <- data.frame(stringsAsFactors = F)
    # library(reticulate)
    # use_python('/home/zy/tools/anaconda3/bin/python3', required = T)
    # metrics <- import('sklearn.metrics')
    true.tag <- meta.tag$ori.tag
    true.tag[true.tag %in% unknow.cell] <- 'unknown'
    our.tag <- meta.tag$scRef.tag
    our.tag[meta.tag$log10Pval <= cutoff] <- 'unknown'
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


