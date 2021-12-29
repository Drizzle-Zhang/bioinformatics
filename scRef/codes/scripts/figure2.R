library(reticulate)
# import python package: sklearn.metrics
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

source('/home/zy/my_git/scRef/main/scRef.v7.R')

# import unlabeled data
# Habib
setwd('/home/zy/scRef/')
dataset <- 'Habib'
file.data.unlabeled <- paste0('./summary/', dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0('./summary/', dataset, '_exp_sc_mat_cluster_original.txt')
# OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, del.label = c('miss'))
# 
# saveRDS(OUT, file = './Habib/Habib.Rdata')
OUT <- readRDS('./Habib/Habib.Rdata')
data.filter <- OUT$data.filter
label.filter <- OUT$label.filter
exp_sc_mat <- data.filter
label.Habib <- label.filter

######## regard counts data as reference
setwd('/home/zy/scRef/')
dataset <- 'Tasic'
file.data.unlabeled <- paste0('./summary/', dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0('./summary/', dataset, '_exp_sc_mat_cluster_merged.txt')
# OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled)
# saveRDS(OUT, file = './Tasic/Tasic.Rdata')
OUT <- readRDS('./Tasic/Tasic.Rdata')
data.filter <- OUT$data.filter
label.filter <- OUT$label.filter
label.in <- data.frame(cell_id = row.names(label.filter), tag = label.filter$label.unlabeled.use.cols...)
exp.Tasic.sum <- .generate_ref(data.filter, label.in, M='SUM')

# run scRef
exp_ref_mat <- exp.Tasic.sum
source('/home/zy/my_git/scRef/main/scRef.v7.R')
setwd('~/my_git/scRef')
result.scref <- SCREF(exp_sc_mat, exp_ref_mat, type_ref = 'count', 
                      cluster.speed = T, cluster.cell = 10,
                      min_cell = 10, CPU = 8)

meta.tag <- merge(result.scref$final.out, label.filter, by = 'row.names')
row.names(meta.tag) <- meta.tag$Row.names
meta.tag$Row.names <- NULL
names(meta.tag) <- c('scRef.tag', 'log10Pval', 'ori.tag')

### evaluation
ori.tag = meta.tag$ori.tag
scRef.tag = meta.tag$scRef.tag

ref.names <- colnames(exp_ref_mat)
# list of cell names
all.cell <- unique(label.filter[,1])
sc.name <- c("Astrocyte", "EndothelialCells",
             "microglia", "Neurons", "Oligodend", "OPC")
unknow.cell <- setdiff(all.cell, sc.name)
df.cell.names <- data.frame(ref.name = ref.names, sc.name = sc.name, idx = 1:length(sc.name))

# uniform tags
for (j in 1:dim(df.cell.names)[1]) {
    scRef.tag[scRef.tag == df.cell.names[j, 'ref.name']] <- 
        df.cell.names[j, 'sc.name']
}
meta.tag$scRef.tag <- scRef.tag

# default cutoff
true.tag <- meta.tag$ori.tag
true.tag[true.tag %in% unknow.cell] <- 'Unassigned'
our.tag <- meta.tag$scRef.tag
metrics$f1_score(true.tag, our.tag, average = 'weighted')
metrics$f1_score(true.tag, our.tag, average = 'macro')
metrics$accuracy_score(true.tag, our.tag)





