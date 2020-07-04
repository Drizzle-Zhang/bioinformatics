library(stringr)
library(reticulate)

Sys.setenv(RETICULATE_PYTHON = "/home/zy/tools/anaconda3/bin/python3")
source('/home/zy/my_git/bioinformatics/scRef/try/scRef.R')

setwd('/home/zy/scRef/try_data')
# input file
file.ref <- './scRef/Reference/MouseBrain_Bulk_Zhang2014/Reference_expression.txt'
# file.data <- "./summary/Tasic_exp_sc_mat.txt"
# file.label.unlabeled <- './summary/Tasic_exp_sc_mat_cluster_merged.txt'
# parameters
num.cpu <- 10

# reference file
file.ref <- file.ref
exp_ref_mat.origin <- read.table(file.ref, header=T, row.names=1, sep='\t', check.name=F)
# MCA cell name
names(exp_ref_mat.origin) <- c("Astrocyte", "Neuron", "Oligodendrocyte precursor cell",
                               "Newly Formed Oligodendrocyte", "Myelinating oligodendrocyte",
                               "Microglia", "Endothelial cell")

######################### delete specific cell
cell.delete <- "Astrocyte"
exp_ref_mat.origin <- exp_ref_mat.origin[, setdiff(names(exp_ref_mat.origin), c(cell.delete))]

out.markers <- find.markers(exp_ref_mat.origin)
list.cell.genes <- out.markers[['list.cell.genes']]
exp_ref_mat.origin <- out.markers[['exp_ref_mat']]

######################### unlabeled data
######################################################################################
######################### Zeisel
file.data <- "./summary/Zeisel_exp_sc_mat.txt"
file.label.unlabeled <- './summary/Zeisel_exp_sc_mat_cluster_merged.txt'
file.data <- file.data
data.unlabeled <- read.delim(file.data, row.names=1)
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
use.cols <- row.names(label.unlabeled)[label.unlabeled[,1] != 'Unclassified']
data.filter <- data.unlabeled[,use.cols]
label.filter <- data.frame(label.unlabeled[use.cols,], row.names = use.cols)

# get overlap genes
exp_sc_mat = data.filter
exp_ref_mat <- exp_ref_mat.origin
exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
gene_sc=rownames(exp_sc_mat)
gene_ref=rownames(exp_ref_mat)
gene_over= gene_sc[which(gene_sc %in% gene_ref)]
exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
print('Number of overlapped genes:')
print(nrow(exp_sc_mat))

######################### run scRef
result.scref <- SCREF(exp_sc_mat, exp_ref_mat, CPU = num.cpu)
tag <- result.scref$tag2

# confirm label
exp_sc_mat <- exp_sc_mat[gene_over,]
ori.tag = label.filter[names(exp_sc_mat), 1]
scRef.tag = tag[,2]
# method.test <- 'wilcox'
method.test <- 't.test'
# method.test <- 'oneway_test'
meta.tag <- comfirm.label(exp_sc_mat, ori.tag, scRef.tag, method.test)

# evaluation
# import python package: sklearn.metrics
use_condaenv("/home/zy/tools/anaconda3")
# py_config()
metrics <- import('sklearn.metrics')
# uniform tags
scRef.tag[scRef.tag == "Astrocyte"] <- "astrocytes_ependymal"
scRef.tag[scRef.tag == "Newly Formed Oligodendrocyte"] <- "oligodendrocytes"
scRef.tag[scRef.tag == "Myelinating oligodendrocyte"] <- "oligodendrocytes"
scRef.tag[scRef.tag == "Endothelial cell"] <- "endothelial-mural" 
scRef.tag[scRef.tag == "Neuron"] <- "neurons" 
scRef.tag[scRef.tag == "Microglia"] <- "microglia" 
cell.delete <- "astrocytes_ependymal"
meta.tag$scRef.tag <- scRef.tag
vec.cutoff <- c(seq(0.005, 0.1, 0.005), seq(0.15, 0.95, 0.05))
df.metrics <- data.frame(cutoff = vec.cutoff, weighted.f1 = rep(0, length(vec.cutoff)))
for (i in 1:length(vec.cutoff)) {
    cutoff <- vec.cutoff[i]
    true.tag <- meta.tag$ori.tag
    true.tag[true.tag == cell.delete] <- 'unknown'
    our.tag <- meta.tag$scRef.tag
    our.tag[meta.tag$qvalue > cutoff] <- 'unknown'
    sub.weighted.f1 <- metrics$f1_score(true.tag, our.tag, average = 'weighted')
    df.metrics[i, 'weighted.f1'] <- sub.weighted.f1
}

# best cutoff
best.cutoff <- df.metrics$cutoff[df.metrics$weighted.f1 == max(df.metrics$weighted.f1)][1]
new.tag <- meta.tag$scRef.tag
new.tag[meta.tag$qvalue > best.cutoff] <- 'unknown'
meta.tag$new.tag <- new.tag
print(best.cutoff)

# no scRef plus
true.tag <- meta.tag$ori.tag
true.tag[true.tag == cell.delete] <- 'unknown'
meta.tag$ori.tag <- true.tag
metrics$f1_score(true.tag, scRef.tag, average = 'weighted')
metrics$f1_score(true.tag, new.tag, average = 'weighted')

# remove neuron
meta.tag.remove <- meta.tag[meta.tag$ori.tag != "neurons", c("ori.tag", "scRef.tag", "qvalue")]
vec.cutoff <- c(seq(0.005, 0.1, 0.005), seq(0.15, 0.95, 0.05))
df.metrics.remove <- data.frame(cutoff = vec.cutoff, weighted.f1 = rep(0, length(vec.cutoff)))
for (i in 1:length(vec.cutoff)) {
    cutoff <- vec.cutoff[i]
    true.tag <- meta.tag.remove$ori.tag
    true.tag[true.tag == cell.delete] <- 'unknown'
    our.tag <- meta.tag.remove$scRef.tag
    our.tag[meta.tag.remove$qvalue > cutoff] <- 'unknown'
    sub.weighted.f1 <- metrics$f1_score(true.tag, our.tag, average = 'weighted')
    df.metrics.remove[i, 'weighted.f1'] <- sub.weighted.f1
}
# best cutoff
best.cutoff.remove <- df.metrics.remove$cutoff[df.metrics.remove$weighted.f1 == max(df.metrics.remove$weighted.f1)][1]
new.tag.remove <- meta.tag.remove$scRef.tag
new.tag.remove[meta.tag.remove$qvalue > best.cutoff.remove] <- 'unknown'
meta.tag.remove$new.tag <- new.tag.remove
print(best.cutoff.remove)

# no scRef plus
true.tag.remove <- meta.tag.remove$ori.tag
true.tag.remove[true.tag.remove == cell.delete] <- 'unknown'
meta.tag.remove$ori.tag <- true.tag.remove
metrics$f1_score(true.tag.remove, meta.tag.remove$scRef.tag, average = 'weighted')
metrics$f1_score(true.tag.remove, new.tag.remove, average = 'weighted')


######################################################################################
######################### Tasic
file.data <- "./summary/Tasic_exp_sc_mat.txt"
file.label.unlabeled <- './summary/Tasic_exp_sc_mat_cluster_merged.txt'
file.data <- file.data
data.unlabeled <- read.delim(file.data, row.names=1)
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
use.cols <- row.names(label.unlabeled)[label.unlabeled[,1] != 'Unclassified']
data.filter <- data.unlabeled[,use.cols]
label.filter <- data.frame(label.unlabeled[use.cols,], row.names = use.cols)

# get overlap genes
exp_sc_mat = data.filter
exp_ref_mat <- exp_ref_mat.origin
exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
gene_sc=rownames(exp_sc_mat)
gene_ref=rownames(exp_ref_mat)
gene_over= gene_sc[which(gene_sc %in% gene_ref)]
exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
print('Number of overlapped genes:')
print(nrow(exp_sc_mat))

######################### run scRef
result.scref <- SCREF(exp_sc_mat, exp_ref_mat, CPU = num.cpu)
tag <- result.scref$tag2

# confirm label
exp_sc_mat <- exp_sc_mat[gene_over,]
ori.tag = label.filter[names(exp_sc_mat), 1]
scRef.tag = tag[,2]
method.test <- 'wilcox'
# method.test <- 't.test'
# method.test <- 'oneway_test'
meta.tag <- comfirm.label(exp_sc_mat, ori.tag, scRef.tag, method.test)

# evaluation
# import python package: sklearn.metrics
use_condaenv("/home/zy/tools/anaconda3")
# py_config()
metrics <- import('sklearn.metrics')
# uniform tags
scRef.tag[scRef.tag == "Newly Formed Oligodendrocyte"] <- "Oligodendrocyte"
scRef.tag[scRef.tag == "Myelinating oligodendrocyte"] <- "Oligodendrocyte"
scRef.tag[scRef.tag == "Endothelial cell"] <- "Endothelial Cell" 
scRef.tag[scRef.tag == "Oligodendrocyte precursor cell"] <- "Oligodendrocyte Precursor Cell" 
cell.delete <- "Astrocyte" 
meta.tag$scRef.tag <- scRef.tag
vec.cutoff <- c(seq(0.005, 0.1, 0.005), seq(0.15, 0.95, 0.05))
df.metrics <- data.frame(cutoff = vec.cutoff, weighted.f1 = rep(0, length(vec.cutoff)))
for (i in 1:length(vec.cutoff)) {
    cutoff <- vec.cutoff[i]
    true.tag <- meta.tag$ori.tag
    true.tag[true.tag == cell.delete] <- 'unknown'
    our.tag <- meta.tag$scRef.tag
    our.tag[meta.tag$qvalue > cutoff] <- 'unknown'
    sub.weighted.f1 <- metrics$f1_score(true.tag, our.tag, average = 'weighted')
    df.metrics[i, 'weighted.f1'] <- sub.weighted.f1
}

# best cutoff
best.cutoff <- df.metrics$cutoff[df.metrics$weighted.f1 == max(df.metrics$weighted.f1)][1]
new.tag <- meta.tag$scRef.tag
new.tag[meta.tag$qvalue > best.cutoff] <- 'unknown'
meta.tag$new.tag <- new.tag
print(best.cutoff)

# no scRef plus
true.tag <- meta.tag$ori.tag
true.tag[true.tag == cell.delete] <- 'unknown'
meta.tag$ori.tag <- true.tag
metrics$f1_score(true.tag, scRef.tag, average = 'weighted')
metrics$f1_score(true.tag, new.tag, average = 'weighted')

# remove neuron
meta.tag.remove <- meta.tag[meta.tag$ori.tag != "Neuron", c("ori.tag", "scRef.tag", "qvalue")]
vec.cutoff <- c(seq(0.005, 0.1, 0.005), seq(0.15, 0.95, 0.05))
df.metrics.remove <- data.frame(cutoff = vec.cutoff, weighted.f1 = rep(0, length(vec.cutoff)))
for (i in 1:length(vec.cutoff)) {
    cutoff <- vec.cutoff[i]
    true.tag <- meta.tag.remove$ori.tag
    true.tag[true.tag == cell.delete] <- 'unknown'
    our.tag <- meta.tag.remove$scRef.tag
    our.tag[meta.tag.remove$qvalue > cutoff] <- 'unknown'
    sub.weighted.f1 <- metrics$f1_score(true.tag, our.tag, average = 'weighted')
    df.metrics.remove[i, 'weighted.f1'] <- sub.weighted.f1
}
# best cutoff
best.cutoff.remove <- df.metrics.remove$cutoff[df.metrics.remove$weighted.f1 == max(df.metrics.remove$weighted.f1)][1]
new.tag.remove <- meta.tag.remove$scRef.tag
new.tag.remove[meta.tag.remove$qvalue > best.cutoff.remove] <- 'unknown'
meta.tag.remove$new.tag <- new.tag.remove
print(best.cutoff.remove)

# no scRef plus
true.tag.remove <- meta.tag.remove$ori.tag
true.tag.remove[true.tag.remove == cell.delete] <- 'unknown'
meta.tag.remove$ori.tag <- true.tag.remove
metrics$f1_score(true.tag.remove, meta.tag.remove$scRef.tag, average = 'weighted')
metrics$f1_score(true.tag.remove, new.tag.remove, average = 'weighted')


######################################################################################
######################### Habib
file.data <- "./summary/Habib_exp_sc_mat.txt"
file.label.unlabeled <- './summary/Habib_exp_sc_mat_cluster_merged.txt'
file.data <- file.data
data.unlabeled <- read.delim(file.data, row.names=1)
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
use.cols <- row.names(label.unlabeled)[label.unlabeled[,1] != 'OTHER']
data.filter <- data.unlabeled[,use.cols]
label.filter <- data.frame(label.unlabeled[use.cols,], row.names = use.cols)

# get overlap genes
exp_sc_mat = data.filter
exp_ref_mat <- exp_ref_mat.origin
exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
gene_sc=rownames(exp_sc_mat)
gene_ref=rownames(exp_ref_mat)
gene_over= gene_sc[which(gene_sc %in% gene_ref)]
exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
print('Number of overlapped genes:')
print(nrow(exp_sc_mat))

######################### run scRef
result.scref <- SCREF(exp_sc_mat, exp_ref_mat, CPU = num.cpu)
tag <- result.scref$tag2

# confirm label
exp_sc_mat <- exp_sc_mat[gene_over,]
ori.tag = label.filter[names(exp_sc_mat), 1]
scRef.tag = tag[,2]
# method.test <- 'wilcox'
# method.test <- 't.test'
method.test <- 'oneway_test'
meta.tag <- comfirm.label(exp_sc_mat, ori.tag, scRef.tag, method.test)

# evaluation
# import python package: sklearn.metrics
use_condaenv("/home/zy/tools/anaconda3")
# py_config()
metrics <- import('sklearn.metrics')
# uniform tags
scRef.tag[scRef.tag == "Oligodendrocyte precursor cell"] <- "OPC" 
scRef.tag[scRef.tag == "Newly Formed Oligodendrocyte"] <- "Oligodend"
scRef.tag[scRef.tag == "Myelinating oligodendrocyte"] <- "Oligodend"
scRef.tag[scRef.tag == "Endothelial cell"] <- "EndothelialCells" 
scRef.tag[scRef.tag == "Neuron"] <- "Neurons" 
scRef.tag[scRef.tag == "Microglia"] <- "microglia" 
cell.delete <- "Astrocyte"
meta.tag$scRef.tag <- scRef.tag
vec.cutoff <- c(seq(0.005, 0.1, 0.005), seq(0.15, 0.95, 0.05))
df.metrics <- data.frame(cutoff = vec.cutoff, weighted.f1 = rep(0, length(vec.cutoff)))
for (i in 1:length(vec.cutoff)) {
    cutoff <- vec.cutoff[i]
    true.tag <- meta.tag$ori.tag
    true.tag[true.tag == cell.delete] <- 'unknown'
    our.tag <- meta.tag$scRef.tag
    our.tag[meta.tag$qvalue > cutoff] <- 'unknown'
    sub.weighted.f1 <- metrics$f1_score(true.tag, our.tag, average = 'weighted')
    df.metrics[i, 'weighted.f1'] <- sub.weighted.f1
}

# best cutoff
best.cutoff <- df.metrics$cutoff[df.metrics$weighted.f1 == max(df.metrics$weighted.f1)][1]
new.tag <- meta.tag$scRef.tag
new.tag[meta.tag$qvalue > best.cutoff] <- 'unknown'
meta.tag$new.tag <- new.tag
print(best.cutoff)

# no scRef plus
true.tag <- meta.tag$ori.tag
true.tag[true.tag == cell.delete] <- 'unknown'
meta.tag$ori.tag <- true.tag
metrics$f1_score(true.tag, scRef.tag, average = 'weighted')
metrics$f1_score(true.tag, new.tag, average = 'weighted')

# remove neuron
meta.tag.remove <- meta.tag[meta.tag$ori.tag != "Neurons", c("ori.tag", "scRef.tag", "qvalue")]
vec.cutoff <- c(seq(0.005, 0.1, 0.005), seq(0.15, 0.95, 0.05))
df.metrics.remove <- data.frame(cutoff = vec.cutoff, weighted.f1 = rep(0, length(vec.cutoff)))
for (i in 1:length(vec.cutoff)) {
    cutoff <- vec.cutoff[i]
    true.tag <- meta.tag.remove$ori.tag
    true.tag[true.tag == cell.delete] <- 'unknown'
    our.tag <- meta.tag.remove$scRef.tag
    our.tag[meta.tag.remove$qvalue > cutoff] <- 'unknown'
    sub.weighted.f1 <- metrics$f1_score(true.tag, our.tag, average = 'weighted')
    df.metrics.remove[i, 'weighted.f1'] <- sub.weighted.f1
}
# best cutoff
best.cutoff.remove <- df.metrics.remove$cutoff[df.metrics.remove$weighted.f1 == max(df.metrics.remove$weighted.f1)][1]
new.tag.remove <- meta.tag.remove$scRef.tag
new.tag.remove[meta.tag.remove$qvalue > best.cutoff.remove] <- 'unknown'
meta.tag.remove$new.tag <- new.tag.remove
print(best.cutoff.remove)

# no scRef plus
true.tag.remove <- meta.tag.remove$ori.tag
true.tag.remove[true.tag.remove == cell.delete] <- 'unknown'
meta.tag.remove$ori.tag <- true.tag.remove
metrics$f1_score(true.tag.remove, meta.tag.remove$scRef.tag, average = 'weighted')
metrics$f1_score(true.tag.remove, new.tag.remove, average = 'weighted')


######################################################################################
######################### Hochgerner
file.data <- "./summary/Hochgerner_exp_sc_mat.txt"
file.label.unlabeled <- './summary/Hochgerner_exp_sc_mat_cluster_merged.txt'
file.data <- file.data
data.unlabeled <- read.delim(file.data, row.names=1)
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
use.cols <- row.names(label.unlabeled)[label.unlabeled[,1] != 'OTHER']
data.filter <- data.unlabeled[,use.cols]
label.filter <- data.frame(label.unlabeled[use.cols,], row.names = use.cols)

# get overlap genes
exp_sc_mat = data.filter
exp_ref_mat <- exp_ref_mat.origin
exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
gene_sc=rownames(exp_sc_mat)
gene_ref=rownames(exp_ref_mat)
gene_over= gene_sc[which(gene_sc %in% gene_ref)]
exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
print('Number of overlapped genes:')
print(nrow(exp_sc_mat))

######################### run scRef
result.scref <- SCREF(exp_sc_mat, exp_ref_mat, CPU = num.cpu)
tag <- result.scref$tag2

# confirm label
exp_sc_mat <- exp_sc_mat[gene_over,]
ori.tag = label.filter[names(exp_sc_mat), 1]
scRef.tag = tag[,2]
# method.test <- 'wilcox'
# method.test <- 't.test'
method.test <- 'oneway_test'
meta.tag <- comfirm.label(exp_sc_mat, ori.tag, scRef.tag, method.test)

# evaluation
# import python package: sklearn.metrics
use_condaenv("/home/zy/tools/anaconda3")
# py_config()
metrics <- import('sklearn.metrics')
# uniform tags
scRef.tag[scRef.tag == "Oligodendrocyte precursor cell"] <- "Oligodendrocyte.Precursor.Cell"
scRef.tag[scRef.tag == "Newly Formed Oligodendrocyte"] <- "Oligodendrocyte"
scRef.tag[scRef.tag == "Myelinating oligodendrocyte"] <- "Oligodendrocyte"
scRef.tag[scRef.tag == "Endothelial cell"] <- "Endothelial.Cells" 
scRef.tag[scRef.tag == "Neuron"] <- "Neuron" 
scRef.tag[scRef.tag == "Microglia"] <- "Microglia" 
cell.delete <- "Astrocytes"
meta.tag$scRef.tag <- scRef.tag
vec.cutoff <- c(seq(0.005, 0.1, 0.005), seq(0.11, 1, 0.01))
df.metrics <- data.frame(cutoff = vec.cutoff, weighted.f1 = rep(0, length(vec.cutoff)))
for (i in 1:length(vec.cutoff)) {
    cutoff <- vec.cutoff[i]
    true.tag <- meta.tag$ori.tag
    true.tag[true.tag == cell.delete] <- 'unknown'
    our.tag <- meta.tag$scRef.tag
    our.tag[meta.tag$qvalue > cutoff] <- 'unknown'
    sub.weighted.f1 <- metrics$f1_score(true.tag, our.tag, average = 'weighted')
    df.metrics[i, 'weighted.f1'] <- sub.weighted.f1
}

# best cutoff
best.cutoff <- df.metrics$cutoff[df.metrics$weighted.f1 == max(df.metrics$weighted.f1)][1]
new.tag <- meta.tag$scRef.tag
new.tag[meta.tag$qvalue > best.cutoff] <- 'unknown'
meta.tag$new.tag <- new.tag
print(best.cutoff)

# no scRef plus
true.tag <- meta.tag$ori.tag
true.tag[true.tag == cell.delete] <- 'unknown'
meta.tag$ori.tag <- true.tag
metrics$f1_score(true.tag, scRef.tag, average = 'weighted')
metrics$f1_score(true.tag, new.tag, average = 'weighted')

# remove neuron
meta.tag.remove <- meta.tag[meta.tag$ori.tag != "Neuron", c("ori.tag", "scRef.tag", "qvalue")]
vec.cutoff <- c(seq(0.005, 0.1, 0.005), seq(0.15, 0.95, 0.05))
df.metrics.remove <- data.frame(cutoff = vec.cutoff, weighted.f1 = rep(0, length(vec.cutoff)))
for (i in 1:length(vec.cutoff)) {
    cutoff <- vec.cutoff[i]
    true.tag <- meta.tag.remove$ori.tag
    true.tag[true.tag == cell.delete] <- 'unknown'
    our.tag <- meta.tag.remove$scRef.tag
    our.tag[meta.tag.remove$qvalue > cutoff] <- 'unknown'
    sub.weighted.f1 <- metrics$f1_score(true.tag, our.tag, average = 'weighted')
    df.metrics.remove[i, 'weighted.f1'] <- sub.weighted.f1
}
# best cutoff
best.cutoff.remove <- df.metrics.remove$cutoff[df.metrics.remove$weighted.f1 == max(df.metrics.remove$weighted.f1)][1]
new.tag.remove <- meta.tag.remove$scRef.tag
new.tag.remove[meta.tag.remove$qvalue > best.cutoff.remove] <- 'unknown'
meta.tag.remove$new.tag <- new.tag.remove
print(best.cutoff.remove)

# no scRef plus
true.tag.remove <- meta.tag.remove$ori.tag
true.tag.remove[true.tag.remove == cell.delete] <- 'unknown'
meta.tag.remove$ori.tag <- true.tag.remove
metrics$f1_score(true.tag.remove, meta.tag.remove$scRef.tag, average = 'weighted')
metrics$f1_score(true.tag.remove, new.tag.remove, average = 'weighted')



