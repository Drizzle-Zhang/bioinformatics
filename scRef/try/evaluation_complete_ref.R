library(stringr)
library(reticulate)

# sample classifier
sample.classifier <- function(exp_sc_mat, list.cell.genes, method.test = 'wilcox', iter.permutation = NULL) {
    library(foreach)
    library(doParallel)
    library(coin)
    # confirm label
    registerDoParallel(10)
    confirm.classify <- function(exp_sc_mat, list.cell.genes, method.test, barcode) {
        expression.barcode <- exp_sc_mat[, barcode]
        bool.mark.gene <- rep(1, dim(exp_sc_mat)[1])
        vector.pval <- c()
        cells <- names(list.cell.genes)
        for (cell in cells) {
            genes.marker <- list.cell.genes[[cell]]
            bool.mark.gene[row.names(exp_sc_mat) %in% genes.marker] <- 2
            test.in <- cbind(expression.barcode, bool.mark.gene)
            test.in <- as.data.frame(test.in)
            names(test.in) <- c('expression.level', 'factor.mark.gene')
            test.in$factor.mark.gene <- as.factor(test.in$factor.mark.gene)
            if (method.test == 'wilcox') {
                out.test <- wilcox.test(expression.level ~ factor.mark.gene, data = test.in, 
                                        alternative = 'less')
                pvalue <- out.test$p.value
            }
            if (method.test == 't.test') {
                out.test <- t.test(expression.level ~ factor.mark.gene, data = test.in, 
                                   alternative = 'less', var.equal = T)
                pvalue <- out.test$p.value
            }
            if (method.test == 'oneway_test') {
                if (is.null(iter.permutation)) {iter.permutation = 1000}
                out.test <- oneway_test(formula = expression.level ~ factor.mark.gene, data = test.in, 
                                        alternative = 'less', distribution = approximate(B=iter.permutation))
                pvalue <- pvalue(out.test)
            }
            vector.pval <- c(vector.pval, pvalue)
        }
        min.pval <- min(vector.pval)
        tag.cell <- cells[vector.pval == min.pval]
        
        return(data.frame(tag.test = tag.cell, pvalue = min.pval))
        
    }
    
    out.par <- foreach(barcode = dimnames(exp_sc_mat)[[2]], .combine = rbind) %dopar% confirm.classify(exp_sc_mat, list.cell.genes, method.test, barcode)
    meta.tag <- out.par
    row.names(meta.tag) <- dimnames(exp_sc_mat)[[2]]
    # meta.tag$qvalue <- p.adjust(meta.tag$pvalue, method = 'BH')
    
    return(meta.tag)
    
}


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

out.markers <- find.markers(exp_ref_mat.origin)
list.cell.genes <- out.markers[['list.cell.genes']]
genes.ref <- dimnames(out.markers[['exp_ref_mat']])[[1]]

df.evaluate <- data.frame()

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
exp_ref_mat <- exp_ref_mat.origin[genes.ref,]
exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
gene_sc=rownames(exp_sc_mat)
gene_ref=rownames(exp_ref_mat)
gene_over= gene_sc[which(gene_sc %in% gene_ref)]
exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
print('Number of overlapped genes:')
print(nrow(exp_sc_mat))

################### simple classifier
res.simple <- sample.classifier(exp_sc_mat, list.cell.genes)

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
meta.tag <- cbind(meta.tag, res.simple)

# evaluation
# import python package: sklearn.metrics
use_condaenv("/home/zy/tools/anaconda3")
# py_config()
metrics <- import('sklearn.metrics')
# uniform tags
tag.test <- meta.tag$tag.test
tag.test[tag.test == "Astrocyte"] <- "astrocytes_ependymal"
tag.test[tag.test == "Newly Formed Oligodendrocyte"] <- "oligodendrocytes"
tag.test[tag.test == "Myelinating oligodendrocyte"] <- "oligodendrocytes"
tag.test[tag.test == "Endothelial cell"] <- "endothelial-mural" 
tag.test[tag.test == "Neuron"] <- "neurons" 
tag.test[tag.test == "Microglia"] <- "microglia" 
meta.tag$tag.test <- tag.test

scRef.tag[scRef.tag == "Astrocyte"] <- "astrocytes_ependymal"
scRef.tag[scRef.tag == "Newly Formed Oligodendrocyte"] <- "oligodendrocytes"
scRef.tag[scRef.tag == "Myelinating oligodendrocyte"] <- "oligodendrocytes"
scRef.tag[scRef.tag == "Endothelial cell"] <- "endothelial-mural" 
scRef.tag[scRef.tag == "Neuron"] <- "neurons" 
scRef.tag[scRef.tag == "Microglia"] <- "microglia" 
meta.tag$scRef.tag <- scRef.tag

# f1 score
f1.scref <- metrics$f1_score(ori.tag, scRef.tag, average = 'weighted')
f1.test <- metrics$f1_score(ori.tag, tag.test, average = 'weighted')

# remove neuron
meta.tag.remove <- meta.tag[meta.tag$ori.tag != "neurons", ]
true.tag.remove <- meta.tag.remove$ori.tag
f1.scref.remove <- metrics$f1_score(true.tag.remove, meta.tag.remove$scRef.tag, average = 'weighted')
f1.test.remove <- metrics$f1_score(true.tag.remove, meta.tag.remove$tag.test, average = 'weighted')

df.evaluate <- 
    rbind(df.evaluate,
          data.frame(test.method = c(rep('Neuron merged', 2), rep('Neuron removed', 2)),
                     dataset = rep('Zeisel', 4),
                     method = c('scRef', 'Simple test', 'scRef', 'Simple test'),
                     macro.f1 = c(f1.scref, f1.test, f1.scref.remove, f1.test.remove)))

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

################### simple classifier
res.simple <- sample.classifier(exp_sc_mat, list.cell.genes)

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
meta.tag <- cbind(meta.tag, res.simple)

# evaluation
# import python package: sklearn.metrics
use_condaenv("/home/zy/tools/anaconda3")
# py_config()
metrics <- import('sklearn.metrics')
# uniform tags
tag.test <- meta.tag$tag.test
tag.test[tag.test == "Newly Formed Oligodendrocyte"] <- "Oligodendrocyte"
tag.test[tag.test == "Myelinating oligodendrocyte"] <- "Oligodendrocyte"
tag.test[tag.test == "Endothelial cell"] <- "Endothelial Cell" 
tag.test[tag.test == "Oligodendrocyte precursor cell"] <- "Oligodendrocyte Precursor Cell" 
meta.tag$tag.test <- tag.test

scRef.tag[scRef.tag == "Newly Formed Oligodendrocyte"] <- "Oligodendrocyte"
scRef.tag[scRef.tag == "Myelinating oligodendrocyte"] <- "Oligodendrocyte"
scRef.tag[scRef.tag == "Endothelial cell"] <- "Endothelial Cell" 
scRef.tag[scRef.tag == "Oligodendrocyte precursor cell"] <- "Oligodendrocyte Precursor Cell" 
meta.tag$scRef.tag <- scRef.tag

# f1 score
f1.scref <- metrics$f1_score(ori.tag, scRef.tag, average = 'weighted')
f1.test <- metrics$f1_score(ori.tag, tag.test, average = 'weighted')

# remove neuron
meta.tag.remove <- meta.tag[meta.tag$ori.tag != "Neuron", ]
true.tag.remove <- meta.tag.remove$ori.tag
f1.scref.remove <- metrics$f1_score(true.tag.remove, meta.tag.remove$scRef.tag, average = 'weighted')
f1.test.remove <- metrics$f1_score(true.tag.remove, meta.tag.remove$tag.test, average = 'weighted')

df.evaluate <- 
    rbind(df.evaluate,
          data.frame(test.method = c(rep('Neuron merged', 2), rep('Neuron removed', 2)),
                     dataset = rep('Tasic', 4),
                     method = c('scRef', 'Simple test', 'scRef', 'Simple test'),
                     macro.f1 = c(f1.scref, f1.test, f1.scref.remove, f1.test.remove)))

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

################### simple classifier
res.simple <- sample.classifier(exp_sc_mat, list.cell.genes)

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
meta.tag <- cbind(meta.tag, res.simple)

# evaluation
# import python package: sklearn.metrics
use_condaenv("/home/zy/tools/anaconda3")
# py_config()
metrics <- import('sklearn.metrics')
# uniform tags
tag.test <- meta.tag$tag.test
tag.test[tag.test == "Oligodendrocyte precursor cell"] <- "OPC" 
tag.test[tag.test == "Newly Formed Oligodendrocyte"] <- "Oligodend"
tag.test[tag.test == "Myelinating oligodendrocyte"] <- "Oligodend"
tag.test[tag.test == "Endothelial cell"] <- "EndothelialCells" 
tag.test[tag.test == "Neuron"] <- "Neurons" 
tag.test[tag.test == "Microglia"] <- "microglia" 
meta.tag$tag.test <- tag.test

scRef.tag[scRef.tag == "Oligodendrocyte precursor cell"] <- "OPC" 
scRef.tag[scRef.tag == "Newly Formed Oligodendrocyte"] <- "Oligodend"
scRef.tag[scRef.tag == "Myelinating oligodendrocyte"] <- "Oligodend"
scRef.tag[scRef.tag == "Endothelial cell"] <- "EndothelialCells" 
scRef.tag[scRef.tag == "Neuron"] <- "Neurons" 
scRef.tag[scRef.tag == "Microglia"] <- "microglia" 
meta.tag$scRef.tag <- scRef.tag

# f1 score
f1.scref <- metrics$f1_score(ori.tag, scRef.tag, average = 'weighted')
f1.test <- metrics$f1_score(ori.tag, tag.test, average = 'weighted')

# remove neuron
meta.tag.remove <- meta.tag[meta.tag$ori.tag != "Neurons", ]
true.tag.remove <- meta.tag.remove$ori.tag
f1.scref.remove <- metrics$f1_score(true.tag.remove, meta.tag.remove$scRef.tag, average = 'weighted')
f1.test.remove <- metrics$f1_score(true.tag.remove, meta.tag.remove$tag.test, average = 'weighted')

df.evaluate <- 
    rbind(df.evaluate,
          data.frame(test.method = c(rep('Neuron merged', 2), rep('Neuron removed', 2)),
                     dataset = rep('Habib', 4),
                     method = c('scRef', 'Simple test', 'scRef', 'Simple test'),
                     macro.f1 = c(f1.scref, f1.test, f1.scref.remove, f1.test.remove)))

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

################### simple classifier
res.simple <- sample.classifier(exp_sc_mat, list.cell.genes)


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
meta.tag <- cbind(meta.tag, res.simple)

# evaluation
# import python package: sklearn.metrics
use_condaenv("/home/zy/tools/anaconda3")
# py_config()
metrics <- import('sklearn.metrics')
# uniform tags
tag.test <- meta.tag$tag.test
tag.test[tag.test == "Oligodendrocyte precursor cell"] <- "Oligodendrocyte.Precursor.Cell"
tag.test[tag.test == "Newly Formed Oligodendrocyte"] <- "Oligodendrocyte"
tag.test[tag.test == "Myelinating oligodendrocyte"] <- "Oligodendrocyte"
tag.test[tag.test == "Endothelial cell"] <- "Endothelial.Cells" 
tag.test[tag.test == "Neuron"] <- "Neuron" 
tag.test[tag.test == "Microglia"] <- "Microglia" 
tag.test[tag.test == "Astrocyte"] <- "Astrocytes" 
meta.tag$tag.test <- tag.test

scRef.tag[scRef.tag == "Oligodendrocyte precursor cell"] <- "Oligodendrocyte.Precursor.Cell"
scRef.tag[scRef.tag == "Newly Formed Oligodendrocyte"] <- "Oligodendrocyte"
scRef.tag[scRef.tag == "Myelinating oligodendrocyte"] <- "Oligodendrocyte"
scRef.tag[scRef.tag == "Endothelial cell"] <- "Endothelial.Cells" 
scRef.tag[scRef.tag == "Neuron"] <- "Neuron" 
scRef.tag[scRef.tag == "Microglia"] <- "Microglia" 
scRef.tag[scRef.tag == "Astrocyte"] <- "Astrocytes" 
meta.tag$scRef.tag <- scRef.tag

# f1 score
f1.scref <- metrics$f1_score(ori.tag, scRef.tag, average = 'weighted')
f1.test <- metrics$f1_score(ori.tag, tag.test, average = 'weighted')

# remove neuron
meta.tag.remove <- meta.tag[meta.tag$ori.tag != "Neuron", ]
true.tag.remove <- meta.tag.remove$ori.tag
f1.scref.remove <- metrics$f1_score(true.tag.remove, meta.tag.remove$scRef.tag, average = 'weighted')
f1.test.remove <- metrics$f1_score(true.tag.remove, meta.tag.remove$tag.test, average = 'weighted')

df.evaluate <- 
    rbind(df.evaluate,
          data.frame(test.method = c(rep('Neuron merged', 2), rep('Neuron removed', 2)),
                     dataset = rep('Hochgerner', 4),
                     method = c('scRef', 'Simple test', 'scRef', 'Simple test'),
                     macro.f1 = c(f1.scref, f1.test, f1.scref.remove, f1.test.remove)))

library(ggplot2)
ggplot(df.evaluate,
       aes(x = dataset, y = macro.f1, fill = method)) +
    geom_bar(position = 'dodge', stat = 'identity') +
    facet_wrap(~ test.method, scales = 'free', ncol = 2) +
    labs(title = "", y = 'Weighted Macro F1', x = '', fill = 'Method') + 
    coord_cartesian(ylim = c(0.5, 1)) + 
    theme(axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1))   


