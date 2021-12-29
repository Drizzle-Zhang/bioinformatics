library(stringr)
library(reticulate)
# library(metap)

source('/home/zy/my_git/scRef/main/scRef.v3.R')
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
ref.names <- c("Astrocyte", "Neuron", "Oligodendrocyte precursor cell",
               "Oligodendrocyte", "Myelinating oligodendrocyte",
               "Microglia", "Endothelial cell")
names(exp_ref_mat.origin) <- ref.names
# folder saving results
path.out <- '/home/zy/scRef/try_data/evaluation_del_cell_scRef3/'
if (!file.exists(path.out)) {
    dir.create(path.out)
}

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

# function of evaluation
sub.evaluation <- function(exp_ref_mat.origin, data.filter, label.filter, path.out,
                           df.cell.names, cells.ref, cells.unlabeled, dataset, 
                           default.cutoff = 25) {
    df.plot <- data.frame(stringsAsFactors = F)
    list.metrics <- list()
    if (dataset == 'Tasic') {name.neuron <- 'Neuron'}
    if (dataset == 'Zeisel') {name.neuron <- 'neurons'}
    if (dataset == 'Habib') {name.neuron <- 'Neurons'}
    
    for (i in 1:length(cells.ref)) {
        del.cell.ref <- cells.ref[i]
        del.cell.unlabeled <- cells.unlabeled[i]
        df.sub <- data.frame(stringsAsFactors = F)
        sub.list <- list()
        
        # find mark genes
        exp_ref_mat <- 
            exp_ref_mat.origin[, setdiff(names(exp_ref_mat.origin), c(del.cell.ref))]
        exp_sc_mat = data.filter
        
        ######################### run scRef
        result.scref <- SCREF3(exp_sc_mat, exp_ref_mat, type_ref = 'fpkm', min_cell=5, CPU = num.cpu)
        meta.tag <- merge(result.scref$final.out, label.filter, by = 'row.names')
        row.names(meta.tag) <- meta.tag$Row.names
        meta.tag$Row.names <- NULL
        names(meta.tag) <- c('scRef.tag', 'log10Pval', 'ori.tag')
        
        # confirm label
        ori.tag = meta.tag$ori.tag
        scRef.tag = meta.tag$scRef.tag
        # save txt file
        meta.tag <- meta.tag[order(meta.tag$log10Pval),]
        write.table(meta.tag, 
                    paste0(path.out, 'tags_', dataset,
                           '_scRef_', del.cell.ref, '.txt'), 
                    sep = '\t', quote = F)
        # df.tags <- read.delim(paste0(path.out, 'tags_', dataset,
        #                              '_scRef_', del.cell.ref, '.txt'), row.names = 1)
        # ori.tag <- df.tags$ori.tag
        # scRef.tag <- df.tags$scRef.tag
        
        # ideal condition
        if (del.cell.ref == '') {
            # uniform tags
            for (j in 1:dim(df.cell.names)[1]) {
                scRef.tag[scRef.tag == df.cell.names[j, 'ref.name']] <- 
                    df.cell.names[j, 'sc.name']
            }
            f1.merged <- metrics$f1_score(ori.tag, scRef.tag, average = 'weighted')
            df.sub[1, 1] <- dataset
            df.sub[1, 2] <- 'Neuron merged'
            df.sub[1, 3] <- 'UpLimit'
            df.sub[1, 4] <- f1.merged
            df.sub[1, 5] <- 'no cutoff'
            df.sub[1, 6] <- -1
            df.sub[1, 7] <- 'weighted macro F1'
            df.sub[1, 8] <- 'UpLimit'
            df.sub[1, 9] <- 'macro F1'
            df.sub[1, 10] <- metrics$f1_score(ori.tag, scRef.tag, average = 'macro')
            df.sub[1, 11] <- 'micro F1'
            df.sub[1, 12] <- metrics$f1_score(ori.tag, scRef.tag, average = 'micro')
            
            # remove neuron
            ori.tag.remove <- ori.tag[ori.tag != name.neuron]
            scRef.tag.remove <- scRef.tag[ori.tag != name.neuron]
            f1.removed <- metrics$f1_score(ori.tag.remove, scRef.tag.remove, average = 'weighted')
            df.sub[2, 1] <- dataset
            df.sub[2, 2] <- 'Neuron removed'
            df.sub[2, 3] <- 'UpLimit'
            df.sub[2, 4] <- f1.removed
            df.sub[2, 5] <- 'no cutoff'
            df.sub[2, 6] <- -1
            df.sub[2, 7] <- 'weighted macro F1'
            df.sub[2, 8] <- 'UpLimit'
            df.sub[2, 9] <- 'macro F1'
            df.sub[2, 10] <- metrics$f1_score(ori.tag, scRef.tag, average = 'macro')
            df.sub[2, 11] <- 'micro F1'
            df.sub[2, 12] <- metrics$f1_score(ori.tag, scRef.tag, average = 'micro')
            
            df.plot <- rbind(df.plot, df.sub)
            next()
        }

        # evaluation
        # uniform tags
        for (j in 1:dim(df.cell.names)[1]) {
            scRef.tag[scRef.tag == df.cell.names[j, 'ref.name']] <- 
                df.cell.names[j, 'sc.name']
        }
        meta.tag$scRef.tag <- scRef.tag
        vec.cutoff <- 1:50
        df.metrics <- data.frame(cutoff = vec.cutoff, 
                                 weighted.f1 = rep(0, length(vec.cutoff)), 
                                 macro.f1 = rep(0, length(vec.cutoff)), 
                                 micro.f1 = rep(0, length(vec.cutoff)), 
                                 REC = rep(0, length(vec.cutoff)))
        for (idx in 1:length(vec.cutoff)) {
            cutoff <- vec.cutoff[idx]
            true.tag <- meta.tag$ori.tag
            true.tag[true.tag == del.cell.unlabeled] <- 'unknown'
            our.tag <- meta.tag$scRef.tag
            our.tag[meta.tag$log10Pval <= cutoff] <- 'unknown'
            df.metrics[idx, 'weighted.f1'] <- 
                metrics$f1_score(true.tag, our.tag, average = 'weighted')
            df.metrics[idx, 'macro.f1'] <- 
                metrics$f1_score(true.tag, our.tag, average = 'macro')
            df.metrics[idx, 'micro.f1'] <- 
                metrics$f1_score(true.tag, our.tag, average = 'micro')
            df.metrics[idx, 'accuracy'] <- 
                metrics$accuracy_score(true.tag, our.tag)
            df.metrics[idx, 'REC.f1'] <- 
                REC(true.tag, meta.tag$scRef.tag, our.tag, method = 'f1')
            df.metrics[idx, 'REC.precision'] <- 
                REC(true.tag, meta.tag$scRef.tag, our.tag, method = 'precision')
            df.metrics[idx, 'REC.recall'] <- 
                REC(true.tag, meta.tag$scRef.tag, our.tag, method = 'recall')
            df.metrics[idx, 'REC.accuracy'] <- 
                REC(true.tag, meta.tag$scRef.tag, our.tag, method = 'accuracy')
        }
        sub.list[['merge']] <- df.metrics
        
        # best cutoff
        best.cutoff <- df.metrics$cutoff[df.metrics$weighted.f1 == max(df.metrics$weighted.f1)][1]
        new.tag <- meta.tag$scRef.tag
        new.tag[meta.tag$log10Pval <= best.cutoff] <- 'unknown'
        meta.tag$new.tag <- new.tag
        # print(best.cutoff)
        
        # no scRef plus
        true.tag <- meta.tag$ori.tag
        true.tag[true.tag == del.cell.unlabeled] <- 'unknown'
        meta.tag$ori.tag <- true.tag
        
        df.sub[1, 1] <- dataset
        df.sub[1, 2] <- 'Neuron merged'
        df.sub[1, 3] <- del.cell.ref
        df.sub[1, 4] <- metrics$f1_score(true.tag, scRef.tag, average = 'weighted')
        df.sub[1, 5] <- 'control'
        df.sub[1, 6] <- 0
        df.sub[1, 7] <- 'weighted macro F1'
        df.sub[1, 8] <- 'scRef'
        df.sub[1, 9] <- 'macro F1'
        df.sub[1, 10] <- metrics$f1_score(ori.tag, scRef.tag, average = 'macro')
        df.sub[1, 11] <- 'micro F1'
        df.sub[1, 12] <- metrics$f1_score(ori.tag, scRef.tag, average = 'micro')
        
        df.sub[2, 1] <- dataset
        df.sub[2, 2] <- 'Neuron merged'
        df.sub[2, 3] <- del.cell.ref
        df.sub[2, 4] <- metrics$f1_score(true.tag, new.tag, average = 'weighted')
        df.sub[2, 5] <- 'best cutoff'
        df.sub[2, 6] <- best.cutoff
        df.sub[2, 7] <- 'weighted macro F1'
        df.sub[2, 8] <- 'scRef plus'
        df.sub[2, 9] <- 'macro F1'
        df.sub[2, 10] <- metrics$f1_score(ori.tag, new.tag, average = 'macro')
        df.sub[2, 11] <- 'micro F1'
        df.sub[2, 12] <- metrics$f1_score(ori.tag, new.tag, average = 'micro')
        
        # rate of error correction
        df.sub[3, 1] <- dataset
        df.sub[3, 2] <- 'Neuron merged'
        df.sub[3, 3] <- del.cell.ref
        df.sub[3, 4] <- REC(true.tag, scRef.tag, new.tag)
        df.sub[3, 5] <- 'best cutoff'
        df.sub[3, 6] <- best.cutoff
        df.sub[3, 7] <- 'Rate of error correction'
        df.sub[3, 8] <- 'scRef plus'
        df.sub[3, 9] <- REC(true.tag, scRef.tag, new.tag, method = 'precision')
        df.sub[3, 10] <- 'precision'
        df.sub[3, 11] <- REC(true.tag, scRef.tag, new.tag, method = 'recall')
        df.sub[3, 12] <- 'recall'
        
        # default cutoff
        best.cutoff <- df.metrics$cutoff[df.metrics$weighted.f1 == max(df.metrics$weighted.f1)][1]
        new.tag <- meta.tag$scRef.tag
        new.tag[meta.tag$log10Pval <= default.cutoff] <- 'unknown'

        df.sub[4, 1] <- dataset
        df.sub[4, 2] <- 'Neuron merged'
        df.sub[4, 3] <- del.cell.ref
        df.sub[4, 4] <- metrics$f1_score(true.tag, new.tag, average = 'weighted')
        df.sub[4, 5] <- 'default cutoff'
        df.sub[4, 6] <- default.cutoff
        df.sub[4, 7] <- 'weighted macro F1'
        df.sub[4, 8] <- 'scRef plus'
        df.sub[4, 9] <- 'macro F1'
        df.sub[4, 10] <- metrics$f1_score(ori.tag, new.tag, average = 'macro')
        df.sub[4, 11] <- 'micro F1'
        df.sub[4, 12] <- metrics$f1_score(ori.tag, new.tag, average = 'micro')
        
        # rate of error correction
        df.sub[5, 1] <- dataset
        df.sub[5, 2] <- 'Neuron merged'
        df.sub[5, 3] <- del.cell.ref
        df.sub[5, 4] <- REC(true.tag, scRef.tag, new.tag)
        df.sub[5, 5] <- 'default cutoff'
        df.sub[5, 6] <- default.cutoff
        df.sub[5, 7] <- 'Rate of error correction'
        df.sub[5, 8] <- 'scRef plus'
        df.sub[5, 9] <- REC(true.tag, scRef.tag, new.tag, method = 'precision')
        df.sub[5, 10] <- 'precision'
        df.sub[5, 11] <- REC(true.tag, scRef.tag, new.tag, method = 'recall')
        df.sub[5, 12] <- 'recall'
        
        # save txt file
        write.table(meta.tag, 
                    paste0(path.out, 'tags_', dataset, '_scRef2_', del.cell.ref, '.txt'), 
                    sep = '\t', quote = F)
        
        
        # remove neuron
        meta.tag.remove <- meta.tag[meta.tag$ori.tag != name.neuron, c("ori.tag", "scRef.tag", "log10Pval")]
        vec.cutoff <- 1:50
        df.metrics.remove <- data.frame(cutoff = vec.cutoff, 
                                        weighted.f1 = rep(0, length(vec.cutoff)), 
                                        macro.f1 = rep(0, length(vec.cutoff)), 
                                        micro.f1 = rep(0, length(vec.cutoff)), 
                                        REC = rep(0, length(vec.cutoff)))
        for (idx in 1:length(vec.cutoff)) {
            cutoff <- vec.cutoff[idx]
            true.tag <- meta.tag.remove$ori.tag
            true.tag[true.tag == del.cell.unlabeled] <- 'unknown'
            our.tag <- meta.tag.remove$scRef.tag
            our.tag[meta.tag.remove$log10Pval <= cutoff] <- 'unknown'
            df.metrics.remove[idx, 'weighted.f1'] <- 
                metrics$f1_score(true.tag, our.tag, average = 'weighted')
            df.metrics.remove[idx, 'macro.f1'] <- 
                metrics$f1_score(true.tag, our.tag, average = 'macro')
            df.metrics.remove[idx, 'micro.f1'] <- 
                metrics$f1_score(true.tag, our.tag, average = 'micro')
            df.metrics.remove[idx, 'accuracy'] <- 
                metrics$accuracy_score(true.tag, our.tag)
            df.metrics.remove[idx, 'REC.f1'] <- 
                REC(true.tag, meta.tag.remove$scRef.tag, our.tag, method = 'f1')
            df.metrics.remove[idx, 'REC.precision'] <- 
                REC(true.tag, meta.tag.remove$scRef.tag, our.tag, method = 'precision')
            df.metrics.remove[idx, 'REC.recall'] <- 
                REC(true.tag, meta.tag.remove$scRef.tag, our.tag, method = 'recall')
            df.metrics.remove[idx, 'REC.accuracy'] <- 
                REC(true.tag, meta.tag.remove$scRef.tag, our.tag, method = 'accuracy')
        }
        sub.list[['remove']] <- df.metrics.remove
        # best cutoff
        best.cutoff.remove <- df.metrics.remove$cutoff[df.metrics.remove$weighted.f1 == max(df.metrics.remove$weighted.f1)][1]
        new.tag.remove <- meta.tag.remove$scRef.tag
        new.tag.remove[meta.tag.remove$log10Pval <= best.cutoff.remove] <- 'unknown'
        meta.tag.remove$new.tag <- new.tag.remove
        true.tag.remove <- meta.tag.remove$ori.tag
        true.tag.remove[true.tag.remove == del.cell.unlabeled] <- 'unknown'
        meta.tag.remove$ori.tag <- true.tag.remove
        scRef.tag.remove <- meta.tag.remove$scRef.tag
        
        df.sub[6, 1] <- dataset
        df.sub[6, 2] <- 'Neuron removed'
        df.sub[6, 3] <- del.cell.ref
        df.sub[6, 4] <- metrics$f1_score(true.tag.remove, scRef.tag.remove, average = 'weighted')
        df.sub[6, 5] <- 'control'
        df.sub[6, 6] <- 0
        df.sub[6, 7] <- 'weighted macro F1'
        df.sub[6, 8] <- 'scRef'
        df.sub[6, 9] <- 'macro F1'
        df.sub[6, 10] <- metrics$f1_score(true.tag.remove, scRef.tag.remove, average = 'macro')
        df.sub[6, 11] <- 'micro F1'
        df.sub[6, 12] <- metrics$f1_score(true.tag.remove, scRef.tag.remove, average = 'micro')
        
        df.sub[7, 1] <- dataset
        df.sub[7, 2] <- 'Neuron removed'
        df.sub[7, 3] <- del.cell.ref
        df.sub[7, 4] <- metrics$f1_score(true.tag.remove, new.tag.remove, average = 'weighted')
        df.sub[7, 5] <- 'best cutoff'
        df.sub[7, 6] <- best.cutoff
        df.sub[7, 7] <- 'weighted macro F1'
        df.sub[7, 8] <- 'scRef plus'
        df.sub[7, 9] <- 'macro F1'
        df.sub[7, 10] <- metrics$f1_score(true.tag.remove, new.tag.remove, average = 'macro')
        df.sub[7, 11] <- 'micro F1'
        df.sub[7, 12] <- metrics$f1_score(true.tag.remove, new.tag.remove, average = 'micro')
        
        # rate of error correction
        df.sub[8, 1] <- dataset
        df.sub[8, 2] <- 'Neuron removed'
        df.sub[8, 3] <- del.cell.ref
        df.sub[8, 4] <- REC(true.tag.remove, scRef.tag.remove, new.tag.remove)
        df.sub[8, 5] <- 'best cutoff'
        df.sub[8, 6] <- best.cutoff
        df.sub[8, 7] <- 'Rate of error correction'
        df.sub[8, 8] <- 'scRef plus'
        df.sub[8, 9] <- REC(true.tag.remove, scRef.tag.remove, new.tag.remove, method = 'precision')
        df.sub[8, 10] <- 'precision'
        df.sub[8, 11] <- REC(true.tag.remove, scRef.tag.remove, new.tag.remove, method = 'recall')
        df.sub[8, 12] <- 'recall'
        
        
        # default cutoff
        new.tag.remove <- meta.tag.remove$scRef.tag
        new.tag.remove[meta.tag.remove$log10Pval <= default.cutoff] <- 'unknown'
        true.tag.remove <- meta.tag.remove$ori.tag

        df.sub[9, 1] <- dataset
        df.sub[9, 2] <- 'Neuron removed'
        df.sub[9, 3] <- del.cell.ref
        df.sub[9, 4] <- metrics$f1_score(true.tag.remove, new.tag.remove, average = 'weighted')
        df.sub[9, 5] <- 'default cutoff'
        df.sub[9, 6] <- default.cutoff
        df.sub[9, 7] <- 'weighted macro F1'
        df.sub[9, 8] <- 'scRef plus'
        df.sub[9, 9] <- 'macro F1'
        df.sub[9, 10] <- metrics$f1_score(true.tag.remove, new.tag.remove, average = 'macro')
        df.sub[9, 11] <- 'micro F1'
        df.sub[9, 12] <- metrics$f1_score(true.tag.remove, new.tag.remove, average = 'micro')
        
        # rate of error correction
        df.sub[10, 1] <- dataset
        df.sub[10, 2] <- 'Neuron removed'
        df.sub[10, 3] <- del.cell.ref
        df.sub[10, 4] <- REC(true.tag.remove, scRef.tag.remove, new.tag.remove)
        df.sub[10, 5] <- 'default cutoff'
        df.sub[10, 6] <- default.cutoff
        df.sub[10, 7] <- 'Rate of error correction'
        df.sub[10, 8] <- 'scRef plus'
        df.sub[10, 9] <- REC(true.tag.remove, scRef.tag.remove, new.tag.remove, method = 'precision')
        df.sub[10, 10] <- 'precision'
        df.sub[10, 11] <- REC(true.tag.remove, scRef.tag.remove, new.tag.remove, method = 'recall')
        df.sub[10, 12] <- 'recall'
        
        df.plot <- rbind(df.plot, df.sub)
        list.metrics[[del.cell.ref]] <- sub.list
        
    }
    
    OUT <- list()
    OUT[['plot']] <- df.plot
    OUT[['metrics']] <- list.metrics
    saveRDS(OUT, file = paste0(path.out, 'plot_', dataset, '_scRef2.Rdata'))
    return(OUT)
    
}

######################################################################################
######################### Zeisel
dataset <- 'Zeisel'
file.data.unlabeled <- paste0('./summary/', dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0('./summary/', dataset, '_exp_sc_mat_cluster_merged.txt')
OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled)
data.filter <- OUT$data.filter
label.filter <- OUT$label.filter
# list of cell names
df.cell.names <- data.frame(
    ref.name = ref.names, 
    sc.name = c("astrocytes_ependymal", "neurons", "Oligodendrocyte precursor cell",
                "oligodendrocytes", "oligodendrocytes",
                "microglia", "endothelial-mural"))

# delete specific cell
cells.ref <- c("", "Astrocyte", "Microglia", "Endothelial cell")
cells.unlabeled <- c("", "astrocytes_ependymal", "microglia", "endothelial-mural")
# evaluation
out <- sub.evaluation(exp_ref_mat.origin, data.filter, label.filter, path.out,
                      df.cell.names, cells.ref, cells.unlabeled, dataset)

######################### Tasic
dataset <- 'Tasic'
file.data.unlabeled <- paste0('./summary/', dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0('./summary/', dataset, '_exp_sc_mat_cluster_merged.txt')
OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled)
data.filter <- OUT$data.filter
label.filter <- OUT$label.filter
# list of cell names
df.cell.names <- data.frame(
    ref.name = ref.names, 
    sc.name = c("Astrocyte", "Neuron", "Oligodendrocyte Precursor Cell",
                "Oligodendrocyte", "Oligodendrocyte",
                "Microglia", "Endothelial Cell"))

# delete specific cell
cells.ref <- c("", "Astrocyte", "Microglia", "Endothelial cell", "Oligodendrocyte precursor cell")
cells.unlabeled <- c("", "Astrocyte", "Microglia", "Endothelial Cell", "Oligodendrocyte Precursor Cell")
# evaluation
out <- sub.evaluation(exp_ref_mat.origin, data.filter, label.filter, path.out,
                      df.cell.names, cells.ref, cells.unlabeled, dataset)


######################### Habib
dataset <- 'Habib'
file.data.unlabeled <- paste0('./summary/', dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0('./summary/', dataset, '_exp_sc_mat_cluster_merged.txt')
OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled)
data.filter <- OUT$data.filter
label.filter <- OUT$label.filter
# list of cell names
df.cell.names <- data.frame(
    ref.name = ref.names, 
    sc.name = c("Astrocyte", "Neurons", "OPC",
                "Oligodend", "Oligodend",
                "microglia", "EndothelialCells"))

# delete specific cell
cells.ref <- c("", "Astrocyte", "Microglia", "Endothelial cell", "Oligodendrocyte precursor cell")
cells.unlabeled <- c("", "Astrocyte", "microglia", "EndothelialCells", "OPC")
# evaluation
out <- sub.evaluation(exp_ref_mat.origin, data.filter, label.filter, path.out,
                      df.cell.names, cells.ref, cells.unlabeled, dataset)



