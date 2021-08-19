
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

path.input <- '/mdshare/node9/zy/scRef/sc_data/'
path.output <- '/mdshare/node9/zy/scRef/atlas_anno/'
dataset <- 'Tasic2018'
file.data.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat_cluster_original.txt')
# OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, del.label = c('miss'))
# saveRDS(OUT, file = paste0(path.output, dataset, '.Rdata'))
OUT <- readRDS(paste0('/mdshare/node9/zy/scRef/Benchmark/mouse_brain/', dataset, '.Rdata'))
exp_sc_mat <- OUT$mat_exp
label_sc <- OUT$label

library(Seurat)
library(scMAGIC)
data("MCA_ref")
set.seed(123)
sel.sample <- rownames(label_sc)[label_sc[,1] %in% 
                                   c('Endothelial_Endo', 'Endothelial_Peri', 'Endothelial_SMC',
                                     'Non-Neuronal_Astro', 'Non-Neuronal_Macrophage',
                                     'Non-Neuronal_Oligo', 'Non-Neuronal_VLMC')]
col.neu <- setdiff(colnames(exp_sc_mat), sel.sample)
sel.sample <- c(sel.sample, sample(col.neu,3500-length(sel.sample)))
exp_input <- exp_sc_mat[, sel.sample]
sel_label <- label_sc[sel.sample, ]

# sel.col <- sample(colnames(exp_sc_mat),3000)
# exp_input <- exp_sc_mat[, sel.col]
# sel_label <- label_sc[sel.col, ]
time1 <- Sys.time()
output.scMAGIC <- scMAGIC(exp_input, MCA_ref,
                          type_ref = 'sum-counts', use_RUVseq = F,
                          min_cell = 5, num_threads = 10)
time2 <- Sys.time()
time.diff <- difftime(time2, time1, units = 'mins')
pred.scMAGIC <- output.scMAGIC$scMAGIC.tag
true.tags <- sel_label
table(true.tags, pred.scMAGIC)
# df.tags <- result.scref$combine.out
# df.view <- merge(label_sc, df.tags, by = 'row.names')
# View(df.view)

library(ggplot2)
path.res <- '/mdshare/node9/zy/scRef/figure/atlas_anno/'
method <- 'scMAGIC'
file.pred <- paste0(path.res, 'MCA_', dataset, '_scMAGIC.Rdata')
# saveRDS(pred.scMAGIC, file.pred)
pred.scMAGIC <- readRDS(file.pred)
true.tags <- sel_label
true.tags[true.tags %in% c('Endothelial_Endo')] <- 'Endothelial cell'
true.tags[true.tags %in% c('Endothelial_Peri')] <- 'Pericyte'
true.tags[true.tags %in% c('Endothelial_SMC')] <- 'Smooth muscle cell'
true.tags[true.tags %in% c('GABAergic_', 'GABAergic_Lamp5', 'GABAergic_Meis2', 
                           'GABAergic_Pvalb', 'GABAergic_Serpinf1', 'GABAergic_Sncg', 
                           'GABAergic_Sst', 'GABAergic_Vip')] <- 'GABAergic Neuron'
true.tags[true.tags %in% c('Glutamatergic_', 'Glutamatergic_CR', 'Glutamatergic_L2/3 IT', 
                           'Glutamatergic_L4', 'Glutamatergic_L5 IT', "Glutamatergic_L5 PT", 
                           'Glutamatergic_L6 CT',  'Glutamatergic_L6 IT', 
                           'Glutamatergic_L6b', 'Glutamatergic_NP')] <- 'Glutamatergic Neuron'
true.tags[true.tags %in% c('Non-Neuronal_Astro')] <- 'Astrocyte'
true.tags[true.tags %in% c('Non-Neuronal_Macrophage')] <- 'PVM & Microglia'
true.tags[true.tags %in% c('Non-Neuronal_Oligo')] <- 'Oligodendrocyte & OPC'
true.tags[true.tags %in% c('Non-Neuronal_VLMC')] <- 'VLMC'
col.names <- colnames(exp_input)
sel.cols <- c(setdiff(col.names, col.names[pred.scMAGIC=='Pyramidal neuron cell(Fetal_Brain)']),
              sample(col.names[pred.scMAGIC=='Pyramidal neuron cell(Fetal_Brain)'], 250))
true.tags <- true.tags[col.names%in%sel.cols]
pred.scMAGIC <- pred.scMAGIC[col.names%in%sel.cols]
table(true.tags, pred.scMAGIC)

# simplify colnames of atlas
df.simple <- data.frame(check.names = F)
for (col in colnames(MCA_ref)) {
    col.split <- strsplit(col, split = "(", fixed = T)[[1]]
    col.simple <- col.split[1]
    # col.simple <- paste(col.split[1:(length(col.split)-1)], collapse = '-')
    df.simple <- rbind(df.simple, data.frame(col = col,
                                             col.simple = col.simple))
}
row.names(df.simple) <- df.simple$col

# heatmap
mytable <- table(true.tags, pred.scMAGIC)
mytable.tag <- data.frame(colnames(mytable), df.simple[colnames(mytable), 'col.simple'])
mytable.tag[mytable.tag[,1] == 'Unassigned', 2] <- 'Unassigned'
mytable.tag[,2] <- gsub('Astrocyte_Atp1b2 high', 'Astrocyte', mytable.tag[,2])
mytable.tag[,2] <- gsub('Endothelial cell_Ly6c1 high', 'Endothelial cell', mytable.tag[,2])
mytable.tag[,2] <- gsub('Endothelial cell_Tm4sf1 high', 'Endothelial cell', mytable.tag[,2])
mytable.tag[,2] <- gsub('Endothelial cells_Vwf high', 'Endothelial cell', mytable.tag[,2])
mytable.tag[,2] <- gsub('Ovarian vascular surface endothelium cell', 'Endothelial cell', mytable.tag[,2])
mytable.tag[,2] <- gsub('Vascular endothelial cell', 'Endothelial cell', mytable.tag[,2])
mytable.tag[,2] <- gsub('Macrophage_C1qc high', 'Macrophage', mytable.tag[,2])
mytable.tag[,2] <- gsub('Macrophage_Pf4high', 'Endothelial cell', mytable.tag[,2])
mytable.tag[,2] <- gsub('Atrial cardiomyocyte_Acta2 high', 'Smooth muscle cell', mytable.tag[,2])
mytable.tag[,2] <- gsub('Muscle cell_Pcp4 high', 'Smooth muscle cell', mytable.tag[,2])
mytable.tag[,2] <- gsub('Smooth muscle cell_Acta2 high', 'Smooth muscle cell', mytable.tag[,2])
mytable.tag[,2] <- gsub('Smooth muscle cell_Rgs5 high', 'Smooth muscle cell', mytable.tag[,2])
mytable.tag[,2] <- gsub('T cell_Gm14303 high', 'Other cells in MCA', mytable.tag[,2])
# mytable.tag[,2] <- gsub('Oligodendrocyte precursor cell', 'OPC', mytable.tag[,2])
# mytable.tag[,2] <- gsub('Smooth muscle cell', 'SMC', mytable.tag[,2])
mytable.tag[,2] <- gsub('Neuron', 'Pyramidal neuron cell', mytable.tag[,2])

mytable.sum <- generate_ref(mytable, mytable.tag)
mydata <- data.frame(stringsAsFactors = F)
table.true <- table(true.tags)
for (label1 in rownames(mytable.sum)) {
    row.sum <- table.true[label1]
    for (label2 in unique(c(colnames(mytable.sum), 'Other cells in MCA'))) {
        if (label2 %in% colnames(mytable.sum)) {
            mydata <- rbind(mydata, data.frame(origin = label1, annotation = label2,
                                               count = mytable.sum[label1, label2],
                                               prop = mytable.sum[label1, label2]/row.sum))
        } else {
            mydata <- rbind(mydata, data.frame(origin = label1, annotation = label2,
                                               count = 0, prop = 0))
        }
    }
}

tag_order <- c("Astrocyte", "Endothelial cell",
               "Dopaminergic neurons ",
               "Hypothalamic ependymal cell", "Macrophage",
               "Microglia", "Smooth muscle cell", "Oligodendrocyte precursor cell", 
               "Unassigned", "Pyramidal neuron cell", "Other cells in MCA")
mydata$annotation <- factor(mydata$annotation, levels = tag_order)
mydata$origin <- factor(mydata$origin,
                        levels = c("Astrocyte", "Endothelial cell",
                                   "GABAergic Neuron", "Glutamatergic Neuron",
                                   "PVM & Microglia", "Smooth muscle cell", "Oligodendrocyte & OPC",
                                   "Pericyte", "VLMC"))

plot.heatmap <-
    ggplot(data = mydata, aes(x = origin, y = annotation)) +
    geom_tile(aes(fill = prop)) +
    # scale_fill_gradient2(low = "#C0C0C0", high = "#FFFF00", mid = "#32CD32", midpoint = 0.5) +
    scale_fill_gradient2(low = "#FFF5EE", mid = '#EE7700', high = "#B22222", midpoint = 0.5) +
    labs(fill = 'Proportion', title = 'MCA -> Mouse Neocortex') +
    theme_bw() +
    theme(
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        title = element_text(size = 10, color = "black", family = 'Arial'),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 9, color = "black", family = 'Arial',
                                   angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 9, color = "black", family = 'Arial'),
        legend.title = element_text(
            size = 7, color = "black", family = 'Arial'),
        legend.text = element_text(
            size = 6, color = "black", family = 'Arial'),
        legend.key.size = unit(0.3, 'cm'),
        legend.position = 'none'
    )
ggsave(filename = paste0('heatmap_MCA_', dataset, '_', method, '1.png'),
       path = path.res, plot = plot.heatmap,
       units = 'cm', height = 8.6, width = 9.5)


# accessment
# import python package: sklearn.metrics
library(reticulate)
use_python('/local/zy/tools/anaconda3/bin/python3', required = T)
# py_config()
py_module_available('sklearn')
metrics <- import('sklearn.metrics')

mytable <- table(true.tags, pred.scMAGIC)

ref.names <- colnames(mytable)
all.cell <- names(table(true.tags))
uniform.names <-
    c('Astrocytes', 'Neurons', 'Endothelial cells', 'Endothelial cells',
      'PVMs & Microglia', 'PVMs & Microglia', 'OPC', 'Pyramidal neuron cell',
      'SMC', 'Unassigned', 'Endothelial cells')
df.ref.names <- data.frame(ref.name = ref.names, name = uniform.names)
uniform.names <-
    c("Astrocytes", "Endothelial cells", "Neurons", "Neurons",
      "OPC", "Unassigned", "PVMs & Microglia", "SMC", "Unassigned")
df.sc.names <- data.frame(sc.name = all.cell, name = uniform.names)

simple.evaluation <- function(true.tag, scRef.tag, df.ref.names, df.sc.names) {
    # uniform tags
    for (j in 1:dim(df.ref.names)[1]) {
        scRef.tag[scRef.tag == df.ref.names[j, 'ref.name']] <- df.ref.names[j, 'name']
    }
    scRef.tag[!(scRef.tag %in% df.ref.names$name)] <- 'Unassigned'
    for (j in 1:dim(df.sc.names)[1]) {
        true.tag[true.tag == df.sc.names[j, 'sc.name']] <- df.sc.names[j, 'name']
    }

    percent.unassigned <- sum(scRef.tag == 'Unassigned')/sum(true.tag == 'Unassigned')
    
    # function of balanced accuracy
    # balanced_accuracy <- function(y_true, y_pred, labels) {
    #     mytable <- table(y_true, y_pred)
    #     acc <- c()
    #     for (label in labels) {
    #         num.true.tag <- sum(mytable[label,])
    #         if (label %in% colnames(mytable)) {
    #             num.match <- mytable[label, label]
    #         } else {
    #             num.match <- 0
    #         }
    #         sub.acc <- num.match/num.true.tag
    #         acc <- c(acc, sub.acc)
    #     }
    #     names(acc) <- labels
    #     print(acc)
    #     return(mean(acc))
    # }


    # true.labels <- setdiff(unique(true.tag), 'Unassigned')
    true.labels <- unique(true.tag)
    our.tag <- scRef.tag
    weighted_macro_f1 <- metrics$f1_score(true.tag, our.tag, average = 'weighted', labels = true.labels)
    macro_f1 <- metrics$f1_score(true.tag, our.tag, average = 'macro', labels = true.labels)
    accuracy <- metrics$accuracy_score(true.tag, our.tag)
    # balanced_acc <- balanced_accuracy(true.tag, our.tag, true.labels)
    balanced_acc <- metrics$balanced_accuracy_score(true.tag, our.tag)
    # rm unassigned in tags
    true.tag.rm <- true.tag[our.tag != 'Unassigned']
    our.tag.rm <- our.tag[our.tag != 'Unassigned']
    accuracy.rm.unassigned <- metrics$accuracy_score(true.tag.rm, our.tag.rm)
    macro_f1.rm.unassigned <- metrics$f1_score(true.tag.rm, our.tag.rm, average = 'macro', labels = unique(true.tag.rm))
    ref.labels <- intersect(unique(our.tag.rm), true.labels)
    # balanced_acc <- balanced_accuracy(our.tag.rm, true.tag.rm, ref.labels)
    balanced.accuracy.rm.unassigned <- metrics$balanced_accuracy_score(our.tag.rm, true.tag.rm)

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

    our.labels <- setdiff(unique(our.tag), 'Unassigned')
    precision <- c()
    for (label in our.labels) {
        tmp.true.tag <- true.tag.rm
        tmp.our.tag <- our.tag.rm
        tmp.true.tag[tmp.true.tag != label] <- '0'
        tmp.our.tag[tmp.our.tag != label] <- '0'
        sub.precision <- metrics$precision_score(tmp.true.tag, tmp.our.tag, average = 'binary', pos_label = label)
        precision <- c(precision, sub.precision)

    }
    names(precision) <- our.labels
    mean.precision <- mean(precision)
    
    # remove unassigend in true tags
    true.tag.rm <- true.tag[true.tag != 'Unassigned']
    our.tag.rm <- our.tag[true.tag != 'Unassigned']
    accuracy.rm.newcell <- metrics$accuracy_score(true.tag.rm, our.tag.rm)
    
    out <- list()
    out$percent.unassigned <- percent.unassigned
    out$weighted_macro_f1 <- weighted_macro_f1
    out$macro_f1 <- macro_f1
    out$accuracy <- accuracy
    out$balanced.accuracy <- balanced_acc
    out$f1 <- f1
    out$accuracy.rm.unassigned <- accuracy.rm.unassigned
    out$macro_f1.rm.unassigned <- macro_f1.rm.unassigned
    out$precision.rm.unassigned <- precision
    out$balanced.accuracy.rm.unassigned <- balanced.accuracy.rm.unassigned
    out$mean.precision.rm.unassigned <- mean.precision
    out$accuracy.rm.newcell <- accuracy.rm.newcell
    out$conf <- table(true.tag, our.tag)

    return(out)

}

res.scMAGIC <- simple.evaluation(true.tags, pred.scMAGIC, df.ref.names, df.sc.names)

file.res.scMAGIC <- paste0(path.res, 'RES_MCA_', dataset, '_scMAGIC.Rdata')
# write.table(res.scMAGIC, file.res.scMAGIC)
saveRDS(res.scMAGIC, file.res.scMAGIC)
readRDS(file.res.scMAGIC)

# $accuracy
# [1] 0.9464175
#
# $balanced.accuracy
# [1] 0.8587342

# $accuracy.rm.unassigned
# [1] 0.9697406

# MCA_ref <- MCA_ref[,!(colnames(MCA_ref) %in% c('Schwann cell(Fetal_Brain)',
#                                                'Astroglial cell(Bergman glia)(Brain)',
#                                                'Macrophage_Klf2 high(Brain)'))]
