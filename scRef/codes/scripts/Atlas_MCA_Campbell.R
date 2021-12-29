
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
dataset <- 'Campbell'
file.data.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat_cluster_original.txt')
# OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, del.label = c('miss'))
# saveRDS(OUT, file = paste0(path.output, dataset, '.Rdata'))
OUT <- readRDS(paste0('/mdshare/node9/zy/scRef/Benchmark/mouse_brain/', dataset, '.Rdata'))
exp_Habib <- OUT$mat_exp
label_Habib <- OUT$label
exp_sc_mat <- exp_Habib
label_sc <- label_Habib

library(Seurat)
library(scMAGIC)
data("MCA_ref")
output.scMAGIC <- scMAGIC(exp_sc_mat, MCA_ref,
                          type_ref = 'sum-counts', use_RUVseq = F,
                          num_threads = 10)
pred.scMAGIC <- output.scMAGIC$scMAGIC.tag
true.tags <- label_sc[,1]
table(true.tags, pred.scMAGIC)
# df.tags <- result.scref$combine.out
# df.view <- merge(label_sc, df.tags, by = 'row.names')
# View(df.view)

library(ggplot2)
path.res <- '/mdshare/node9/zy/scRef/figure/atlas_anno/'
method <- 'scMAGIC'
file.pred <- paste0(path.res, 'MCA_', dataset, '_scMAGIC.Rdata')
saveRDS(pred.scMAGIC, file.pred)
pred.scMAGIC <- readRDS(file.pred)
true.tags <- label_sc[,1]
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
mytable.sum <- generate_ref(mytable, mytable.tag)
mydata <- data.frame(stringsAsFactors = F)
table.true <- table(true.tags)
for (label1 in rownames(mytable.sum)) {
    row.sum <- table.true[label1]
    for (label2 in c(unique(df.simple$col.simple), "Unassigned")) {
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

ref.cells <- setdiff(colnames(mytable.sum), c("Ganglion cell_Cartpt high",
                                              "Pan-GABAergic"))
set.seed(1234)
all.cells <- sort(unique(c(sample(setdiff(unique(df.simple$col.simple), ref.cells), 100),
                           ref.cells, "Unassigned")))
pos <- c(-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0)
pos.cells <- c()
for (idx in 1:length(ref.cells)) {
    pos.idx <- pos[idx]
    pos.cells <- c(pos.cells, all.cells[which(all.cells==ref.cells[idx])+pos.idx])
}
mydata <- mydata[mydata$annotation %in% all.cells,]
mydata$annotation <- factor(mydata$annotation, levels = all.cells)
mydata$origin <- factor(mydata$origin,
                        levels = c("Astrocytes", "Oligodendrocytes",
                                   "Endothelial cells",
                                   "Neurons", "Ependymocytes",
                                   "PVMs & Microglia", "OPC",
                                   "Mural cells", "Pars tuberalis",
                                   "Tanycytes", "VLMCs"))

plot.heatmap <-
    ggplot(data = mydata, aes(x = origin, y = annotation)) +
    geom_tile(aes(fill = prop)) +
    # scale_fill_gradient2(low = "#C0C0C0", high = "#FFFF00", mid = "#32CD32", midpoint = 0.5) +
    scale_fill_gradient2(low = "#FFF5EE", mid = '#EE7700', high = "#B22222", midpoint = 0.5) +
    labs(fill = 'Proportion') +
    theme_bw() +
    theme(
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 6.5, color = "black", family = 'Arial',
                                   angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 6.5, color = "black", family = 'Arial'),
        legend.title = element_text(
            size = 7, color = "black", family = 'Arial'),
        legend.text = element_text(
            size = 6, color = "black", family = 'Arial'),
        legend.key.size = unit(0.3, 'cm')
    ) +
    scale_y_discrete(breaks = pos.cells, labels = ref.cells)
ggsave(filename = paste0('heatmap_MCA_', dataset, '_', method, '.png'),
       path = path.res, plot = plot.heatmap,
       units = 'cm', height = 10, width = 10)


# accessment
# import python package: sklearn.metrics
library(reticulate)
use_python('/home/zy/tools/anaconda3/bin/python3', required = T)
# py_config()
py_module_available('sklearn')
metrics <- import('sklearn.metrics')

mytable <- table(true.tags, pred.scMAGIC)

ref.names <- colnames(mytable)
all.cell <- names(table(true.tags))
uniform.names <-
    c('Astrocytes', 'Astrocytes', 'Astroglial cell', 'Endothelial cells',
      'Endothelial cells', 'Neurons', 'Ependymocytes', 'PVMs & Microglia',
      'Neurons', 'OPC', 'Neurons', 'Unassigned', 'Endothelial cells')
df.ref.names <- data.frame(ref.name = ref.names, name = uniform.names)
uniform.names <-
    c("Astrocytes", "Endothelial cells", "Ependymocytes", "Unassigned",
      "Neurons", "Oligodendrocytes", "OPC", "Unassigned", "PVMs & Microglia",
      "Unassigned", "Unassigned")
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


    # true.labels <- setdiff(unique(true.tag), 'Unassigned')
    true.labels <- unique(true.tag)
    our.tag <- scRef.tag
    weighted_macro_f1 <- metrics$f1_score(true.tag, our.tag, average = 'weighted', labels = true.labels)
    macro_f1 <- metrics$f1_score(true.tag, our.tag, average = 'macro', labels = true.labels)
    accuracy <- metrics$accuracy_score(true.tag, our.tag)
    balanced_acc <- metrics$balanced_accuracy_score(true.tag, our.tag)
    # rm unassigned in tags
    true.tag.rm <- true.tag[our.tag != 'Unassigned']
    our.tag.rm <- our.tag[our.tag != 'Unassigned']
    accuracy.rm.unassigned <- metrics$accuracy_score(true.tag.rm, our.tag.rm)
    macro_f1.rm.unassigned <- metrics$f1_score(true.tag.rm, our.tag.rm, average = 'macro', labels = unique(true.tag.rm))
    balanced.accuracy.rm.unassigned <-
        metrics$balanced_accuracy_score(true.tag.rm, our.tag.rm)

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
    out$conf <- table(true.tag, our.tag)

    return(out)

}

res.scMAGIC <- simple.evaluation(true.tags, pred.scMAGIC, df.ref.names, df.sc.names)

file.res.scMAGIC <- paste0(path.res, 'RES_MCA_', dataset, '_scMAGIC.txt')
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
