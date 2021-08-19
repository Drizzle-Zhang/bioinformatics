# import python package: sklearn.metrics
# library(reticulate)
# use_python('/home/zy/tools/anaconda3/bin/python3', required = T)
# # py_config()
# py_module_available('sklearn')
# metrics <- import('sklearn.metrics')


library(Seurat)
library(SeuratData)
data("panc8")
exp_sc_mat <- as.matrix(panc8@assays$RNA@counts[, panc8$dataset %in%
                                                    c('indrop1', 'indrop2', 'indrop3', 'indrop4')])
label_sc <- data.frame(
    annotations = as.character(panc8$celltype)[panc8$dataset %in%
                                                   c('indrop1', 'indrop2', 'indrop3', 'indrop4')],
    row.names = colnames(exp_sc_mat))
dataset <- 'panc8_indrop'

library(Seurat)
library(scMAGIC)
data("HCL_ref")
time1 <- Sys.time()
output.scMAGIC <- scMAGIC(exp_sc_mat, HCL_ref, type_ref = 'sum-counts', atlas = 'HCL',
                          opt_speed = T, combine_num_cell = 5, min_cell = 25,
                          GMM.floor_cutoff = 3, GMM.ceiling_cutoff = 20,
                          use_RUVseq = F,  num_threads = 15)
time2 <- Sys.time()
time.diff <- difftime(time2, time1, units = 'mins')
pred.scMAGIC <- output.scMAGIC$scMAGIC.tag
true.tags <- label_sc[,1]
table(true.tags, pred.scMAGIC)


library(ggplot2)
path.res <- '/mdshare/node9/zy/scRef/figure/atlas_anno/'
method <- 'scMAGIC'
file.pred <- paste0(path.res, 'HCA_', dataset, '_scMAGIC.Rdata')
# saveRDS(pred.scMAGIC, file.pred)
pred.scRef <- readRDS(file.pred)
true.tags <- label_sc$annotations
table(true.tags, pred.scRef)

# simplify colnames of atlas
df.simple <- data.frame(check.names = F)
for (col in colnames(HCL_ref)) {
    col.split <- strsplit(col, split = '_')[[1]]
    col.simple <- paste(col.split[1:(length(col.split)-1)], collapse = '-')
    df.simple <- rbind(df.simple, data.frame(col = col,
                                             col.simple = col.simple))
}
row.names(df.simple) <- df.simple$col

# heatmap
true.tags <- true.tags[pred.scRef != 'Beta cell_AdultPancreas']
pred.scRef <- pred.scRef[pred.scRef != 'Beta cell_AdultPancreas']
mytable <- table(true.tags, pred.scRef)
mytable.tag <- data.frame(colnames(mytable), df.simple[colnames(mytable), 'col.simple'])
mytable.tag[mytable.tag[,1] == 'Unassigned', 2] <- 'Unassigned'
mytable.tag[,2] <- gsub('Acniar cell-ANXA4 high', 'Acniar cell', mytable.tag[,2])
# mytable.tag[,2] <- gsub('Basal cell', 'Other cells in HCL', mytable.tag[,2])
# mytable.tag[,2] <- gsub('D cell/ X/A cell', 'Other cells in HCL', mytable.tag[,2])
mytable.tag[,2] <- gsub('Endothelial cell in EMT', 'Endothelial cell', mytable.tag[,2])
mytable.tag[,2] <- gsub('Endothelial cell-FABP4 high', 'Endothelial cell', mytable.tag[,2])
mytable.tag[,2] <- gsub('Endothelial cell-IGFBP3 high', 'Endothelial cell', mytable.tag[,2])
mytable.tag[,2] <- gsub('Exocrine cell-SAA1 high', 'Exocrine cell', mytable.tag[,2])
mytable.tag[,2] <- gsub('Fibroblast', 'Stromal cell', mytable.tag[,2])
# mytable.tag[,2] <- gsub('Smooth muscle cell-PDK4 high', 'Smooth muscle cell', mytable.tag[,2])
mytable.tag[,2] <- gsub('Stromal cell-ERRFI1 high', 'Stromal cell', mytable.tag[,2])

mytable.sum <- generate_ref(mytable, mytable.tag)
mydata <- data.frame(stringsAsFactors = F)
table.true <- table(true.tags)
for (label1 in rownames(mytable.sum)) {
    row.sum <- table.true[label1]
    for (label2 in unique(c(colnames(mytable.sum), 'Other cells in HCL'))) {
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

tag_order <- c("Acniar cell", "Exocrine cell",
               "Alpha cell", "Endocrine cell",
               "Ductal cell", "Endothelial cell", "Stromal cell",
               "Unassigned", "Other cells in HCL")
mydata$annotation <- factor(mydata$annotation, levels = tag_order)
mydata$origin <- factor(mydata$origin,
                        levels = c("acinar", "alpha",
                                   "beta", "delta", "epsilon", "gamma",
                                   "ductal", "endothelial",
                                   "activated_stellate", "quiescent_stellate",
                                   "macrophage", "mast", "schwann"))

plot.heatmap <-
    ggplot(data = mydata, aes(x = origin, y = annotation)) +
    geom_tile(aes(fill = prop)) +
    scale_fill_gradient2(low = "#FFF5EE", mid = '#EE7700', high = "#B22222", midpoint = 0.5) +
    labs(fill = 'Proportion', title = 'HCL -> Human Pancreas') +
    scale_x_discrete(breaks = c("acinar", "alpha",
                                "beta", "delta", "epsilon", "gamma",
                                "ductal", "endothelial",
                                "activated_stellate", "quiescent_stellate",
                                "macrophage", "mast", "schwann"),
                     labels = c("Acinar cell", "Alpha cell", "Beta cell",
                                "Delta cell", "Epsilon cell", "Gamma cell",
                                "Dutal cell", "Endothelial Cell",
                                "ASC", "QSC",
                                "Macrophage", "Mast cell", "Schwann cell")) +
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
            size = 9, color = "black", family = 'Arial'),
        legend.text = element_text(
            size = 7.5, color = "black", family = 'Arial'),
        legend.key.size = unit(0.5, 'cm'),
        legend.position = 'bottom'
    )
ggsave(filename = paste0('heatmap_HCL_', dataset, '_', method, '.png'),
       path = path.res, plot = plot.heatmap,
       units = 'cm', height = 9, width = 9.5)


# accessment
# import python package: sklearn.metrics
library(reticulate)
use_python('/local/zy/tools/anaconda3/bin/python3', required = T)
# py_config()
py_module_available('sklearn')
metrics <- import('sklearn.metrics')

mytable <- table(true.tags, pred.scRef)

ref.names <- colnames(mytable)
all.cell <- names(table(true.tags))
uniform.names <-
    c("Exocrine cell", "Endocrine cell", "Endocrine cell",
      "ductal", "Endocrine cell", 'Endothelial cell',
      'Exocrine cell', 'Stromal cell', 'Stromal cell',
      'Stromal cell', 'Unassigned')
df.ref.names <- data.frame(ref.name = ref.names, name = uniform.names)
uniform.names <-
    c("Exocrine cell", "Stromal cell", "Endocrine cell", "Endocrine cell",
      "Endocrine cell", "ductal", "Endothelial cell", "Endocrine cell",
      'Endocrine cell', 'macrophage', 'Mast cell', 'Stromal cell', 'schwann')
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
  ref.labels <- intersect(unique(our.tag.rm), true.labels)
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

res.scMAGIC <- simple.evaluation(true.tags, pred.scRef, df.ref.names, df.sc.names)

file.res.scMAGIC <- paste0(path.res, 'RES_HCL_', dataset, '_scMAGIC.Rdata')
# write.table(res.scMAGIC, file.res.scMAGIC)
saveRDS(res.scMAGIC, file.res.scMAGIC)


# $accuracy
# [1] 0.9400163
#
# $accuracy.rm.unassigned
# [1] 0.9541077
