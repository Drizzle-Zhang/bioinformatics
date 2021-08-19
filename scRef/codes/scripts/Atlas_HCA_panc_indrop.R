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


source('/home/zy/my_git/scRef/main/scRef.v20.R')
setwd('~/my_git/scRef')
# df.atlas <- .imoprt_outgroup('HCA', normalization = F)
df.out.group <- 
    read.table('/home/disk/scRef/HumanAtlas_SingleCell_Han2020/combinedHCA/HCA_combined.txt', 
               header = T, row.names = 1, sep = '\t', check.names = F)
df.out.group[is.na(df.out.group)] <- 0
seurat.out.group <- 
    CreateSeuratObject(counts = df.out.group, project = "out.group", 
                       min.cells = 1, min.features = 5000)
seurat.out.group <- 
    NormalizeData(seurat.out.group, normalization.method = "LogNormalize", 
                  scale.factor = 1e6, verbose = F)
df.atlas <- data.frame(seurat.out.group@assays$RNA@counts, check.names = F)
df.atlas <- df.atlas[, !(colnames(df.atlas) %in%
                             c('Smooth muscle cell_AdultPancreas'))]


result.scref <- SCREF(exp_sc_mat, df.atlas,
                      type_ref = 'sum-counts', use.RUVseq = F, out.group = 'HCA',
                      GMM.floor_cutoff = 3, GMM.ceiling_cutoff = 20,
                      cluster.speed = T, CPU = 4,
                      cluster.cell = 5, min_cell = 5)
pred.scRef <- result.scref$final.out$scRef.tag

true.tags <- label_sc$annotation
table(true.tags, pred.scRef)
# df.tags <- result.scref$combine.out
# df.view <- merge(label_sc, df.tags, by = 'row.names')
# View(df.view)

library(ggplot2)
path.res <- '/home/zy/scRef/figure/atlas_anno/'
method <- 'scMAGIC'
file.pred <- paste0(path.res, 'HCA_', dataset, '_scMAGIC.Rdata')
# saveRDS(pred.scRef, file.pred)
pred.scRef <- readRDS(file.pred)
true.tags <- label_sc$annotations
table(true.tags, pred.scRef)

# simplify colnames of atlas
df.simple <- data.frame(check.names = F)
for (col in colnames(df.atlas)) {
    col.split <- strsplit(col, split = '_')[[1]]
    col.simple <- paste(col.split[1:(length(col.split)-1)], collapse = '-')
    df.simple <- rbind(df.simple, data.frame(col = col,
                                             col.simple = col.simple))
}
row.names(df.simple) <- df.simple$col

# heatmap
mytable <- table(true.tags, pred.scRef)
mytable.tag <- data.frame(colnames(mytable), df.simple[colnames(mytable), 'col.simple'])
mytable.tag[mytable.tag[,1] == 'Unassigned', 2] <- 'Unassigned'
mytable.sum <- .generate_ref(mytable, mytable.tag)
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

ref.cells <- setdiff(colnames(mytable.sum), 
                     c("Basal cell", "Beta cell", 
                       'D cell/ X/A cell', 'Endothelial cell in EMT',
                       'Endothelial cell-FABP4 high', 'Fibroblast',
                       'Smooth muscle cell-PDK4 high'))
set.seed(1234)
all.cells <- sort(unique(c(sample(setdiff(unique(df.simple$col.simple), ref.cells), 100), 
                           ref.cells, "Unassigned")))
pos <- c(-1, 1, -5, -2, 1, 3, 0, -1, 1, 0)
pos.cells <- c()
for (idx in 1:length(ref.cells)) {
    pos.idx <- pos[idx]
    pos.cells <- c(pos.cells, all.cells[which(all.cells==ref.cells[idx])+pos.idx])
}
mydata$annotation <- factor(mydata$annotation, levels = c(all.cells))
mydata$origin <- factor(mydata$origin, 
                        levels = c("acinar", "alpha", "ductal",
                                   "beta", "delta", "epsilon", "gamma",
                                   "endothelial", 
                                   "activated_stellate", "quiescent_stellate", 
                                   "macrophage", "mast", "schwann"))

plot.heatmap <- 
    ggplot(data = mydata, aes(x = origin, y = annotation)) + 
    geom_tile(aes(fill = prop)) + 
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
ggsave(filename = paste0('heatmap_HCA_', dataset, '_', method, '.png'), 
       path = path.res, plot = plot.heatmap,
       units = 'cm', height = 10, width = 10)


# accessment
# import python package: sklearn.metrics
library(reticulate)
use_python('/home/zy/tools/anaconda3/bin/python3', required = T)
# py_config()
py_module_available('sklearn')
metrics <- import('sklearn.metrics')

mytable <- table(true.tags, pred.scRef)

ref.names <- colnames(mytable)
all.cell <- names(table(true.tags))
uniform.names <- 
    c("Exocrine cell", "Endocrine cell", "Basal cell", "Endocrine cell",
      "D cell/ X/A cell", "ductal", "Endocrine cell", "Endocrine cell", 
      'Endothelial cell', 'Endothelial cell', 'Endothelial cell', 'Endothelial cell', 
      'Exocrine cell', 'Fibroblast', 'Smooth muscle cell', 'Stromal cell',
      'Stromal cell', 'Unassigned')
df.ref.names <- data.frame(ref.name = ref.names, name = uniform.names)
uniform.names <- 
    c("Exocrine cell", "Stromal cell", "Endocrine cell", "Endocrine cell",
      "Endocrine cell", "ductal", "endothelial", "Endocrine cell", 
      'Endocrine cell', 'macrophage', 'Unassigned', 'Stromal cell', 'schwann')
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

res.scMAGIC <- simple.evaluation(true.tags, pred.scRef, df.ref.names, df.sc.names)

file.res.scMAGIC <- paste0(path.res, 'RES_MCA_', dataset, '_scMAGIC.txt')
# write.table(res.scMAGIC, file.res.scMAGIC)
saveRDS(res.scMAGIC, file.res.scMAGIC)


# $accuracy
# [1] 0.9400163
# 
# $accuracy.rm.unassigned
# [1] 0.9541077