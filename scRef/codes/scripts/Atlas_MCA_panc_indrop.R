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
saveRDS(pred.scRef, file.pred)

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

all.cells <- sort(c(unique(df.simple$col.simple), "Unassigned"))
ref.cells <- setdiff(colnames(mytable.sum), 
                     c("Basal cell", "Beta cell", 
                       'D cell/ X/A cell', 'Endothelial cell in EMT',
                       'Endothelial cell-FABP4 high', 'Fibroblast',
                       'Smooth muscle cell-PDK4 high'))
pos <- c(-1, 1, -9, 0, 9, 7, 0, -2, 2, 0)
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
    scale_fill_gradient2(low = "#000000", high = "#FFFF00", mid = "#32CD32", midpoint = 0.5) + 
    labs(fill = 'Proportion') + 
    theme_bw() +
    theme(
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 12, color = "black", face = "bold", 
                                   angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 11, color = "black", face = "bold"),
        legend.title = element_text(
            size = 12, color = "black", face = "bold"),
        legend.text = element_text(
            size = 10, color = "black", face = "bold")
    ) + 
    scale_y_discrete(breaks = pos.cells, labels = ref.cells)
ggsave(filename = paste0('heatmap_HCA_', dataset, '_', method, '.png'), 
       path = path.res, plot = plot.heatmap,
       units = 'cm', height = 20, width = 18)




