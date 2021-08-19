# import python package: sklearn.metrics
# library(reticulate)
# use_python('/home/zy/tools/anaconda3/bin/python3', required = T)
# # py_config()
# py_module_available('sklearn')
# metrics <- import('sklearn.metrics')


library(Seurat)
library(SeuratData)
data("ifnb")
dataset <- 'hPBMC_10Xv1'
ifnb <- subset(ifnb, stim == 'CTRL')
exp_sc_mat <- as.matrix(ifnb@assays$RNA@counts)
label_sc <- data.frame(
    annotation = as.character(ifnb$seurat_annotations),
    row.names = colnames(exp_sc_mat))


data("pbmcsca")
dataset <- 'hPBMC_10Xv2'
exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method %in%
                                                      c('10x Chromium (v2)')])
label_sc <- data.frame(
    annotation = as.character(pbmcsca$CellType)[pbmcsca$Method %in%
                                                    c('10x Chromium (v2)')],
    row.names = colnames(exp_sc_mat))

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
df.atlas <- as.data.frame(seurat.out.group@assays$RNA@counts)

source('/home/zy/my_git/scRef/main/scRef.v21.R')
setwd('~/my_git/scRef')
result.scref <- SCREF(exp_sc_mat, df.atlas,
                      type_ref = 'sum-counts', use.RUVseq = F, out.group = 'HCA',
                      # method1 = 'spearman',
                      GMM.floor_cutoff = 3, GMM.ceiling_cutoff = 20,
                      cluster.speed = T, CPU = 8,
                      cluster.cell = 4, min_cell = 4)
pred.scRef <- result.scref$final.out$scRef.tag


true.tags <- label_sc$annotation
table(true.tags, pred.scRef)
# df.tags <- result.scref$combine.out
# df.view <- merge(label_sc, df.tags, by = 'row.names')
# View(df.view)

library(ggplot2)
path.res <- '/home/zy/scRef/figure/atlas_anno'

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
method <- 'scMAGIC'
mytable <- table(true.tags, pred.scRef)
mytable.tag <- data.frame(colnames(mytable), df.simple[colnames(mytable), 'col.simple'])
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

# mydata$origin <- factor(mydata$origin, levels = c(rownames(mytable)))
# ref.cells <- c("Unassigned", 
#                "Astrocyte", "Astroglial cell", 
#                "Vascular endothelial cell", 
#                "Ovarian vascular surface endothelium cell",
#                "Hypothalamic ependymal cell", 
#                "Dopaminergic neurons", "Ganglion cell", "Granule neuron", "Hippocampus neuron",
#                "Interstitial macrophage", "Microglia",
#                "Oligodendrocyte precursor cell")
ref.cells <- setdiff(colnames(mytable.sum), 
                     c("activative T cell", "Proliferating  B cell", 
                       'T cell-GNLY high', 'CD8-T cell',
                       'Monocyte-IGHG4 high'))
show.labels <- c(sort(unique(df.simple$col.simple)), "Unassigned")
show.labels[!(show.labels %in% ref.cells)] <- ''
show.labels <- gsub('CD4-T cell', 'CD4+ T cell', show.labels)
mydata$annotation <- factor(mydata$annotation, 
                            levels = c(sort(unique(df.simple$col.simple)), "Unassigned"))
mydata$origin <- factor(mydata$origin, 
                        levels = c("B cell", "CD4+ T cell", "Dendritic cell",
                                   "CD14+ monocyte", "CD16+ monocyte", "Megakaryocyte", 
                                   "Cytotoxic T cell", "Natural killer cell",
                                   "Plasmacytoid dendritic cell"))

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
        axis.text.y = element_text(size = 12, color = "black", face = "bold"),
        legend.title = element_text(
            size = 12, color = "black", face = "bold"),
        legend.text = element_text(
            size = 10, color = "black", face = "bold")
    ) + 
    scale_y_discrete(labels = show.labels)
    # geom_text(aes(label = round(prop, 2)), family = "Arial", size = 2.5)
ggsave(filename = paste0('heatmap_allMCA_', dataset, '_', method, '.pdf'), 
       path = path.res, plot = plot.heatmap,
       units = 'cm', height = 20, width = 15)




