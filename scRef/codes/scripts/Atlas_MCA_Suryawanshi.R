# import python package: sklearn.metrics
library(reticulate)
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
    OUT$mat_exp <- data.filter
    OUT$label <- label.filter
    return(OUT)
    
}

path.input <- '/home/disk/scRef/HumanReference_v1/HumanPlacentaDecidua_SingleCell_Suryawanshi2018'
path.output <- '/home/zy/scRef/atlas_anno/'
dataset <- 'Suryawanshi_Placenta'
# df.meta <- read.delim(paste0(path.input, '/E-MTAB-6701.processed.2'), 
#                       row.names=1)
# df.mat <- read.delim(paste0(path.input, '/E-MTAB-6701.processed.1.geneTransed.dedupGene'), 
#                      row.names=1)
# meta.pla <- df.meta[df.meta$location == 'Placenta',]
# mat.pla <- df.mat[,rownames(meta.pla)]
# dataset <- 'Tormo_Decidua'
# df.meta <- read.delim(paste0(path.input, '/E-MTAB-6678.processed.2'), 
#                       row.names=1)
# df.mat <- read.delim(paste0(path.input, '/E-MTAB-6678.processed.1.geneTransed.dedupGene'), 
#                      row.names=1)
# meta.pla <- df.meta
# mat.pla <- df.mat
OUT <- list()
OUT$mat_exp <- mat.pla
OUT$label <- meta.pla
file.rds <- paste0(path.output, dataset, '.Rdata')
saveRDS(OUT, file = file.rds)
file.data.unlabeled <- paste0(path.input, '/raw_count/cell_exp.txt')
file.label.unlabeled <- paste0(path.input, '/raw_count/cell_meta.txt')
OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, 
                    del.label = c('Unclassified'))
# # saveRDS(OUT, file = paste0(path.output, dataset, '.Rdata'))
exp_sc_mat <- OUT$mat_exp
label_sc <- OUT$label
exp_pla <- exp_sc_mat[,is.na(exp_sc_mat['A1CF',])]
exp_pla <- exp_pla[!(is.na(exp_pla[,1])),]
label_pla <- label_sc[colnames(exp_pla),]

# OUT <- readRDS(file.rds)
# exp_sc_mat <- OUT$mat_exp
# label_sc <- OUT$label
# label_sc <- label_sc[label_sc$annotation %in% 
#                          c('Endo (f)', 'Endo (m)',
#                            'Endo L', 'Epi2', 'EVT', 'fFB1', 
#                            'Granulocytes', 'HB', 'SCT', 'VCT'),]
# vec.anno <- label_sc$annotation
# # vec.anno[vec.anno %in% c('dM1', 'dM3')] <- 'dM'
# vec.anno[vec.anno %in% c('Epi2')] <- 'Epi'
# vec.anno[vec.anno %in% c('fFB1')] <- 'fFB'
# label_sc$annotation2 <- vec.anno
# exp_sc_mat <- exp_sc_mat[,rownames(label_sc)]

# library(SeuratData)
# data("panc8")
# exp_sc_mat <- as.matrix(panc8@assays$RNA@counts[, panc8$dataset %in%
#                                                          c('indrop1', 'indrop2', 'indrop3', 'indrop4')])
# label_sc <- data.frame(
#     annotations = as.character(panc8$celltype)[panc8$dataset %in%
#                                                    c('indrop1', 'indrop2', 'indrop3', 'indrop4')],
#     row.names = colnames(OUT$data.filter))
# dataset <- 'panc8_indrop'


source('/home/zy/my_git/scRef/main/scRef.v20.R')
setwd('~/my_git/scRef')
df.atlas <- .imoprt_outgroup('HCA', normalization = F)
# df.atlas <- df.atlas[, colnames(df.atlas) != 'CD4']

library(Seurat)

result.scref <- SCREF(exp_pla, df.atlas,
                      type_ref = 'sum-counts', use.RUVseq = F, out.group = 'HCA',
                      # method1 = 'spearman',
                      GMM.floor_cutoff = 3, GMM.ceiling_cutoff = 20,
                      cluster.speed = T, 
                      cluster.cell = 5, min_cell = 5)
pred.scRef <- result.scref$final.out$scRef.tag

true.tags <- label_pla
table(true.tags, pred.scRef)
# df.tags <- result.scref$combine.out
# df.view <- merge(label_sc, df.tags, by = 'row.names')
# View(df.view)

library(ggplot2)
path.res <- '/home/zy/scRef/figure/atlas_anno'

# heatmap
method <- 'scMAGIC'
mytable <- table(true.tags, pred.scRef)
mydata <- data.frame(stringsAsFactors = F)
table.true <- table(true.tags)
for (label1 in rownames(mytable)) {
    row.sum <- table.true[label1]
    for (label2 in c(colnames(df.atlas), "Unassigned")) {
        if (label2 %in% colnames(mytable)) {
            mydata <- rbind(mydata, data.frame(origin = label1, annotation = label2, 
                                               count = mytable[label1, label2], 
                                               prop = mytable[label1, label2]/row.sum))
        } else {
            mydata <- rbind(mydata, data.frame(origin = label1, annotation = label2, 
                                               count = 0, prop = 0))
        }
    }
}
mydata$origin <- factor(mydata$origin, levels = c(rownames(mytable)))
# ref.cells <- c("Unassigned", 
#                "Astrocyte", "Astroglial cell", 
#                "Vascular endothelial cell", 
#                "Ovarian vascular surface endothelium cell",
#                "Hypothalamic ependymal cell", 
#                "Dopaminergic neurons", "Ganglion cell", "Granule neuron", "Hippocampus neuron",
#                "Interstitial macrophage", "Microglia",
#                "Oligodendrocyte precursor cell")
colnames(mytable) <- gsub('Vascular epithelial cell', 
                   'Vascular endothelial cell', colnames(mytable))
ref.cells <- setdiff(colnames(mytable), 
                     c("Megakaryocyte"))
show.labels <- c(colnames(df.atlas), "Unassigned")
show.labels[!(show.labels %in% ref.cells)] <- ''
mydata$annotation <- factor(mydata$annotation, levels = c(colnames(df.atlas), "Unassigned"))
mydata$origin <- factor(mydata$origin, 
                        levels = c("VCT", "Epi", "EVT",
                                   "HB", "Endo L", "SCT", "Endo (f)",
                                   "acinar", "macrophage",
                                   "endothelial", "activated_stellate",
                                   "quiescent_stellate", "schwann"))

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
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) + 
    scale_y_discrete(labels = show.labels)
    # geom_text(aes(label = round(prop, 2)), family = "Arial", size = 2.5)
ggsave(filename = paste0('heatmap_allMCA_', dataset, '_', method, '.png'), 
       path = path.res, plot = plot.heatmap,
       units = 'cm', height = 20, width = 15)




