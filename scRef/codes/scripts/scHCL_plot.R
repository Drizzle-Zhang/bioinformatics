library(ggplot2)
library(scMAGIC)
data("HCL_ref")
data("MCA_ref")

path.scHCL <- '/mdshare/node9/zy/scRef/scHCL/'
method <- 'scHCL'

###################### Neocortex
list.target <- readRDS('/local/zy/my_git/scMAGIC_scripts/data/MouseNeocortex.Rdata')
dataset <- 'Tasic2018'
exp_sc_mat <- list.target$mat_exp
label_sc <-list.target$label

file.brain <- paste0(path.scHCL, 'scHCL_', dataset, '.Rdata')
hcl_result_brain <- readRDS(file.brain)
true.tags <- label_sc

df.simple <- data.frame(check.names = F)
for (col in colnames(MCA_ref)) {
    col.split <- strsplit(col, split = "(", fixed = T)[[1]]
    col.simple <- col.split[1]
    # col.simple <- paste(col.split[1:(length(col.split)-1)], collapse = '-')
    df.simple <- rbind(df.simple, data.frame(col = col,
                                             col.simple = col.simple))
}
row.names(df.simple) <- df.simple$col
pred.scHCL <- hcl_result_brain$scHCL
for (cell in unique(pred.scHCL)) {
    pred.scHCL[pred.scHCL == cell] <- df.simple[cell, 'col.simple']
}
pred.scHCL[pred.scHCL %in% c('Astrocyte_Atp1b2 high')] <- 'Astrocyte'
pred.scHCL[pred.scHCL %in% c('Atrial cardiomyocyte_Acta2 high', 
                             'Ventricle cardiomyocyte_Kcnj8 high')] <- 'Cardiomyocyte'
# pred.scHCL[pred.scHCL %in% c('Dopaminergic neurons', 'Granule neurons',
#                              'Hippocampus neurons_Asic4 high', 'Neuron', 'Purkinje cell')] <- 'Neuron'
# pred.scHCL[pred.scHCL %in% c('Pyramidal neuron cell', 'Schwann cell')] <- 'Other Neuron'
pred.scHCL[pred.scHCL %in% c('Endothelial cell_Fabp4 high', 'Endothelial cell_Ly6c1 high',
                             'Endothelial cell_Tm4sf1 high', 'Endothelial cells_Vwf high')] <- 'Endothelial cell'
pred.scHCL[pred.scHCL %in% c('Muscle cell_Lrrc15 high', 'Muscle cell_Mgp high',
                             'Muscle cell_Myl9 high')] <- 'Muscle cell'
pred.scHCL[pred.scHCL %in% c('Mesenchymal stromal cell', 'Stromal cell_Car3 high', 
                             'Stromal cell_Col3a1 high', 'Stromal cell_Cxcl14 high',
                             'Stromal cell_Dcn high', 'Stromal cell_Fmod high', 'Stromal cell_Gas6 high',
                             'Stromal cell_Mfap4 high', 'Stromal cell_Ptgds high',
                             'Stromal cell_Smoc2 high')] <- 'Stromal cell'
pred.scHCL[pred.scHCL %in% c('Osteoblast_Dlk1 high')] <- 'Osteoblast'
pred.scHCL[pred.scHCL %in% c('Macrophage', 'Macrophage_C1qc high',
                             'Macrophage_Pf4 high')] <- 'Macrophage'
pred.scHCL[pred.scHCL %in% c('Smooth muscle cell_Mgp high', 
                             'Smooth muscle cell_Rgs5 high')] <- 'Smooth muscle cell'
pred.scHCL[pred.scHCL %in% c('Erythroblast_Hbb-y high', 'MSC',
                             'T cell_Pclaf high', 'Purkinje cell',
                             'Pyramidal neuron cell', 'Schwann cell',
                             'Hippocampus neurons_Asic4 high')] <- 'Other cells in MCA'

mytable.sum <- table(true.tags, pred.scHCL)
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
               "Dopaminergic neurons ", 'Granule neurons', 'Neuron',
               "Macrophage", "Microglia", "Smooth muscle cell", "Oligodendrocyte precursor cell",
               "Unassigned", 'Cardiomyocyte', "Muscle cell", 'Osteoblast', 'Stromal cell',
               "Other cells in MCA")
mydata$annotation <- factor(mydata$annotation, levels = tag_order)
mydata$origin <- factor(mydata$origin,
                        levels = c("Astrocyte", "Endothelial cell",
                                   "GABAergic Neuron", "Glutamatergic Neuron",
                                   "PVM & Microglia", "Smooth muscle cell", "Oligodendrocyte & OPC",
                                   "Pericyte", "VLMC"))

plot.heatmap_Tasic2018 <-
    ggplot(data = mydata, aes(x = origin, y = annotation)) +
    geom_tile(aes(fill = prop)) +
    # scale_fill_gradient2(low = "#C0C0C0", high = "#FFFF00", mid = "#32CD32", midpoint = 0.5) +
    scale_fill_gradient2(low = "#FFF5EE", mid = '#EE7700', high = "#B22222", midpoint = 0.5, limits = c(0,1)) +
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
ggsave(filename = paste0('heatmap_MCA_', dataset, '_', method, '.png'),
       path = path.scHCL, plot = plot.heatmap_Tasic2018,
       units = 'cm', height = 10.5, width = 9.5)


###################### Duodenum
path.input <- '/mdshare/node9/zy/scRef/sc_data/MouseSmallIntestinalEpithelium_SingleCell_Haber2018'
dataset <- 'Haber_Duodenum'
file.data.unlabeled <- paste0(path.input, '/raw_count/cell_exp.txt')
file.label.unlabeled <- paste0(path.input, '/raw_count/cell_meta.txt')
data.unlabeled <- read.delim(file.data.unlabeled, row.names=1)
label.unlabeled <- read.delim(file.label.unlabeled, row.names=1)
label.unlabeled$location <- unlist(lapply(strsplit(row.names(label.unlabeled), '_'),
                                          function(x) {x[2]}))
label.unlabeled$location[!(label.unlabeled$location %in%
                               c('Ileum', 'Duodenum', 'Jejunum'))] <- 'Other'
label.Duodenum <- label.unlabeled[label.unlabeled$location == 'Duodenum',]
mat.Duodenum <- data.unlabeled[, rownames(label.Duodenum)]
mat.Duodenum <- mat.Duodenum[!(is.na(mat.Duodenum[,1])),]
exp_sc_mat <- mat.Duodenum
label_sc <- label.Duodenum

file.Duodenum <- paste0(path.scHCL, 'scHCL_', dataset, '.Rdata')
hcl_result_Duodenum <- readRDS(file.Duodenum)
true.tags <- label_sc$Cluster
true.tags[true.tags %in% c('Goblet')] <- 'Goblet cell'
true.tags[true.tags %in% c('Enteroendocrine')] <- 'Enteroendocrine cell'
true.tags[true.tags %in% c('Paneth')] <- 'Paneth cell'
true.tags[true.tags %in% c('Tuft')] <- 'Tuft cell'
true.tags[true.tags %in% c('EP')] <- 'Enterocyte progenitor'
true.tags[true.tags %in% c('TA')] <- 'Transit-amplifying cell'
true.tags[true.tags %in% c('Stem')] <- 'Intestinal stem cell'

df.simple <- data.frame(check.names = F)
for (col in colnames(MCA_ref)) {
    col.split <- strsplit(col, split = "(", fixed = T)[[1]]
    col.simple <- col.split[1]
    # col.simple <- paste(col.split[1:(length(col.split)-1)], collapse = '-')
    df.simple <- rbind(df.simple, data.frame(col = col,
                                             col.simple = col.simple))
}
row.names(df.simple) <- df.simple$col
pred.scHCL <- hcl_result_Duodenum$scHCL
for (cell in unique(pred.scHCL)) {
    pred.scHCL[pred.scHCL == cell] <- df.simple[cell, 'col.simple']
}
pred.scHCL[pred.scHCL %in% c('Enteroendocrine')] <- 'Enteroendocrine cell'
pred.scHCL[pred.scHCL %in% c('S cell_Chgb high', 'S cell_Gip high')] <- 'S cell'
pred.scHCL[pred.scHCL %in% c('Epithelial cell_Bex1 high', 'Epithelial cell_Lgals2 high',
                             'Epithelial cell_Gkn3 high')] <- 'Epithelial cell'
pred.scHCL[pred.scHCL %in% c('Epithelium of small intestinal villi_Fabp1 high', 
                             'Epithelium of small intestinal villi_mt-Nd1 high',
                             'Epithelium of small intestinal villi_S100g high')] <- 'Epithelium of small intestinal villi'
pred.scHCL[pred.scHCL %in% c('Erythroblast_Mt2 high')] <- 'Erythroblast'
pred.scHCL[pred.scHCL %in% c('Glandular epithelium_Ltf high')] <- 'Glandular epithelium'
pred.scHCL[pred.scHCL %in% c('G cell', 'Dendritic cell', 'B cell')] <- 'Other cells in MCA'

mytable.sum <- table(true.tags, pred.scHCL)

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

tag_order <- c("Paneth cell", "Tuft cell",
               "Enteroendocrine cell", 'S cell', 'Columnar epithelium',
               "Epithelium of small intestinal villi", 
               "Acinar cell", "Endocrine progenitor cell", "Epithelial cell",
               "Erythroblast", 'Gastric mucosal cell', "Glandular epithelium", 
               "Other cells in MCA")
mydata$annotation <- factor(mydata$annotation, levels = tag_order)
mydata$origin <- factor(mydata$origin,
                        levels = c("Paneth cell", "Tuft cell", "Enteroendocrine cell",
                                   "Goblet cell", "Enterocyte",
                                   "Enterocyte progenitor", 
                                   "Intestinal stem cell", "Transit-amplifying cell"))

plot.heatmap_Haber <-
    ggplot(data = mydata, aes(x = origin, y = annotation)) +
    geom_tile(aes(fill = prop)) +
    scale_fill_gradient2(low = "#FFF5EE", mid = '#EE7700', high = "#B22222", midpoint = 0.5, limits = c(0,1)) +
    labs(fill = 'Proportion', title = 'MCA -> Mouse Duodenum') +
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
        legend.key.size = unit(0.3, 'cm')
    )
ggsave(filename = paste0('heatmap_MCA_', dataset, '_', method, '.png'),
       path = path.scHCL, plot = plot.heatmap_Haber,
       units = 'cm', height = 10, width = 11.5)



######################### panc8
library(Seurat)
library(SeuratData)
data("panc8")
dataset <- 'panc8_indrop'
exp_sc_mat <- as.matrix(panc8@assays$RNA@counts[, panc8$dataset %in%
                                                    c('indrop1', 'indrop2', 'indrop3', 'indrop4')])
label_sc <- data.frame(
    annotations = as.character(panc8$celltype)[panc8$dataset %in%
                                                   c('indrop1', 'indrop2', 'indrop3', 'indrop4')],
    row.names = colnames(exp_sc_mat))

file.panc <- paste0(path.scHCL, 'scHCL_', dataset, '.Rdata')
hcl_result_panc <- readRDS(file.panc)
true.tags <- label_sc$annotation
true.tags[true.tags %in% c('acinar')] <- 'Acinar cell'
true.tags[true.tags %in% c('activated_stellate')] <- 'ASC'
true.tags[true.tags %in% c('alpha')] <- 'Alpha cell'
true.tags[true.tags %in% c('beta')] <- 'Beta cell'
true.tags[true.tags %in% c('delta')] <- 'Delta cell'
true.tags[true.tags %in% c('ductal')] <- 'Dutal cell'
true.tags[true.tags %in% c('endothelial')] <- 'Endothelial Cell'
true.tags[true.tags %in% c('epsilon')] <- 'Epsilon cell'
true.tags[true.tags %in% c('gamma')] <- 'Gamma cell'
true.tags[true.tags %in% c('macrophage')] <- 'Macrophage'
true.tags[true.tags %in% c('mast')] <- 'Mast cell'
true.tags[true.tags %in% c('quiescent_stellate')] <- 'QSC'
true.tags[true.tags %in% c('schwann')] <- 'Schwann cell'

df.simple <- data.frame(check.names = F)
for (col in colnames(HCL_ref)) {
    col.split <- strsplit(col, split = '_')[[1]]
    col.simple <- paste(col.split[1:(length(col.split)-1)], collapse = '-')
    df.simple <- rbind(df.simple, data.frame(col = col,
                                             col.simple = col.simple))
}
row.names(df.simple) <- df.simple$col
pred.scHCL <- hcl_result_pbmc$scHCL
for (cell in unique(pred.scHCL)) {
    pred.scHCL[pred.scHCL == cell] <- df.simple[cell, 'col.simple']
}
pred.scHCL[pred.scHCL %in% c('Acinar cell-CPA1 high', 'Acinar cell-REG1B high',
                             'Acniar cell-ANXA4 high')] <- 'Acinar cell'
pred.scHCL[pred.scHCL %in% c('Adipocyte-SPP1 high')] <- 'Adipocyte'
pred.scHCL[pred.scHCL %in% c('Arterial endothelial cell', 'Endothelial cell in EMT',
                             'Endothelial cell-ACKR1 high', 'Endothelial cell-APC',
                             'Endothelial cell-CCL21 high', 'Endothelial cell-ESM1 high',
                             'Endothelial cell-FABP4 high', 'Endothelial cell-IGFBP3 high',
                             'Endothelial cell-IGFBP5 high', 'Endothelial cell-SELE high',
                             'Endothelial cell-SOCS3 high', 'Endothelial cell-SPARCL1 high',
                             'Endothelial cell-TMEM100 high', 'Endothelial cell-VWF high',
                             'Vascular endothelial cell-IGFBP3 high')] <- 'Endothelial cell'
pred.scHCL[pred.scHCL %in% c('Basal cell-KRT6A high')] <- 'Basal  cell'
pred.scHCL[pred.scHCL %in% c('Chromaffin cell-VIP high')] <- 'Chromaffin cell'
pred.scHCL[pred.scHCL %in% c('Actived T cell', 'CD8 T cell',
                             'Proliferating T cell', 'T cell-GNLY high', 'T cell-IL7R high',
                             'T cell-TRAC high')] <- 'T cell'
pred.scHCL[pred.scHCL %in% c('Epithelial cell-KRT17 high', 'Luminal epithelium ', 
                             'Mucous Epithelial cell-REG1A high', 
                             'Mucous Epithelial cell-TFF1 high', 
                             'Secretory epithelial cell',
                             'Unknown Epithelial cell-FOS high',
                             'Enterocyte-SLC26A3 high')] <- 'Epithelial cell'
pred.scHCL[pred.scHCL %in% c('Exocrine cell-SAA1 high')] <- 'Exocrine cell'
pred.scHCL[pred.scHCL %in% c('Fibroblast-APOD high', 'Fibroblast-PTX3 high', 
                             'Myofibroblast-POSTN high')] <- 'Fibroblast'
pred.scHCL[pred.scHCL %in% c('M1 Macrophage', 'Macrophage-C1QB high',
                             'Macrophage-RNASE1 high')] <- 'Macrophage'
pred.scHCL[pred.scHCL %in% c('Neutrophil-IL1B high')] <- 'Neutrophil'
pred.scHCL[pred.scHCL %in% c('Neuroendocrine cell-SST high')] <- 'Neuroendocrine cell'
pred.scHCL[pred.scHCL %in% c('Pit cell-FOXQ1 high')] <- 'Pit cell'
pred.scHCL[pred.scHCL %in% c('Smooth muscel cell', 'Smooth muscle cell-CCL19 high',
                             'Smooth muscle cell-MYL9 high', 
                             'Smooth muscle cell-PDK4 high')] <- 'Smooth muscle cell'
pred.scHCL[pred.scHCL %in% c('Stromal cell-ERRFI1 high', 'Stromal cell-LUM high',
                             'Stromal cell-PLA2G2A high', 
                             'Smooth muscle cell-PDK4 high')] <- 'Stromal cell'
pred.scHCL[pred.scHCL %in% c('Alveolar bipotent/intermediate cell', 'Cervical Mesothelial cell', 
                             'Chondrocyte', 'Club cell', 'Dendritic cell', 'Pit cell',
                             'Gastric mucosa cell', 'Inflammatory cell')] <- 'Other cells in HCL'
mytable.sum <- table(true.tags, pred.scHCL)

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

tag_order <- c("Acinar cell", "Exocrine cell",
               "Alpha cell", "Beta cell", "Endocrine cell",
               "Ductal cell", "Endothelial cell", 
               "Stromal cell", "Macrophage", "Mast cell",
               "Chromaffin cell", "Neuroendocrine cell", "D cell/ X/A cell", 
               "Endocervix Mesothelial cell", "Endometrial cell", "Epithelial cell", 
               "Goblet cell", "Basal  cell",
               "Fibroblast", "Smooth muscle cell", "Adipocyte", "Monocyte", "Neutrophil",
               "Antigen-presenting cell", "T cell", "Unknown",
               "Other cells in HCL")
mydata$annotation <- factor(mydata$annotation, levels = tag_order)
mydata$origin <- factor(mydata$origin,
                        levels = c("Acinar cell", "Alpha cell", "Beta cell",
                                   "Delta cell", "Epsilon cell", "Gamma cell",
                                   "Dutal cell", "Endothelial Cell", "Macrophage", "Mast cell",
                                   "ASC", "QSC", "Schwann cell"))

plot.heatmap_panc <-
    ggplot(data = mydata, aes(x = origin, y = annotation)) +
    geom_tile(aes(fill = prop)) +
    scale_fill_gradient2(low = "#FFF5EE", mid = '#EE7700', high = "#B22222", 
                         midpoint = 0.5, limits = c(0,1)) +
    labs(fill = 'Proportion', title = 'HCL -> Human Pancreas') +
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
        legend.position = 'none'
    )
ggsave(filename = paste0('heatmap_HCL_', dataset, '_', method, '.png'),
       path = path.scHCL, plot = plot.heatmap_panc,
       units = 'cm', height = 15, width = 11)


#################### pbmcsca
library(Seurat)
library(SeuratData)
data("pbmcsca")
dataset <- 'hPBMC_10Xv2'
exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method %in%
                                                      c('10x Chromium (v2)')])
label_sc <- data.frame(
    annotation = as.character(pbmcsca$CellType)[pbmcsca$Method %in%
                                                    c('10x Chromium (v2)')],
    row.names = colnames(exp_sc_mat))

file.pbmc <- paste0(path.scHCL, 'scHCL_', dataset, '.Rdata')
hcl_result_pbmc <- readRDS(file.pbmc)
true.tags <- label_sc$annotation
true.tags[true.tags %in% c('CD14+ monocyte', 'CD16+ monocyte')] <- 'Monocyte'
true.tags[true.tags %in% c('CD4+ T cell', 'Cytotoxic T cell')] <- 'T cell'
true.tags[true.tags %in% c('Plasmacytoid dendritic cell')] <- 'Dendritic cell'

df.simple <- data.frame(check.names = F)
for (col in colnames(HCL_ref)) {
    col.split <- strsplit(col, split = '_')[[1]]
    col.simple <- paste(col.split[1:(length(col.split)-1)], collapse = '-')
    df.simple <- rbind(df.simple, data.frame(col = col,
                                             col.simple = col.simple))
}
row.names(df.simple) <- df.simple$col
pred.scHCL <- hcl_result_pbmc$scHCL
for (cell in unique(pred.scHCL)) {
    pred.scHCL[pred.scHCL == cell] <- df.simple[cell, 'col.simple']
}
pred.scHCL[pred.scHCL %in% c('B cell (Centrocyte)', 'B cell(Centrocyte)',
                             'B cell(Plasmocyte)', 'Proliferating  B cell')] <- 'B cell'
pred.scHCL[pred.scHCL %in% c('Dendritic cell-FCER1A high')] <- 'Dendritic cell'
pred.scHCL[pred.scHCL %in% c("Megakaryocyte/Erythtoid progenitor cell")] <- 'MEPC'
pred.scHCL[pred.scHCL %in% c('Monocyte-FCGR3A high', 'Monocyte-IGHG4 high',
                             'Monocyte-TPPP3 high')] <- 'Monocyte'
# pred.scHCL[pred.scHCL %in% c('NK cell', 'CD16+ monocyte')] <- 'Natural killer cell'
pred.scHCL[pred.scHCL %in% c('activative T cell', 'Actived T cell', 'Proliferating T cell',
                             'T cell-CCL5 high', 'T cell-GNLY high', 'T cell-IL7R high',
                             'T cell-TRAC high')] <- 'T cell'
pred.scHCL[pred.scHCL %in% c('Arterial endothelial cell-GJA5 high', 'Multipotential  progenitor cell',
                             'Proliferating cell')] <- 'Other cells in HCL'
pred.scHCL[pred.scHCL %in% c('Macrophage-CXCL2 high', 'Macrophage-FCGR3A high',
                             'Macrophage-HLA-DRA high', 'Motile liver macrophage')] <- 'Macrophage'
pred.scHCL[pred.scHCL %in% c('Neutrophil-ELANE high', 'Neutrophil-MMP high',
                             'Neutrophil-MMP9 high')] <- 'Neutrophil'
mytable.sum <- table(true.tags, pred.scHCL)
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

tag_order <- c("B cell", "T cell", "Dendritic cell", 
               "Conventional dendritic cell", "Plasmacytoid dendritic cell", 
               "Megakaryocyte", "Monocyte", "NK cell",
               "Macrophage", "MEPC", "Neutrophil",
               "Other cells in HCL")
mydata$annotation <- factor(mydata$annotation, levels = tag_order)
mydata$origin <- factor(mydata$origin,
                        levels = c("B cell", "T cell",
                                   "Dendritic cell", "Megakaryocyte",
                                   "Monocyte",
                                   "Natural killer cell"))

plot.heatmap_pbmc <-
    ggplot(data = mydata, aes(x = origin, y = annotation)) +
    geom_tile(aes(fill = prop)) +
    scale_fill_gradient2(low = "#FFF5EE", mid = '#EE7700', high = "#B22222", midpoint = 0.5) +
    labs(fill = 'Proportion', title = 'HCL -> Human PBMC') +
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
        legend.position = 'right'
    )
ggsave(filename = paste0('heatmap_HCL_', dataset, '_', method, '.png'),
       path = path.scHCL, plot = plot.heatmap_pbmc,
       units = 'cm', height = 8.5, width = 9.2)


