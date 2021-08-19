path.noatlas <- '/mdshare/node9/zy/scRef/no_atlas/'
if (!file.exists(path.noatlas)) {
    dir.create(path.noatlas)
}

source('/local/zy/my_git/scRef/main/scRef.v22.R')
setwd('~/my_git/scRef')

### cross platform
library(Seurat)
library(SeuratData)
data("panc8")

ref.dataset <- 'celseq2'
ref.mtx <- as.matrix(panc8@assays$RNA@counts[, panc8$dataset %in% c('celseq2')])
ref.labels <- as.character(panc8$celltype)[panc8$dataset %in% c('celseq2')]
dataset <- 'indrop'
exp_sc_mat <- as.matrix(panc8@assays$RNA@counts[, panc8$dataset %in%
                                                    c('indrop1', 'indrop2', 'indrop3', 'indrop4')])
label_sc <- as.character(panc8$celltype)[panc8$dataset %in%
                                             c('indrop1', 'indrop2', 'indrop3', 'indrop4')]

saveRDS(ref.labels, paste0(path.noatlas, 'ref_labels_', ref.dataset, '.Rdata'))
saveRDS(label_sc, paste0(path.noatlas, 'true_labels_', dataset, '.Rdata'))

result.scref <- SCREF(exp_sc_mat, ref.mtx, ref.labels,
                      type_ref = 'sc-counts', use.RUVseq = T, 
                      out.group = 'HCA',
                      cluster.speed = T, cluster.cell = 3,
                      GMM.floor_cutoff = 3, GMM.ceiling_cutoff = 20,
                      threshold.recall = 0.5,
                      min_cell = 1, CPU = 10)
pred.scMAGIC.noatlas <- result.scref$final.out$scRef.tag
table(label_sc, pred.scMAGIC.noatlas)
rda.scMAGIC.noatlas <- paste0(path.noatlas, ref.dataset, '_', dataset, '_scMAGIC_noatlas.Rdata')
saveRDS(pred.scMAGIC.noatlas, file = rda.scMAGIC.noatlas)

### 2nd situation
path.output <- '/mdshare/node9/zy/scRef/Benchmark/mouse_brain/'
dataset <- 'Tasic'
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_Tasic <- OUT$mat_exp
label_Tasic <- OUT$label
ref.labels <-label_Tasic[,1]
ref.mtx <- exp_Tasic
ref.dataset <- 'Tasic'

file.in <- '/mdshare/node9/zy/scRef/sc_data/Hochgerner/GSE95315_10X_expression_data_v2.tab'
path.output <- '/mdshare/node9/zy/scRef/Benchmark/mouse_brain/'
dataset <- 'HochgernerA'
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_Hochgerner <- OUT$mat_exp
label_Hochgerner <- OUT$label
exp_sc_mat <- exp_Hochgerner
label_sc <- label_Hochgerner

saveRDS(ref.labels, paste0(path.noatlas, 'ref_labels_', ref.dataset, '.Rdata'))
saveRDS(label_sc[,2], paste0(path.noatlas, 'true_labels_', dataset, '.Rdata'))


result.scref <- SCREF(exp_sc_mat, ref.mtx, ref.labels,
                      type_ref = 'sc-counts', CPU = 8)
pred.scMAGIC.noatlas <- result.scref$final.out$scRef.tag
rda.scMAGIC.noatlas <- paste0(path.noatlas, ref.dataset, '_', dataset, '_scMAGIC_noatlas.Rdata')
saveRDS(pred.scMAGIC.noatlas, file = rda.scMAGIC.noatlas)
table(label_sc[,2], pred.scMAGIC.noatlas)


### 3rd situation
library(stringr)
file.mtx <- '/mdshare/node9/zy/scRef/sc_data/MouseAtlas_SingleCell_Han2018/MCA_by_tissue/Brain/Count_all_batch.txt'
df.mtx <- read.delim(file.mtx, stringsAsFactors = F, row.names = 1)
file.cellid <- '/mdshare/node9/zy/scRef/sc_data/MouseAtlas_SingleCell_Han2018/MCA_by_tissue/Brain/CellAssignments_all_batch.txt'
df.cellid <- read.delim(file.cellid, stringsAsFactors = F, row.names = 1)
df.labels <- df.cellid[colnames(df.mtx),]
ref.labels <- df.labels$CellType
ref.mtx <- df.mtx
ref.labels <- df.labels$CellType[df.labels$CellType != 'Pan-GABAergic']
ref.mtx <- df.mtx[, df.labels$CellType != 'Pan-GABAergic']
ref.dataset <- 'MCA_brain'

path.output <- '/mdshare/node9/zy/scRef/Benchmark/mouse_brain/'
dataset <- 'Campbell'
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_sc_mat <- OUT$mat_exp
label_sc <- OUT$label

saveRDS(ref.labels, paste0(path.noatlas, 'ref_labels_', ref.dataset, '.Rdata'))
saveRDS(label_sc$CellType, paste0(path.noatlas, 'true_labels_', dataset, '.Rdata'))

result.scref <- SCREF(exp_sc_mat, ref.mtx, ref.labels,
                      type_ref = 'sc-counts', use.RUVseq = T, 
                      cluster.speed = T, cluster.cell = 5, 
                      GMM.ceiling_cutoff = 10,
                      min_cell = 10, CPU = 8)
pred.scMAGIC.noatlas <- result.scref$final.out$scRef.tag
rda.scMAGIC.noatlas <- paste0(path.noatlas, ref.dataset, '_', dataset, '_scMAGIC_noatlas.Rdata')
saveRDS(pred.scMAGIC.noatlas, file = rda.scMAGIC.noatlas)
table(label_sc$CellType, pred.scMAGIC.noatlas)

### cross-species
path.output <- '/mdshare/node9/zy/scRef/Benchmark/cross_species/'
ref.dataset <- 'BaronM'
file.save <- paste0(path.output, ref.dataset, '.Rdata')
OUT <- readRDS(file.save)
ref.labels <- OUT$label[,1]
ref.mtx <- OUT$mat_exp

dataset <- 'panc8_celseq2'
library(Seurat)
file.save <- paste0(path.output, dataset, '.Rdata')
OUT <- readRDS(file.save)
exp_sc_mat <- OUT$mat_exp
label_sc <- OUT$label$annotations
exp_sc_mat <- transform.HomoloGene(exp_sc_mat)

saveRDS(ref.labels, paste0(path.noatlas, 'ref_labels_', ref.dataset, '.Rdata'))
saveRDS(label_sc, paste0(path.noatlas, 'true_labels_', dataset, '.Rdata'))

result.scref <- SCREF(exp_sc_mat, ref.mtx, ref.labels, use.RUVseq = T,
                      cluster.speed = F, cluster.resolution = 1,
                      GMM.floor_cutoff = 2, GMM.ceiling_cutoff = 20, CPU = 8)
pred.scMAGIC.noatlas <- result.scref$final.out$scRef.tag
rda.scMAGIC.noatlas <- paste0(path.noatlas, ref.dataset, '_', dataset, '_scMAGIC_noatlas.Rdata')
saveRDS(pred.scMAGIC.noatlas, file = rda.scMAGIC.noatlas)
table(label_sc, pred.scMAGIC.noatlas)

path.output <- '/mdshare/node9/zy/scRef/Benchmark/cross_species/'
ref.dataset <- 'MCA'
file.save <- paste0(path.output, ref.dataset, '1.Rdata')
OUT <- readRDS(file.save)
ref.labels <- OUT$label[,1]
ref.mtx <- OUT$mat_exp
dataset <- 'pbmcsca_10Xv2'
library(Seurat)
library(SeuratData)
data("pbmcsca")
exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
exp_sc_mat <- transform.HomoloGene(exp_sc_mat)
result.scref <- SCREF(exp_sc_mat, ref.mtx, ref.labels, use.RUVseq = T,
                      GMM.floor_cutoff = 3, GMM.ceiling_cutoff = 20, CPU = 10)
pred.scMAGIC.noatlas <- result.scref$final.out$scRef.tag
rda.scMAGIC.noatlas <- paste0(path.noatlas, ref.dataset, '_', dataset, '_scMAGIC_noatlas.Rdata')
saveRDS(pred.scMAGIC.noatlas, file = rda.scMAGIC.noatlas)
# table(label_sc, pred.scMAGIC.noatlas)



#### assessment
# functions of evaluation
simple.evaluation <- function(true.tag, scRef.tag, df.ref.names, df.sc.names) {
    # import python package
    library(reticulate)
    use_python('/local/zy/tools/anaconda3/bin/python3', required = T)
    # py_config()
    py_module_available('sklearn')
    metrics <- import('sklearn.metrics')
    # uniform tags
    for (j in 1:dim(df.ref.names)[1]) {
        scRef.tag[scRef.tag == df.ref.names[j, 'ref.name']] <- df.ref.names[j, 'name']
    }
    scRef.tag[!(scRef.tag %in% df.ref.names$name)] <- 'Unassigned'
    for (j in 1:dim(df.sc.names)[1]) {
        true.tag[true.tag == df.sc.names[j, 'sc.name']] <- df.sc.names[j, 'name']
    }

    percent.unassigned <- sum(scRef.tag == 'Unassigned')/sum(true.tag == 'Unassigned')
    tmp.true.tag <- true.tag
    tmp.our.tag <- scRef.tag
    tmp.true.tag[tmp.true.tag != 'Unassigned'] <- '0'
    tmp.our.tag[tmp.our.tag != 'Unassigned'] <- '0'
    f1.unassigned <- metrics$f1_score(tmp.true.tag, tmp.our.tag, average = 'binary', pos_label = 'Unassigned')
    

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
    out$f1.unassigned <- f1.unassigned
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

ref.datasets <- c('celseq2', 'Tasic', 'MCA_brain', 'MCA')
target.datasets <- c('indrop', 'HochgernerA', 'Campbell', 'pbmcsca_10Xv2')
folder.two <- c('human_panc/', 'mouse_brain/', 'mouse_brain/', 'cross_species/')
path.ben <- '/mdshare/node9/zy/scRef/Benchmark/'
path.noatlas <- '/mdshare/node9/zy/scRef/no_atlas/'
df.plot <- data.frame()

# 1st situation
i <- 1
ref.dataset <- ref.datasets[i]
dataset <- target.datasets[i]

ref.labels <- readRDS(paste0(path.noatlas, 'ref_labels_', ref.dataset, '.Rdata'))
true.labels <- readRDS(paste0(path.noatlas, 'true_labels_', dataset, '.Rdata'))

ref.names <- names(table(ref.labels))
all.cell <- names(table(true.labels))
df.ref.names <- data.frame(ref.name = ref.names, name = ref.names)
df.sc.names <- data.frame(sc.name = all.cell, name = all.cell)

# no atlas
rda.scMAGIC.noatlas <- paste0(path.noatlas, ref.dataset, '_', dataset, '_scMAGIC_noatlas.Rdata')
pred.scMAGIC.noatlas <- readRDS(rda.scMAGIC.noatlas)
res.noatlas.1 <- simple.evaluation(true.labels, pred.scMAGIC.noatlas, df.ref.names, df.sc.names)
df.plot <- rbind(df.plot, data.frame(strategy = 'no atlas',
                                     data = '1st situation',
                                     Dataset = 'Muraro -> Baron, Human pancreas',
                                     accuracy = res.noatlas.1$accuracy,
                                     f1.unassigned = res.noatlas.1$f1.unassigned))

# use atlas
path.output <- paste0(path.ben, folder.two[i])
rda.scMAGIC <- paste0(path.output, ref.dataset, '_', dataset, '_scRef.Rdata')
pred.scMAGIC <- readRDS(rda.scMAGIC)
res.two.1 <- simple.evaluation(true.labels, pred.scMAGIC, df.ref.names, df.sc.names)
df.plot <- rbind(df.plot, data.frame(strategy = 'use atlas', data = '1st situation',
                                     Dataset = 'Muraro -> Baron, Human pancreas',
                                     accuracy = res.two.1$accuracy,
                                     f1.unassigned = res.two.1$f1.unassigned))


# 2nd situation
i <- 2
ref.dataset <- ref.datasets[i]
dataset <- target.datasets[i]

ref.labels <- readRDS(paste0(path.noatlas, 'ref_labels_', ref.dataset, '.Rdata'))
true.labels <- readRDS(paste0(path.noatlas, 'true_labels_', dataset, '.Rdata'))

ref.names <- unique(ref.labels)
all.cell <- unique(true.labels)
uniform.names <- c("Neuron", "Endothelial Cell", "Astrocyte", "Microglia",
                   "Oligodendrocyte", "OPC")
df.ref.names <- data.frame(ref.name = ref.names, name = uniform.names)
uniform.names <- c("Endothelial Cell", "Unassigned", "Unassigned", "Microglia", "Unassigned",
                   "Oligodendrocyte", "Unassigned", "OPC", "Astrocyte", "Unassigned",
                   "Unassigned", "Unassigned", "Neuron", "Neuron",
                   "Neuron", "Neuron", "Neuron", "Neuron", "Neuron", "Unassigned")
df.sc.names <- data.frame(sc.name = all.cell, name = uniform.names)

# noatlas
rda.scMAGIC.noatlas <- paste0(path.noatlas, ref.dataset, '_', dataset, '_scMAGIC_noatlas.Rdata')
pred.scMAGIC.noatlas <- readRDS(rda.scMAGIC.noatlas)
res.noatlas.2 <- simple.evaluation(true.labels, pred.scMAGIC.noatlas, df.ref.names, df.sc.names)
df.plot <- rbind(df.plot, data.frame(strategy = 'no atlas', data = '2nd situation',
                                     Dataset = 'Tasic -> Hochgerner, Mouse brain',
                                     accuracy = res.noatlas.2$accuracy,
                                     f1.unassigned = res.noatlas.2$f1.unassigned))

# two-round
path.output <- paste0(path.ben, folder.two[i])
rda.scMAGIC <- paste0(path.output, ref.dataset, '_', dataset, '_scMAGIC.Rdata')
pred.scMAGIC <- readRDS(rda.scMAGIC)
res.two.2 <- simple.evaluation(true.labels, pred.scMAGIC, df.ref.names, df.sc.names)
df.plot <- rbind(df.plot, data.frame(strategy = 'use atlas', data = '2nd situation',
                                     Dataset = 'Tasic -> Hochgerner, Mouse brain',
                                     accuracy = res.two.2$accuracy,
                                     f1.unassigned = res.two.2$f1.unassigned))


# 3rd situation
i <- 3
ref.dataset <- ref.datasets[i]
dataset <- target.datasets[i]

ref.labels <- readRDS(paste0(path.noatlas, 'ref_labels_', ref.dataset, '.Rdata'))
true.labels <- readRDS(paste0(path.noatlas, 'true_labels_', dataset, '.Rdata'))

ref.names <- unique(ref.labels)
uniform.names <- c("Oligodend", "PVM/Microglia", "Astrocyte", "Neurons",
                   "PVM/Microglia", "Granulocyte", "OPC",
                   "Schwann cell", "Astrocyte", "Ependymocytes")
df.ref.names <- data.frame(ref.name = ref.names, name = uniform.names)
all.cell <- unique(true.labels)
uniform.names <- c("Neurons", "Unassigned", "Unassigned", "Ependymocytes", "Oligodend",
                   "Unassigned", "Unassigned", "Unassigned", "Astrocyte", "PVM/Microglia", "OPC")
df.sc.names <- data.frame(sc.name = all.cell, name = uniform.names)

# noatlas
rda.scMAGIC.noatlas <- paste0(path.noatlas, ref.dataset, '_', dataset, '_scMAGIC_noatlas.Rdata')
pred.scMAGIC.noatlas <- readRDS(rda.scMAGIC.noatlas)
res.noatlas.3 <- simple.evaluation(true.labels, pred.scMAGIC.noatlas, df.ref.names, df.sc.names)
df.plot <- rbind(df.plot, data.frame(strategy = 'no atlas', data = '2nd situation',
                                     Dataset = 'Han -> Campbell, Mouse brain',
                                     accuracy = res.noatlas.3$accuracy,
                                     f1.unassigned = res.noatlas.3$f1.unassigned))

# two-round
path.output <- paste0(path.ben, folder.two[i])
rda.scMAGIC <- paste0(path.output, 'MCA_', dataset, '_scMAGIC.Rdata')
pred.scMAGIC <- readRDS(rda.scMAGIC)
res.two.3 <- simple.evaluation(true.labels, pred.scMAGIC, df.ref.names, df.sc.names)
df.plot <- rbind(df.plot, data.frame(strategy = 'use atlas', data = '2nd situation',
                                     Dataset = 'Han -> Campbell, Mouse brain',
                                     accuracy = res.two.3$accuracy,
                                     f1.unassigned = res.two.3$f1.unassigned))


# 3rd situation
i <- 4
ref.dataset <- ref.datasets[i]
dataset <- target.datasets[i]

ref.labels <- readRDS(paste0(path.noatlas, 'ref_labels_', ref.dataset, '.Rdata'))
true.labels <- readRDS(paste0(path.noatlas, 'true_labels_', dataset, '.Rdata'))

ref.names <- names(table(ref.labels))
all.cell <- names(table(true.labels))
df.ref.names <- data.frame(ref.name = ref.names, name = ref.names)
uniform.names <- c("B cell", "Monocyte", "Monocyte", "T cell",
                   "T cell", "Dendritic cell", "Unassigned", "NK cell",
                   "Dendritic cell")
df.sc.names <- data.frame(sc.name = all.cell, name = uniform.names)

# noatlas
rda.scMAGIC.noatlas <- paste0(path.noatlas, ref.dataset, '_', dataset, '_scMAGIC_noatlas.Rdata')
pred.scMAGIC.noatlas <- readRDS(rda.scMAGIC.noatlas)
res.noatlas.4 <- simple.evaluation(true.labels, pred.scMAGIC.noatlas, df.ref.names, df.sc.names)
df.plot <- rbind(df.plot, data.frame(strategy = 'no atlas', data = '3rd situation',
                                     Dataset = 'Han -> Ding, Human PBMC',
                                     accuracy = res.noatlas.4$accuracy,
                                     f1.unassigned = res.noatlas.4$f1.unassigned))

# two-round
path.output <- paste0(path.ben, folder.two[i])
rda.scMAGIC <- paste0(path.output, ref.dataset, '_', dataset, '_scMAGIC.Rdata')
pred.scMAGIC <- readRDS(rda.scMAGIC)
res.two.4 <- simple.evaluation(true.labels, pred.scMAGIC, df.ref.names, df.sc.names)
df.plot <- rbind(df.plot, data.frame(strategy = 'use atlas', data = '3rd situation',
                                     Dataset = 'Han -> Ding, Human PBMC',
                                     accuracy = res.two.4$accuracy,
                                     f1.unassigned = res.two.4$f1.unassigned))
df.plot$label <- paste(df.plot$Dataset, df.plot$data, sep = ' (')
df.plot$label <- paste0(df.plot$label, ')')

write.table(df.plot, file = paste0(path.noatlas, 'no_atlas.txt'), sep = '\t')
path.tworound <- '/mdshare/node9/zy/scRef/no_atlas/'
df.plot <- read.table(paste0(path.noatlas, 'no_atlas.txt'), sep = '\t')
df.plot$label <- factor(df.plot$label, levels = rev(c('Muraro -> Baron, Human pancreas (1st situation)',
                                                      'Tasic -> Hochgerner, Mouse brain (2nd situation)',
                                                      'Han -> Campbell, Mouse brain (2nd situation)',
                                                      'Han -> Ding, Human PBMC (3rd situation)')))

# bar plot
library(ggplot2)
library(RColorBrewer)
# brewer.pal(11,"Set3")
plot.bar <-
    ggplot(df.plot, aes(x = label, y = accuracy, color = strategy, fill = strategy)) +
    geom_bar(position = 'dodge', stat = 'identity') +
    scale_color_manual(breaks = c('no atlas', 'use atlas'),
                       values = c("#80B1D3", "#BC80BD"),
                       labels = c(expression(paste("scMAGIC"['no-Atlas'])),
                                  "scMAGIC")) +
    scale_fill_manual(breaks = c('no atlas', 'use atlas'),
                       values = c("#80B1D3", "#BC80BD"),
                      labels = c(expression(paste("scMAGIC"['no-Atlas'])),
                                 "scMAGIC")) +
    labs(title = "", y = 'Accuracy', x = '', color = '', fill = '') +
    coord_flip() +
    theme_bw() +
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 7, color = 'black', family = 'Arial'),
          axis.text.y = element_text(size = 10, color = 'black', family = 'Arial'),
          axis.title = element_text(size = 10, family = 'Arial'),
          legend.text = element_text(size = 10, family = 'Arial'),
          legend.position = 'bottom',
          # legend.key.size = unit(1, 'cm'),
          legend.key = element_blank())
ggsave(filename = 'No_atlas_accuracy.png',
       path = path.noatlas, plot = plot.bar,
       units = 'cm', height = 8, width = 15)

# df.plot.f1 <- df.plot[df.plot$Dataset != 'Muraro -> Baron, Human pancreas',]
df.plot.f1 <- df.plot
plot.bar <-
    ggplot(df.plot.f1, aes(x = label, y = f1.unassigned, color = strategy, fill = strategy)) +
    geom_bar(position = 'dodge', stat = 'identity') +
    scale_color_manual(breaks = c('no atlas', 'use atlas'),
                       values = c("#80B1D3", "#FFA07A")) +
    scale_fill_manual(breaks = c('no atlas', 'use atlas'),
                      values = c("#80B1D3", "#FFA07A")) +
    labs(title = "", y = 'F1 of Unassigned', x = '', color = '', fill = '') +
    coord_flip() +
    theme_bw() +
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          axis.text.y = element_text(size = 9, family = 'Arial'),
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 7, color = 'black', family = 'Arial'),
          axis.title = element_text(size = 10, family = 'Arial'),
          legend.text = element_text(size = 8, family = 'Arial'),
          legend.position = 'bottom',
          # legend.key.size = unit(1, 'cm'),
          legend.key = element_blank())
ggsave(filename = 'No_atlas_F1.png',
       path = path.noatlas, plot = plot.bar,
       units = 'cm', height = 8, width = 14)

