library(Seurat)

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

file.data.unlabeled <- '/home/yzj/JingMA/Tang/data/refEXP.txt'
file.label.unlabeled <- '/home/yzj/JingMA/Tang/data/META.txt'
# OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, del.label = c('Unclassified'))
# saveRDS(OUT, file = '/home/zy/scRef/ear/ear.Rdata')
OUT <- readRDS('/home/zy/scRef/ear/ear.Rdata')
exp_Tasic <- OUT$mat_exp
label_Tasic <- OUT$label
ref.labels <- label_Tasic[, 1]
ref.mtx <- exp_Tasic

file.ear <- '/home/disk/drizzle/wgk/data/ear.filtered.Rdata'
seurat.ear.liger <- readRDS(file.ear)
cells <- colnames(seurat.ear.liger@assays$RNA@counts)

# first batch
vec.first <- c('C1', 'C2', 'C3', 'M1', 'M2')
first.cell <- cells[seurat.ear.liger$sample %in% vec.first]
exp_sc_mat <- subset(seurat.ear.liger, cells = first.cell)@assays$RNA@counts

source('/home/zy/my_git/scRef/main/scRef.v20.R')
setwd('~/my_git/scRef')
result.scref <- SCREF(exp_sc_mat, ref.mtx, ref.labels,
                      type_ref = 'sc-counts', use.RUVseq = T, out.group = 'HCA',
                      cluster.speed = T, cluster.cell = 20,
                      method1 = 'spearman', corr_use_HVGene1 = 5000, corr_use_HVGene2 = 5000,
                      GMM.ceiling_cutoff = 20,
                      min_cell = 10, CPU = 4)
pred.scRef <- result.scref$final.out
file.label.1 <- '/home/zy/scRef/ear/scMAGIC.3.Rdata'
saveRDS(pred.scRef, file = file.label.1)

result.scref.2 <- SCREF(exp_sc_mat, ref.mtx, ref.labels,
                        type_ref = 'sc-counts', use.RUVseq = T, out.group = 'HCA',
                        cluster.speed = T, cluster.cell = 10,
                        GMM.floor_cutoff = 3, GMM.ceiling_cutoff = 20,
                        min_cell = 10, CPU = 8)
pred.scRef.2 <- result.scref.2$final.out$scRef.tag

result.scref.3 <- SCREF(exp_sc_mat, ref.mtx, ref.labels,
                        type_ref = 'sc-counts', use.RUVseq = T, out.group = 'HCA',
                        cluster.speed = T, cluster.cell = 10,
                        GMM.floor_cutoff = 2, GMM.ceiling_cutoff = 15,
                        min_cell = 10, CPU = 8)
pred.scRef.3 <- result.scref.3$final.out$scRef.tag

### sciBet
suppressMessages(library(tidyverse))
suppressMessages(library(scibet))
suppressMessages(library(viridis))
suppressMessages(library(ggsci))
train_set <- as.data.frame(t(ref.mtx))
train_set$label <- ref.labels
test_set <- as.data.frame(t(exp_sc_mat))
sciBet <- SciBet(train_set, test_set)

## plot
seurat.unlabeled <- seurat.ear
seurat.unlabeled@meta.data$scRef.tag <- pred.scRef.3
seurat.unlabeled@meta.data$sciBet <- sciBet

DimPlot(seurat.unlabeled, reduction = "umap", label = T, repel = T, group.by = 'scRef.tag')

DimPlot(seurat.unlabeled, reduction = "umap", label = T, repel = T, group.by = 'sciBet')

plot.umap.scRef <- 
    DimPlot(seurat.unlabeled, reduction = "umap", label = T, repel = T, group.by = 'scRef.tag') + 
    scale_color_manual(values = c(hue_pal()(10), 'gray'),
                       breaks = c(names(table(pred.scRef)))) + 
    theme_bw() + 
    theme(axis.text = element_text(size = 9),
          panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 11))


