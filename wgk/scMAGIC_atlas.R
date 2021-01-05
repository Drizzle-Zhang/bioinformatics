library(Seurat)

# file.ear <- '/home/disk/drizzle/wgk/data/ear.filtered.Rdata'
# seurat.ear <- readRDS(file.ear)
# exp_sc_mat <- seurat.ear@assays$RNA@counts

file.seurat <- '/home/disk/drizzle/wgk/data/ear.seurat.first.Rdata'
seurat.first <- readRDS(file.seurat)
exp_sc_mat <- seurat.first@assays$RNA@counts

source('/home/zy/my_git/scRef/main/scRef.v20.R')
setwd('~/my_git/scRef')
df.atlas <- .imoprt_outgroup('HCA', normalization = F)
result.scref <- SCREF(exp_sc_mat, df.atlas,
                      type_ref = 'sum-counts', use.RUVseq = T, out.group = 'HCA',
                      cluster.speed = T, cluster.cell = 10,
                      min_cell = 10, CPU = 4)
pred.scRef <- result.scref$final.out
df.tags <- result.scref$combine.out
df.tags$scRef.tag <- df.tags$scRef.tag.12
df.tags[df.tags$log10Pval < 20, 'scRef.tag'] <- 'Unassigned'

seurat.first$scref <- df.tags[colnames(seurat.first@assays$RNA@counts),'scRef.tag']
DimPlot(seurat.first, group.by = c("scref"), label = T)


result.scref.2 <- SCREF(exp_sc_mat, df.atlas,
                        type_ref = 'sum-counts', use.RUVseq = T, out.group = 'HCA',
                        cluster.speed = T, cluster.cell = 5,
                        GMM.floor_cutoff = 3, GMM.ceiling_cutoff = 20,
                        min_cell = 5, CPU = 4)
pred.scRef.2 <- result.scref.2$final.out
file.2 <- '/home/disk/drizzle/wgk/data/scMAGIC.2.Rdata'
saveRDS(result.scref.2, file.2)
result.scref.2 <- readRDS(file.2)
df.tags <- result.scref.2$combine.out
df.tags$scRef.tag <- df.tags$scRef.tag.12
df.tags[df.tags$log10Pval < 20, 'scRef.tag'] <- 'Unassigned'

# seurat.first$scref <- df.tags[colnames(seurat.first@assays$RNA@counts),'scRef.tag']
seurat.first$scref <- df.tags[colnames(seurat.first@assays$RNA@counts),'scRef.tag.12']
DimPlot(seurat.first, group.by = c("scref"), label = T)

result.scref.3 <- SCREF(exp_sc_mat, df.atlas,
                        type_ref = 'sum-counts', use.RUVseq = T, out.group = 'HCA',
                        method1 = 'spearman',
                        cluster.speed = T, cluster.cell = 5,
                        GMM.floor_cutoff = 4, GMM.ceiling_cutoff = 20,
                        min_cell = 5, CPU = 4)
pred.scRef.3 <- result.scref.3$final.out
file.3 <- '/home/disk/drizzle/wgk/data/scMAGIC.3.Rdata'
saveRDS(result.scref.3, file.3)
result.scref.3 <- readRDS(file.3)
seurat.first$scref <- pred.scRef.3[colnames(seurat.first@assays$RNA@counts),'scRef.tag']
DimPlot(seurat.first, group.by = c("scref"), label = T)

result.scref.4 <- SCREF(exp_sc_mat, df.atlas,
                        type_ref = 'sum-counts', use.RUVseq = T, out.group = 'HCA',
                        cluster.speed = T, cluster.cell = 5,
                        GMM.floor_cutoff = 3, GMM.ceiling_cutoff = 20,
                        min_cell = 5, CPU = 4)
pred.scRef.4 <- result.scref.4$final.out
file.4 <- '/home/disk/drizzle/wgk/data/scMAGIC.4.Rdata'
saveRDS(result.scref.4, file.4)
df.tags <- result.scref.2$combine.out
df.tags$scRef.tag <- df.tags$scRef.tag.12
df.tags[df.tags$log10Pval < 20, 'scRef.tag'] <- 'Unassigned'

seurat.first$scref <- pred.scRef.4[colnames(seurat.first@assays$RNA@counts),'scRef.tag']
DimPlot(seurat.first, group.by = c("scref"), label = T)


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


