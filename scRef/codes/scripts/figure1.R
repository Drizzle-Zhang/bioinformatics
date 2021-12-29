library(Seurat)
library(ggplot2)
library(SeuratData)
data("pbmc3k")

pbmc3k.tag <- pbmc3k$seurat_annotations
table(pbmc3k.tag)
new.tag <- rep('0', length(pbmc3k.tag))
# new.tag[pbmc3k.tag %in% c('Naive CD4 T', 'Memory CD4 T', 'CD8 T')] <- 'Cell 1'
# new.tag[pbmc3k.tag %in% c('CD14+ Mono', 'FCGR3A+ Mono')] <- 'Cell 2'
# new.tag[pbmc3k.tag %in% c('NK')] <- 'Cell 3'
# new.tag[pbmc3k.tag %in% c('Platelet', 'DC', 'B')] <- 'Unassigned'
new.tag[pbmc3k.tag %in% c('Naive CD4 T', 'Memory CD4 T', 'CD8 T')] <- 'Cell 1'
new.tag[pbmc3k.tag %in% c('CD14+ Mono', 'FCGR3A+ Mono')] <- 'Cell 2'
new.tag[pbmc3k.tag %in% c('NK')] <- 'Cell 3'
new.tag[pbmc3k.tag %in% c('B')] <- 'Cell 4'
new.tag[pbmc3k.tag %in% c('Platelet', 'DC')] <- 'Unassigned'


pbmc <- CreateSeuratObject(counts = pbmc3k@assays$RNA@counts[, new.tag != '0'])
pbmc@meta.data$new.tag <- new.tag[new.tag != '0']

# add unassigned
cell_ids <- colnames(pbmc@assays$RNA@counts)
cell_ids_sample <- sample(cell_ids, 200)
unknown.tag <- new.tag[new.tag != '0']
unknown.tag[(cell_ids %in% cell_ids_sample) & (unknown.tag %in% c('Cell 1', 'Cell 2'))] <- 'Unassigned'
pbmc@meta.data$unknown.tag <- unknown.tag

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
plot.1 <- 
    DimPlot(pbmc, reduction = "umap", label = F, pt.size = 0.15, group.by = 'new.tag') + 
    labs(x = 'Dim-1', y = 'Dim-2') + 
    scale_color_manual(values = c(hue_pal()(4), 'DimGray'),
                       breaks = c(names(table(pbmc$new.tag)))) + 
    theme_bw() + 
    theme(axis.text = element_text(size = 9),
          panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 11))
plot.2 <- 
    DimPlot(pbmc, reduction = "umap", label = F, pt.size = 0.15, group.by = 'unknown.tag') + 
    labs(x = 'Dim-1', y = 'Dim-2') + 
    scale_color_manual(values = c(hue_pal()(4), 'DimGray'),
                       breaks = c(names(table(pbmc$unknown.tag)))) + 
    theme_bw() + 
    theme(axis.text = element_text(size = 9),
          panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 11))
ggsave(plot = plot.1, path = '/home/zy/scRef/figure', filename = 'scatter_4.png',
       units = 'cm', height = 10, width = 14)
ggsave(plot = plot.2, path = '/home/zy/scRef/figure', filename = 'scatter_unassigned_4.png',
       units = 'cm', height = 10, width = 14)
