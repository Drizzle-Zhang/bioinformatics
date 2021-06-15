setwd('/home/zy/my_git/bioinformatics/wgk')
.libPaths('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)
library(ggplot2)
require("RColorBrewer")

# file.all <- '/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/PBMC_harmony.RDS'
# seurat.all <- readRDS(file.all)

path.data <- '/home/disk/drizzle/wgk/data/AllSample_2_merge/'
file.merge_2 <- paste0(path.data, 'seurat_celltype.Rdata')
seurat.all.filter <- readRDS(file.merge_2)
seurat.chon <- subset(seurat.all.filter, subset = celltype %in% 
                        c('Chondrocyte1', 'Chondrocyte2', 'Chondral stem cell', 
                          'Transitional chondrocyte'))
# seurat.chon <- NormalizeData(seurat.chon)
# seurat.chon <- FindVariableFeatures(seurat.chon, nfeatures = 3000)
# seurat.chon <- ScaleData(seurat.chon, split.by = "batch")
# seurat.chon <- RunPCA(seurat.chon, verbose = F, npcs = 2)

path.lineage <- paste0(path.data, 'chon_lineage/')
if (!dir.exists(path.lineage)) {
    dir.create(path.lineage)
}


seurat.chon$celltype.abbr <- factor(seurat.chon$celltype.abbr, 
                                    levels = c('CSC', 'TC', 'Chondrocyte1', 'Chondrocyte2'))
plot.lineage <- 
    DimPlot(seurat.chon, reduction = 'pca', dims = c(2, 1),
        group.by = "celltype.abbr", label = T,
        pt.size = 0.5, label.size = 6,
        cols = c("#33A02C","#B2DF8A" ,"#A6CEE3","#1F78B4")) + 
    ylim(-12, 15) + labs(x = 'PC1', y = 'PC2')
ggsave(plot.lineage, path = path.lineage, 
       filename = 'chon_lineage.png',
       height = 15, width = 25, units = 'cm')

# DimPlot(seurat.all.filter, group.by = "celltype", label = T, reduction = 'pca')
# DimPlot(seurat.chon, group.by = "RNA_snn_res.3", label = T, reduction = 'pca')
# DimPlot(subset(seurat.chon, subset = RNA_snn_res.2 %in% c(2, 4, 5, 14, 26)), 
#         group.by = "RNA_snn_res.2", label = T, reduction = 'pca')
DimPlot(seurat.chon, group.by = "celltype.abbr", label = T, reduction = 'harmony')
DimPlot(seurat.chon, group.by = "RNA_snn_res.2", label = T, reduction = 'harmony')

# high correlated genes
PC2 <- seurat.chon@reductions$pca@cell.embeddings[, 'PC_2']
Harmony2 <- seurat.chon@reductions$harmony@cell.embeddings[, 'harmony_2']
mat.gene <- seurat.chon@assays$RNA@data
corr.PC.exp <- function(mat.gene, PC2, gene) {
    vec.exp <- mat.gene[gene,]
    corr <- cor(PC2, vec.exp, method = 'spearman')
    return(data.frame(gene = gene, corr = corr))
}
library(foreach)
library(doParallel)
registerDoParallel(cores = 10)
highvar.genes <- VariableFeatures(seurat.chon)
df.corr <- foreach(gene = highvar.genes, .combine = rbind) %dopar% corr.PC.exp(mat.gene, PC2, gene)
rownames(df.corr) <- highvar.genes
df.corr['COL2A1',]

# mat.gene <- seurat.chon@assays$RNA@data[c('COL2A1', 'VIT', 'CYTL1', 'FRZB', 'CTGF',
#                                           'HES1', 'EGR1'),]
df.pc.gene <- data.frame(t(as.matrix(mat.gene)))
df.pc.gene$PC2 <- PC2
df.pc.gene$Harmony2 <- Harmony2
df.pc.gene$celltype <- seurat.chon$celltype.abbr
df.pc.gene$status <- seurat.chon$type

ggplot(data = df.pc.gene, aes(x = PC2, y = EGR1)) + 
    geom_point(aes(color = celltype), size = 1) + 
    scale_color_manual(labels = c('CSC','TC','Chondrocyte1','Chondrocyte2'),
                       values = c("#33A02C","#B2DF8A" ,"#A6CEE3","#1F78B4")) + 
    # geom_smooth(method = 'lm', formula = y~poly(x,5), color = '#696969') + 
    geom_smooth(color = '#696969') + 
    labs(x = 'PC1') + 
    theme(panel.background=element_rect(fill='transparent', color='black',size = 1),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'none')

ggplot(data = df.pc.gene, aes(x = Harmony2, y = EGR1)) + 
    geom_point(aes(color = celltype), size = 1) + 
    scale_color_manual(labels = c('CSC','TC','Chondrocyte1','Chondrocyte2'),
                       values = c("#33A02C","#B2DF8A" ,"#A6CEE3","#1F78B4")) + 
    # geom_smooth(method = 'lm', formula = y~poly(x,5), color = '#696969') + 
    geom_smooth(color = '#696969') + 
    labs(x = 'PC1') + 
    theme(panel.background=element_rect(fill='transparent', color='black',size = 1),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'none')


ggplot(data = df.pc.gene, aes(x = PC2, y = COL2A1, linetype = status)) + 
    geom_point(aes(color = celltype), size = 1) + 
    scale_color_manual(labels = c('CSC','TC','Chondrocyte1','Chondrocyte2'),
                       values = c("#33A02C","#B2DF8A" ,"#A6CEE3","#1F78B4")) + 
    # geom_smooth(method = 'lm', formula = y~poly(x,5), color = '#696969') + 
    geom_smooth(color = '#696969') + 
    theme(panel.background=element_rect(fill='transparent', color='black',size = 1),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'none')
ggplot(data = df.pc.gene, aes(x = PC2, y = VIT, color = celltype)) + 
    geom_point(size = 0.5)

ggplot(data = df.pc.gene, aes(x = PC2, y = CYTL1, linetype = status)) + 
    geom_point(aes(color = celltype), size = 1) + 
    scale_color_manual(labels = c('CSC','TC','Chondrocyte1','Chondrocyte2'),
                       values = c("#33A02C","#B2DF8A" ,"#A6CEE3","#1F78B4")) + 
    # geom_smooth(method = 'lm', formula = y~poly(x,3), color = '#696969') +
    geom_smooth(color = '#696969')
