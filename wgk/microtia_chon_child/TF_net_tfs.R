.libPaths('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)
library(ggplot2)
library(SCENIC)
require("RColorBrewer")
library(maSigPro)


path.data <- '/home/disk/drizzle/wgk/data/AllSample_2_merge/'
path.lineage <- paste0(path.data, 'chon_lineage/')
file.chon <- paste0(path.lineage, 'seurat_celltype.Rdata')
seurat.chon <- readRDS(file.chon)
seurat.child <- subset(seurat.chon, subset = batch %in% c('C4', 'C6', 'M1', 'M2', 'M3'))
# seurat.child <- NormalizeData(seurat.child)
# seurat.child <- FindVariableFeatures(seurat.child, nfeatures = 5000)
# highvar.genes <- VariableFeatures(seurat.child)
# seurat.child <- ScaleData(seurat.child, split.by = "batch", 
#                           features = rownames(seurat.child@assays$RNA@counts))

path.M123 <- '/home/disk/drizzle/wgk/microtia_chon_child_M1M2M3/'
path.change <- paste0(path.M123, 'chon_lineage_1/')
if (!file.exists(path.change)) {
    dir.create(path.change)
}

# plot single gene
Harmony2 <- seurat.child@reductions$harmony@cell.embeddings[, 'harmony_2']
mat.gene <- seurat.child@assays$RNA@data
# AUC
regulonAUC <- readRDS(file='/home/yzj/JingMA_NEW/res/SCENIC_main/int/3.4_regulonAUC.Rds')
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), colnames(mat.gene)]
mat.auc <- as.matrix(regulonAUC@assays@data@listData$AUC)
# mat.gene.TF <- rbind(as.matrix(mat.gene), mat.auc)
df.pc.gene <- data.frame(t(rbind(as.matrix(mat.gene), mat.auc)), check.names = F)
df.pc.gene$Harmony2 <- Harmony2
df.pc.gene$celltype <- seurat.child$celltype
df.pc.gene$status <- seurat.child$type
df.pc.gene <- df.pc.gene[order(Harmony2, decreasing = F),]
df.pc.gene$idx <- 1:nrow(df.pc.gene)
colors <- c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
names(colors) <- c('CSC', 'TC', 'C1', 'C2')

gene <- 'DDIT3'
# gene <- 'SOX5 (218g)'
df.plot <- df.pc.gene[, c('idx', 'Harmony2', 'celltype', 'status', gene)]
names(df.plot) <- c('idx', 'Harmony2', 'celltype', 'status', 'gene')
p.gene <-
    ggplot(data = df.plot, aes(x = idx, 
                               linetype = status, 
                               y = gene)) + 
    geom_point(aes(color = celltype), size = 0.3) + 
    scale_color_manual(labels = c('CSC', 'TC', 'C1', 'C2'),
                       values = colors) + 
    # xlim(-30, 10) + 
    geom_smooth(color = '#696969') +
    # geom_smooth(color = '#696969', formula = y~poly(x, 3), method = lm) + 
    labs(x = '', y = '') + 
    theme(panel.background=element_rect(fill='transparent', color='black',size = 1),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'none') +
    theme(panel.background=element_rect(fill='transparent', color='black',size = 1)) +
    annotate('text', label = gene, x = 22000, y = max(df.plot$gene), 
             hjust = 1, vjust = 1, size = 7)

# TF change
# gene <- 'SOX5 (218g)'
vec.TF <- c('SOX5 (218g)', 'DBP (45g)', 'SOX8 (158g)',
            'EGR1 (264g)', 'EGR3 (53g)', 'KLF10 (16g)',
            'ATF3 (83g)', 'REL (834g)', 'BCL3 (162g)',
            'JUNB (158g)', 'CEBPB (258g)', 'CEBPD (199g)')
vec.TF <- c('SOX5 (218g)', 'DBP (45g)', 'SOX8 (158g)',
            'EGR1 (264g)', 'EGR3 (53g)', 'KLF10 (16g)',
            'ATF3 (83g)')
for (i in 1:length(vec.TF)) {
    TF <- vec.TF[i]
    df.plot <- df.pc.gene[, c('idx', 'Harmony2', 'celltype', 'status', TF)]
    names(df.plot) <- c('idx', 'Harmony2', 'celltype', 'status', 'TF')
    p1 <- ggplot(data = df.plot, aes(x = Harmony2, 
                                     linetype = status, 
                                     y = TF)) + 
        geom_point(aes(color = celltype), size = 0.1) + 
        scale_color_manual(labels = c('CSC', 'TC', 'C1', 'C2'),
                           values = colors) + 
        xlim(-30, 10) +
        geom_smooth(color = '#696969') +
        # geom_smooth(color = '#696969', formula = y~poly(x, 3), method = lm) + 
        labs(x = '', y = '') + 
        theme(panel.background=element_rect(fill='transparent', color='black',size = 1),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.position = 'none') +
        theme(panel.background=element_rect(fill='transparent', color='black',size = 1)) +
        annotate('text', label = TF, x = 9, y = max(df.plot$TF), 
                 hjust = 1, vjust = 1, size = 7)
    if (i == 1) {
        p <- p1
    } else {
        p <- p / p1
    }
}
ggsave(plot = p, path = path.change, 
       filename = 'TF_AUC.png',
       height = 35, width = 10, units = 'cm')

vec.TF.exp <- c('SOX5', 'DBP', 'SOX8',
                'EGR1', 'EGR3', 'KLF10',
                'ATF3')
for (i in 1:length(vec.TF.exp)) {
    TF <- vec.TF.exp[i]
    df.plot <- df.pc.gene[, c('idx', 'Harmony2', 'celltype', 'status', TF)]
    names(df.plot) <- c('idx', 'Harmony2', 'celltype', 'status', 'TF')
    p1 <- ggplot(data = df.plot, aes(x = Harmony2, 
                                     linetype = status, 
                                     y = TF)) + 
        geom_point(aes(color = celltype), size = 0.1) + 
        scale_color_manual(labels = c('CSC', 'TC', 'C1', 'C2'),
                           values = colors) + 
        xlim(-30, 10) +
        geom_smooth(color = '#696969') +
        # geom_smooth(color = '#696969', formula = y~poly(x, 3), method = lm) + 
        labs(x = '', y = '') + 
        theme(panel.background=element_rect(fill='transparent', color='black',size = 1),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.position = 'none') +
        theme(panel.background=element_rect(fill='transparent', color='black',size = 1)) +
        annotate('text', label = TF, x = 9, y = max(df.plot$TF), 
                 hjust = 1, vjust = 1, size = 7)
    if (i == 1) {
        p <- p1
    } else {
        p <- p / p1
    }
}
ggsave(plot = p, path = path.change, 
       filename = 'TF_EXP.png',
       height = 35, width = 10, units = 'cm')

