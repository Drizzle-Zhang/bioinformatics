setwd('/home/zy/my_git/bioinformatics/wgk')
.libPaths('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)
library(ggplot2)
require("RColorBrewer")

file.all <- '/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat.all.filter_harmony.RDS'
seurat.all <- readRDS(file.all)

path.data <- '/home/disk/drizzle/wgk/data/AllSample_3_merge/'
if (!dir.exists(path.data)) {
    dir.create(path.data)
}

DimPlot(seurat.all, group.by = "batch")
status <- rep('0', length(seurat.all$batch))
status[seurat.all$batch %in% c('C1', 'C2', 'C3', 'C4', 'C5')] <- 'Normal'
status[seurat.all$batch %in% c('M1', 'M2', 'M3')] <- 'Microtia'
seurat.all$status <- status
DimPlot(seurat.all, group.by = "status")

seurat.all <- FindNeighbors(seurat.all, reduction = "harmony", dims = 1:50)
seurat.all <- FindClusters(seurat.all, resolution = 3)
DimPlot(seurat.all, group.by = "RNA_snn_res.3", label = T)


FeaturePlot(seurat.all, features = c('ID3', 'HES1', 'COL1A1', 'CYTL1'))
FeaturePlot(seurat.all, features = c('DCN', 'STC1', 'COL9A3', 'FRZB'))

# unknown 23 28
FeaturePlot(seurat.all, features = c('IFIT1', 'OAS3', 'XAF1', 'OAS1'))

# neuron 27
FeaturePlot(seurat.all, features = c('PLP1', 'SLITRK6', 'ANK3', 'CIT'))

# commu stromal 10
FeaturePlot(seurat.all, features = c('CFD', 'APOE', 'VCAN', 'CSF3'))

# matrix stromal 3
FeaturePlot(seurat.all, features = c('MMP10', 'COCH', 'OGN', 'COMP'))

# stromal stem cell 22
FeaturePlot(seurat.all, features = c('HES1', 'ID3', 'COL1A1', 'LUM'))

# Chondral stem cell 9 21
FeaturePlot(seurat.all, features = c('HES1', 'ID3', 'CYTL1', 'COL2A1'))

# Transitional chondrocyte 13 17
FeaturePlot(seurat.all, features = c('HES1', 'ID3', 'CYTL1', 'COL2A1'))

# Chondrocyte1  1 2 5 8 12 14 15 26
FeaturePlot(seurat.all, features = c('COL1A1', 'COL1A2', 'COL2A1'), ncol=3)

# Chondrocyte2  0 4 6 7 11 16
FeaturePlot(seurat.all, features = c('STC1'))

# filter
seurat.all.filter <- subset(seurat.all, 
                            subset = RNA_snn_res.3 %in% 
                                setdiff(unique(seurat.all$RNA_snn_res.3), c(26, 32, 38, 37)))
DimPlot(seurat.all.filter, group.by = "RNA_snn_res.3", label = T)

cluster_3 <- seurat.all.filter$RNA_snn_res.3
celltypes <- rep('_', length(cluster_3))
celltypes[cluster_3 %in% c(6)] <- 'Stromal cell1'
celltypes[cluster_3 %in% c(7, 31)] <- 'Stromal cell2'
celltypes[cluster_3 %in% c(29, 34)] <- 'Stromal stem cell'
celltypes[cluster_3 %in% c(17, 23, 27)] <- 'Chondral stem cell'
celltypes[cluster_3 %in% c(12, 16, 25, 9)] <- 'Transitional chondrocyte'
celltypes[cluster_3 %in% c(0, 1, 3, 5, 10, 11, 14, 18, 20, 22, 39)] <- 'Chondrocyte1'
celltypes[cluster_3 %in% c(2, 4, 8, 13, 15, 19, 21)] <- 'Chondrocyte2'
celltypes[cluster_3 %in% c(28, 33, 24)] <- 'Perivascular cell'
celltypes[cluster_3 %in% c(30, 35)] <- 'Endothelial cell'
celltypes[cluster_3 %in% c(36)] <- 'Immune cell'
seurat.all.filter$celltype <- celltypes
DimPlot(seurat.all.filter, group.by = "celltype", label = T)


seurat.all.filter$type[seurat.all.filter$type=='Control'] <- 'Normal'
seurat.all.filter$type <- factor(seurat.all.filter$type,levels = c('Normal','Microtia'))
seurat.all.filter$celltype.abbr <- seurat.all.filter$celltype
seurat.all.filter$celltype.abbr[seurat.all.filter$celltype=="Chondral stem cell"] <- 'CSC'
seurat.all.filter$celltype.abbr[seurat.all.filter$celltype == "Transitional chondrocyte"] <- 'TC'
seurat.all.filter$celltype.abbr[seurat.all.filter$celltype =="Stromal stem cell"] <- 'SSC'
seurat.all.filter$celltype.abbr[seurat.all.filter$celltype == "Stromal cell1"] <- 'SC1'
seurat.all.filter$celltype.abbr[seurat.all.filter$celltype == "Stromal cell2"] <- 'SC2'
seurat.all.filter$celltype.abbr[seurat.all.filter$celltype== "Immune cell"] <- 'IC'
seurat.all.filter$celltype.abbr[seurat.all.filter$celltype == "Endothelial cell"] <- 'EC'
seurat.all.filter$celltype.abbr[seurat.all.filter$celltype == "Perivascular cell"] <- 'PVC'

require('RColorBrewer')
library("scales")
library(ggsci)
# show_col(pal_aaas()(10))
# mypalette<- c(brewer.pal(11,"Paired"))
CT <- c('CSC','TC','Chondrocyte1','Chondrocyte2','SSC','SC1','SC2','IC','PVC','EC')
# Color <- c("#6A3D9A","#1F78B4","#B2DF8A","#FB9A99",
#            "#FF7F00","#E31A1C","#33A02C",
#            "#FDBF6F","#A6CEE3","#CAB2D6")
# Color <- c("#9932CC","#4169E1","#CD853F","#20854EFF",
#            "#FF7F00","#C71585","#33A02C",
#            "#0000CD","#40E0D0","#FF0000")
Color <- c("#33A02C","#B2DF8A" ,"#1F78B4","#A6CEE3" ,
           "#FB9A99", "#FDBF6F","#CAB2D6", 
           "#E31A1C" ,"#FF7F00","#6A3D9A")
names(Color) <- CT

seurat.all.filter$color <- Color[seurat.all.filter$celltype.abbr]
seurat.all.filter$celltype.abbr <- factor(seurat.all.filter$celltype.abbr,levels = CT)
Idents(seurat.all.filter) <- seurat.all.filter$celltype.abbr
plot.celltype <- 
    DimPlot(seurat.all.filter, group.by='celltype.abbr', label=T, repel = T, 
            pt.size = 0.1, label.size = 6, cols = Color) +
    theme(panel.background=element_rect(fill='transparent', color='black',size = 1),
          legend.key=element_rect(fill='transparent', color='transparent'), 
          legend.text = element_text(size=18), 
          axis.title = element_text(size = 16),
          axis.text = element_blank(), 
          axis.ticks = element_blank()) + 
    guides(colour = guide_legend(override.aes = list(size=10))) + 
    labs(x = 'UMAP1', y = 'UMAP2')
ggsave(plot = plot.celltype, path = path.data, 
       filename = 'celltype.png',
       height = 16, width = 23, units = 'cm')

# saveRDS(seurat.all.filter,'/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')

file.merge_3 <- paste0(path.data, 'seurat_celltype.Rdata')
saveRDS(seurat.all.filter, file.merge_3)
seurat.all.filter <- readRDS(file.merge_3)


# violin plot
library(ggplot2)
marker.genes <- c('APOE', 'VCAN', 'OGN', 'COCH', 
                  'COL1A1', 'LUM', 'HES1', 'ID3', 'EGR1', 
                  'CYTL1', 'COL2A1', 'COL9A2', 'ACAN',
                  'ACTA2', 'PDGFRB', 'CDH5', 'IL1B')
# marker.genes <- c('CFD', 'APOE', 'VCAN', 'MMP10', 'COCH', 
#                   'ASPN', 'OGN')
df.gene <- data.frame(stringsAsFactors = F)
for (gene in marker.genes) {
    df.sub <- data.frame(expvalue = seurat.all.filter@assays$RNA@data[gene,],
                         gene = rep(gene, ncol(seurat.all.filter@assays$RNA@data)),
                         celltype = seurat.all.filter$celltype.abbr)
    df.gene <- rbind(df.gene, df.sub)
}
df.plot <- df.gene
df.plot$gene <- factor(df.gene$gene, levels = marker.genes)
df.plot$celltype <- factor(df.gene$celltype, 
                           levels = c('SC1', 'SC2', 'SSC', 
                                      'CSC', 'TC',
                                      'Chondrocyte1', 'Chondrocyte2','PVC',
                                      'EC', 'IC'))
plot.vln <- 
    ggplot(data = df.plot, aes(x = gene, y = expvalue, color = celltype, fill = celltype)) + 
    geom_violin(trim = T, scale = 'width') + 
    scale_color_manual(labels = c('CSC','TC','Chondrocyte1','Chondrocyte2',
                                  'SSC','SC1','SC2','IC','PVC','EC'),
                       values = Color) + 
    scale_fill_manual(labels = c('CSC','TC','Chondrocyte1','Chondrocyte2',
                                  'SSC','SC1','SC2','IC','PVC','EC'),
                       values = Color) + 
    facet_grid( ~ celltype) + 
    theme_classic() + coord_flip() +
    stat_summary(fun= mean, geom = "point",
                 shape = 23, size = 2, color = "black") + 
    labs(x = 'Gene', y = 'Expression Level') + 
    theme(axis.text.y = element_text(size = 13, face = 'bold'), 
          axis.title = element_text(size = 15, face = 'bold'), 
          legend.position = 'none')

ggsave(plot = plot.vln, path = path.data, 
       filename = 'Vln.png',
       height = 15, width = 35, units = 'cm')


# cell marker
clusters <- unique(seurat.all.filter$celltype.abbr)
list.marker.all <- list()
for (cluster in clusters) {
    sub.markers <- FindMarkers(seurat.all.filter, ident.1 = cluster, group.by = 'celltype.abbr',
                               logfc.threshold = 0.25, min.diff.pct = 0.05, only.pos = T)
    list.marker.all[[cluster]] <- sub.markers
}
file.marker_3 <- paste0(path.data, 'marker_genes.Rdata')
saveRDS(list.marker.all, file.marker_3)
list.marker.all <- readRDS(file.marker_3)

# sort.cells <- c('SC1', 'SC2', 'SSC', 
#                 'CSC', 'TC', 'Chondrocyte1', 'Chondrocyte2',
#                 'PVC', 'EC', 'IC')
sort.cells <- c('IC', 'EC', 'PVC', 
                'SC1', 'SC2', 'SSC', 
                'CSC', 'TC', 'Chondrocyte1', 'Chondrocyte2')
color.cell <- c("#E31A1C","#6A3D9A","#FF7F00",
              "#CAB2D6","#FDBF6F","#FB9A99",
              "#33A02C","#B2DF8A" ,"#A6CEE3","#1F78B4")
sel.genes <- c()
for (cell in sort.cells) {
    sub.markers <- list.marker.all[[cell]]
    sub.markers$diff.pct <- sub.markers$pct.1 - sub.markers$pct.2
    sub.markers <- sub.markers[sub.markers$p_val_adj < 0.01 & sub.markers$avg_logFC > 0.5 & 
                                   sub.markers$diff.pct > 0.1, ]
    sub.markers <- sub.markers[order(sub.markers$avg_logFC, decreasing = T),]
    sel.genes <- c(sel.genes, rownames(sub.markers)[1:min(nrow(sub.markers), 50)])
}
sel.genes <- unique(sel.genes)

tmp <- AverageExpression(seurat.all.filter, return.seurat = TRUE)
tmp@active.ident <- factor(tmp@active.ident, levels = sort.cells)
# mat.scale <- scale(t(as.matrix(tmp@assays$RNA@data[sel.genes,])))
# tmp@assays$RNA@scale.data <- t(mat.scale)
# pdf(paste0(path.data, "DoHeatmap_markers_CT_1.pdf"),width = 10,height = 20)
plot.heatmap <- 
    DoHeatmap(tmp, features = sel.genes,draw.lines = FALSE, 
          label = T, group.colors = color.cell, size = 6,
          slot = 'scale.data') +
    guides(color = F) + 
    labs(fill = 'Scaled Expression') + 
    theme(axis.text = element_blank(),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 15),
          legend.position = 'bottom', legend.direction = 'horizontal',
          plot.margin = unit(c(0.5,2,0,0), "cm"))+
    scale_fill_gradientn( colors = c("navy", "white", "firebrick3"))
# dev.off()
ggsave(plot = plot.heatmap, path = path.data, 
       filename = 'DoHeatmap_markers.png',
       height = 30, width = 12, units = 'cm')

# # prop
# seurat.4 <- subset(seurat.all, subset = RNA_snn_res.0.6 == 4)
# table(seurat.4$batch)/table(seurat.all$batch)
# df.plot <- data.frame(sample = names(table(seurat.4$batch)), 
#                       prop = as.numeric(table(seurat.4$batch)/table(seurat.all$batch)))
# ggplot(df.plot, aes(x = sample, y = prop)) + geom_bar(stat = 'identity')
# 
# seurat.5 <- subset(seurat.all, subset = RNA_snn_res.0.6 == 5)
# table(seurat.5$batch)
# seurat.9 <- subset(seurat.all, subset = RNA_snn_res.0.6 == 9)
# table(seurat.9$batch)
# 
# seurat.11 <- subset(seurat.all, subset = RNA_snn_res.2 == 11)
# DimPlot(seurat.11, group.by = 'RNA_snn_res.2')
# seurat.0 <- subset(seurat.all, subset = RNA_snn_res.2 == 0)
# DimPlot(seurat.0, group.by = 'RNA_snn_res.2')
# 
# seurat.M3 <- subset(seurat.all, subset = batch == 'M3')
# FeaturePlot(seurat.M3, features = c('HES1', 'ID3', 'COL1A1', 'CYTL1'))
# 
# seurat.M1 <- subset(seurat.all, subset = batch == 'M1')
# FeaturePlot(seurat.M1, features = c('HES1', 'ID3', 'COL1A1', 'CYTL1'))
# FeaturePlot(seurat.M1, features = c('TNF'))
# 
# seurat.C4 <- subset(seurat.all, subset = batch == 'C4')
# FeaturePlot(seurat.C4, features = c('HES1', 'ID3', 'COL1A1', 'CYTL1'))
# 
# seurat.C6 <- subset(seurat.all, subset = batch == 'C6')
# FeaturePlot(seurat.C6, features = c('HES1', 'ID3', 'COL1A1', 'CYTL1'))
# FeaturePlot(seurat.C6, features = c('TNF'))
# 
# seurat.C1 <- subset(seurat.all, subset = batch == 'C1')
# FeaturePlot(seurat.C1, features = c('TNF', 'TNFSF10'))
# 
# (table(seurat.4$batch) + table(seurat.5$batch))/table(seurat.all$batch)
# 
# FeaturePlot(seurat.first, features = c('FRZB', 'CTGF', 
#                                        'SERPINA1', "SCRG1", 'COL9A3', 'FGFBP2'), ncol = 3)

