setwd('/home/zy/my_git/bioinformatics/wgk')
.libPaths('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)
library(ggplot2)
require("RColorBrewer")

file.all <- '/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/PBMC_harmony.RDS'
seurat.all <- readRDS(file.all)

path.data <- '/home/disk/drizzle/wgk/data/AllSample_2_merge/'

DimPlot(seurat.all, group.by = "batch")
status <- rep('0', length(seurat.all$batch))
status[seurat.all$batch %in% c('C1', 'C2', 'C3', 'C4', 'C5')] <- 'Normal'
status[seurat.all$batch %in% c('M1', 'M2', 'M3')] <- 'Microtia'
seurat.all$status <- status
DimPlot(seurat.all, group.by = "status")

seurat.all <- FindNeighbors(seurat.all, reduction = "pca", dims = 1:100)
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
                            subset = RNA_snn_res.2 %in% 
                                setdiff(unique(seurat.all$RNA_snn_res.2), c(23, 27, 28)))
seurat.all.filter <- FindNeighbors(seurat.all.filter, reduction = "pca", dims = 1:100)
seurat.all.filter <- FindClusters(seurat.all.filter, resolution = 1.5)
DimPlot(seurat.all.filter, group.by = "RNA_snn_res.2", label = T)

cluster_3 <- seurat.all$RNA_snn_res.3
celltypes <- rep('_', length(cluster_3))
celltypes[cluster_3 %in% c(22, 25)] <- 'Stromal cell1'
celltypes[cluster_3 %in% c(20, 21, 29)] <- 'Stromal cell2'
celltypes[cluster_3 %in% c(27)] <- 'Stromal stem cell'
celltypes[cluster_3 %in% c(18, 23, 33, 35)] <- 'Chondral stem cell'
celltypes[cluster_3 %in% c(12, 16, 25, 9, 18)] <- 'Transitional chondrocyte'
celltypes[cluster_3 %in% c(10, 16, 1, 5, 15, 12, 17, 4, 6, 13, 8)] <- 'Chondrocyte1'
celltypes[cluster_3 %in% c(0, 3, 7, 11, 14, 9, 19, 26, 30, 38)] <- 'Chondrocyte2'
celltypes[cluster_3 %in% c(28, 36, 24)] <- 'Perivascular cell'
celltypes[cluster_3 %in% c(32, 37, 41)] <- 'Endothelial cell'
celltypes[cluster_3 %in% c(39)] <- 'Immune cell'
seurat.all$celltype_3 <- celltypes
DimPlot(seurat.all, group.by = "celltype_3", label = T)


# cell type
cluster_2 <- seurat.all.filter$RNA_snn_res.2
celltypes <- rep('_', length(cluster_2))
celltypes[cluster_2 %in% c(10)] <- 'Stromal cell1'
celltypes[cluster_2 %in% c(3)] <- 'Stromal cell2'
celltypes[cluster_2 %in% c(22)] <- 'Stromal stem cell'
celltypes[cluster_2 %in% c(9, 21)] <- 'Chondral stem cell'
celltypes[cluster_2 %in% c(13, 17, 14, 5)] <- 'Transitional chondrocyte'
celltypes[cluster_2 %in% c(1, 2, 5, 8, 12, 14, 15, 26)] <- 'Chondrocyte1'
celltypes[cluster_2 %in% c(0, 4, 6, 7, 11, 16)] <- 'Chondrocyte2'
celltypes[cluster_2 %in% c(19, 20, 24)] <- 'Perivascular cell'
celltypes[cluster_2 %in% c(18)] <- 'Endothelial cell'
celltypes[cluster_2 %in% c(25)] <- 'Immune cell'
table(seurat.all.filter$batch, celltypes)/as.vector(table(seurat.all.filter$batch))
seurat.all.filter$celltype <- celltypes
DimPlot(seurat.all.filter, group.by = "celltype", label = T)

file.merge_2 <- paste0(path.data, 'seurat_celltype.Rdata')
saveRDS(seurat.all.filter, file.merge_2)
seurat.all.filter <- readRDS(file.merge_2)

# cell marker
clusters <- unique(seurat.all.filter$celltype)
list.marker.all <- list()
for (cluster in clusters) {
    sub.markers <- FindMarkers(seurat.all.filter, ident.1 = cluster, group.by = 'celltype',
                               logfc.threshold = 0.25, min.diff.pct = 0.05, only.pos = T)
    list.marker.all[[cluster]] <- sub.markers
}


# surface marker
file.surface <- '/home/disk/drizzle/wgk/public_data/surface_marker_wlab.txt'
df.surface <- read.delim2(file.surface)
surface.genes <- unique(toupper(df.surface$ENTREZ.gene.symbol))

# chon and stroma lineage marker genes
seurat.CS <- subset(seurat.all.filter, subset = 
                        celltype %in% c('Stromal cell1', 'Stromal cell2', 'Stromal stem cell',
                                        'Chondral stem cell', 'Transitional chondrocyte', 
                                        'Chondrocyte1', 'Chondrocyte2'))
clusters <- unique(seurat.CS$celltype)
list.marker <- list()
for (cluster in clusters) {
    sub.markers <- FindMarkers(seurat.CS, ident.1 = cluster, group.by = 'celltype',
                               logfc.threshold = 0.3, min.diff.pct = 0, only.pos = T)
    list.marker[[cluster]] <- sub.markers
}

View(list.marker$`Chondral stem cell`)
View(list.marker$`Stromal stem cell`)
intersect(rownames(list.marker$`Chondral stem cell`), surface.genes)
intersect(rownames(list.marker$`Stromal stem cell`), surface.genes)
genes.chon <- c("CD83", "TSPAN6", 'ENPP1', 'ITM2B', 'CD99', 'TSPAN4', 
                'SCARA3', 'BOC', 'PTPRZ1', 'MXRA8', 'FGFR3', 'AQP1',
                'CNTFR', 'CADM1', 'SLC29A1', 'LRP1', 'CD46', 'A2M', 'ITGA10')
genes.stroma <- c("GPNMB", "LRRC32", 'ABCA8', 'SCARA5', 'THY1', 'SLC2A3', 
                'MXRA8', 'PDGFRB', 'LEPR', 'SLC38A2', 'ANK2', 'PCDH9',
                'CLEC2B', 'GPRC5A', 'PRNP', 'SLC20A1', 'CERCAM', 'BOC', 'FAP',
                'PLP2', 'CD276', 'AQP1', 'ABCA1', 'ITM2B', 'BMPR2', 'TSPAN4',
                'SGCE', 'CD83', 'EMP1', 'SLC19A2', 'F3', 'TMEM2', 'FCGRT', 
                'LRP1', 'IL1R1')
genes.CS <- unique(c(genes.chon, genes.stroma))
celltypes <- unique(seurat.all.filter$celltype)
df.bubble <- data.frame(stringsAsFactors = F)
for (cell in celltypes) {
    sub.seurat <- subset(seurat.all.filter, subset = celltype == cell)
    sub.mat <- sub.seurat@assays$RNA@data
    for (gene in genes.CS) {
        vec.exp <- sub.mat[gene,]
        mean.exp <- mean(vec.exp)
        prop.exp <- sum(vec.exp != 0)/length(vec.exp)
        df.bubble <- rbind(df.bubble, 
                           data.frame(Cell = cell, Gene = gene,
                                      MeanExp = mean.exp, Prop = prop.exp))
    }
}
mat.exp <- reshape2::dcast(df.bubble, Cell ~ Gene, value.var = 'MeanExp')
row.names(mat.exp) <- mat.exp$Cell
mat.exp$Cell <- NULL
mat.exp.scale <- scale(mat.exp)
df.exp.scale <- reshape2::melt(mat.exp.scale)
names(df.exp.scale) <- c('Cell', 'Gene', 'ScaleExp')
df.bubble <- merge(df.bubble, df.exp.scale, by = c('Cell', 'Gene'))
df.bubble$Cell <- factor(df.bubble$Cell, 
                           levels = c('Stromal cell1', 'Stromal cell2', 'Stromal stem cell', 
                                      'Chondral stem cell', 'Transitional chondrocyte',
                                      'Chondrocyte1', 'Chondrocyte2','Perivascular cell',
                                      'Endothelial cell', 'Immune cell'))

plot.bubble <- 
    ggplot(data = df.bubble, 
           aes(x = Gene, y = Cell, size = Prop, color = ScaleExp)) + 
    geom_point(fill = 'cornsilk') + 
    scale_colour_gradientn(
        colours = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) + 
    # facet_grid( ~ Status) + 
    labs(x = 'Gene', y = 'Cell Type', color = 'Scaled expression', 
         size = 'Proportion') + 
    theme(panel.background = element_rect(color = 'gray',
                                          fill = 'transparent'),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(plot = plot.bubble, path = path.data, 
       filename = 'marker_surface.png',
       height = 15, width = 35, units = 'cm')

# no filter
# cell type
cluster_2 <- seurat.all$RNA_snn_res.2
celltypes <- as.character(cluster_2)
celltypes[cluster_2 %in% c(10)] <- 'Stromal cell1'
celltypes[cluster_2 %in% c(3)] <- 'Stromal cell2'
celltypes[cluster_2 %in% c(22)] <- 'Stromal stem cell'
celltypes[cluster_2 %in% c(9, 21)] <- 'Chondral stem cell'
celltypes[cluster_2 %in% c(13, 17)] <- 'Transitional chondrocyte'
celltypes[cluster_2 %in% c(1, 2, 5, 8, 12, 14, 15, 26)] <- 'Chondrocyte1'
celltypes[cluster_2 %in% c(0, 4, 6, 7, 11, 16)] <- 'Chondrocyte2'
celltypes[cluster_2 %in% c(19, 20, 24)] <- 'Perivascular cell'
celltypes[cluster_2 %in% c(18)] <- 'Endothelial cell'
celltypes[cluster_2 %in% c(25)] <- 'Immune cell'
# table(seurat.all$batch, celltypes)/as.vector(table(seurat.all$batch))
seurat.all$celltype <- celltypes
DimPlot(seurat.all, group.by = "celltype", label = T)

celltypes <- unique(seurat.all$celltype)
df.bubble <- data.frame(stringsAsFactors = F)
for (cell in celltypes) {
    sub.seurat <- subset(seurat.all, subset = celltype == cell)
    sub.mat <- sub.seurat@assays$RNA@data
    for (gene in genes.CS) {
        vec.exp <- sub.mat[gene,]
        mean.exp <- mean(vec.exp)
        prop.exp <- sum(vec.exp != 0)/length(vec.exp)
        df.bubble <- rbind(df.bubble, 
                           data.frame(Cell = cell, Gene = gene,
                                      MeanExp = mean.exp, Prop = prop.exp))
    }
}
mat.exp <- reshape2::dcast(df.bubble, Cell ~ Gene, value.var = 'MeanExp')
row.names(mat.exp) <- mat.exp$Cell
mat.exp$Cell <- NULL
mat.exp.scale <- scale(mat.exp)
df.exp.scale <- reshape2::melt(mat.exp.scale)
names(df.exp.scale) <- c('Cell', 'Gene', 'ScaleExp')
df.bubble <- merge(df.bubble, df.exp.scale, by = c('Cell', 'Gene'))
df.bubble$Cell <- factor(df.bubble$Cell, 
                         levels = c('Stromal cell1', 'Stromal cell2', 'Stromal stem cell', 
                                    'Chondral stem cell', 'Transitional chondrocyte',
                                    'Chondrocyte1', 'Chondrocyte2','Perivascular cell',
                                    'Endothelial cell', 'Immune cell', 
                                    '23', '27', '28'))

plot.bubble <- 
    ggplot(data = df.bubble, 
           aes(x = Gene, y = Cell, size = Prop, color = ScaleExp)) + 
    geom_point(fill = 'cornsilk') + 
    scale_colour_gradientn(
        colours = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) + 
    # facet_grid( ~ Status) + 
    labs(x = 'Gene', y = 'Cell Type', color = 'Scaled expression', 
         size = 'Proportion') + 
    theme(panel.background = element_rect(color = 'gray',
                                          fill = 'transparent'),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(plot = plot.bubble, path = path.data, 
       filename = 'marker_surface_no_filter.png',
       height = 18, width = 35, units = 'cm')

# sort
df.chon <- df.bubble[df.bubble$Cell == 'Chondral stem cell', ]
df.bubble$Gene <- factor(df.bubble$Gene, 
                         levels = df.chon[order(df.chon$ScaleExp, decreasing = T), 'Gene'])
plot.bubble <- 
    ggplot(data = df.bubble, 
           aes(x = Gene, y = Cell, size = Prop, color = ScaleExp)) + 
    geom_point(fill = 'cornsilk') + 
    scale_colour_gradientn(
        colours = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) + 
    # facet_grid( ~ Status) + 
    labs(x = 'Gene', y = 'Cell Type', color = 'Scaled expression', 
         size = 'Proportion') + 
    theme(panel.background = element_rect(color = 'gray',
                                          fill = 'transparent'),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(plot = plot.bubble, path = path.data, 
       filename = 'marker_surface_no_filter_chon.png',
       height = 18, width = 35, units = 'cm')

df.stroma <- df.bubble[df.bubble$Cell == 'Stromal stem cell', ]
df.bubble$Gene <- factor(df.bubble$Gene, 
                         levels = df.stroma[order(df.stroma$ScaleExp, decreasing = T), 'Gene'])
plot.bubble <- 
    ggplot(data = df.bubble, 
           aes(x = Gene, y = Cell, size = Prop, color = ScaleExp)) + 
    geom_point(fill = 'cornsilk') + 
    scale_colour_gradientn(
        colours = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) + 
    # facet_grid( ~ Status) + 
    labs(x = 'Gene', y = 'Cell Type', color = 'Scaled expression', 
         size = 'Proportion') + 
    theme(panel.background = element_rect(color = 'gray',
                                          fill = 'transparent'),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(plot = plot.bubble, path = path.data, 
       filename = 'marker_surface_no_filter_stroma.png',
       height = 18, width = 35, units = 'cm')


# violin plot
library(ggplot2)
marker.genes <- c('APOE', 'VCAN', 'OGN', 'COCH', 
                  'COL1A1', 'LUM', 'HES1', 'ID3', 'SOX9', 
                  'CYTL1', 'COL2A1', 'COL9A2', 'ACAN',
                  'ACTA2', 'PDGFRB', 'CDH5', 'IL1B')
# marker.genes <- c('CFD', 'APOE', 'VCAN', 'MMP10', 'COCH', 
#                   'ASPN', 'OGN')
df.gene <- data.frame(stringsAsFactors = F)
for (gene in marker.genes) {
    df.sub <- data.frame(expvalue = seurat.all.filter@assays$RNA@data[gene,],
                         gene = rep(gene, ncol(seurat.all.filter@assays$RNA@data)),
                         celltype = seurat.all.filter$celltype)
    df.gene <- rbind(df.gene, df.sub)
}
df.plot <- df.gene
df.plot$gene <- factor(df.gene$gene, levels = marker.genes)
df.plot$celltype <- factor(df.gene$celltype, 
                           levels = c('Stromal cell1', 'Stromal cell2', 'Stromal stem cell', 
                                      'Chondral stem cell', 'Transitional chondrocyte',
                                      'Chondrocyte1', 'Chondrocyte2','Perivascular cell',
                                      'Endothelial cell', 'Immune cell'))
plot.vln <- 
    ggplot(data = df.plot, aes(x = gene, y = expvalue, color = gene, fill = gene)) + 
    geom_violin(trim = T, scale = 'width') + 
    facet_grid( ~ celltype) + 
    theme_classic() + coord_flip() +
    stat_summary(fun= mean, geom = "point",
                 shape = 23, size = 2, color = "black") + 
    labs(x = 'Gene', y = 'Expression Level') + theme(legend.position = 'none')

ggsave(plot = plot.vln, path = path.data, 
       filename = 'Vln.png',
       height = 15, width = 35, units = 'cm')


# prop
seurat.4 <- subset(seurat.all, subset = RNA_snn_res.0.6 == 4)
table(seurat.4$batch)/table(seurat.all$batch)
df.plot <- data.frame(sample = names(table(seurat.4$batch)), 
                      prop = as.numeric(table(seurat.4$batch)/table(seurat.all$batch)))
ggplot(df.plot, aes(x = sample, y = prop)) + geom_bar(stat = 'identity')

seurat.5 <- subset(seurat.all, subset = RNA_snn_res.0.6 == 5)
table(seurat.5$batch)
seurat.9 <- subset(seurat.all, subset = RNA_snn_res.0.6 == 9)
table(seurat.9$batch)

seurat.11 <- subset(seurat.all, subset = RNA_snn_res.2 == 11)
DimPlot(seurat.11, group.by = 'RNA_snn_res.2')
seurat.0 <- subset(seurat.all, subset = RNA_snn_res.2 == 0)
DimPlot(seurat.0, group.by = 'RNA_snn_res.2')

seurat.M3 <- subset(seurat.all, subset = batch == 'M3')
FeaturePlot(seurat.M3, features = c('HES1', 'ID3', 'COL1A1', 'CYTL1'))

seurat.M1 <- subset(seurat.all, subset = batch == 'M1')
FeaturePlot(seurat.M1, features = c('HES1', 'ID3', 'COL1A1', 'CYTL1'))
FeaturePlot(seurat.M1, features = c('TNF'))

seurat.C4 <- subset(seurat.all, subset = batch == 'C4')
FeaturePlot(seurat.C4, features = c('HES1', 'ID3', 'COL1A1', 'CYTL1'))

seurat.C6 <- subset(seurat.all, subset = batch == 'C6')
FeaturePlot(seurat.C6, features = c('HES1', 'ID3', 'COL1A1', 'CYTL1'))
FeaturePlot(seurat.C6, features = c('TNF'))

seurat.C1 <- subset(seurat.all, subset = batch == 'C1')
FeaturePlot(seurat.C1, features = c('TNF', 'TNFSF10'))

(table(seurat.4$batch) + table(seurat.5$batch))/table(seurat.all$batch)

FeaturePlot(seurat.first, features = c('FRZB', 'CTGF', 
                                       'SERPINA1', "SCRG1", 'COL9A3', 'FGFBP2'), ncol = 3)

