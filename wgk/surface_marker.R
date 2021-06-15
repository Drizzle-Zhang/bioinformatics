setwd('/home/zy/my_git/bioinformatics/wgk')
.libPaths('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)
library(ggplot2)
require("RColorBrewer")

path.data <- '/home/disk/drizzle/wgk/data/AllSample_3_merge/'
file.merge_3 <- paste0(path.data, 'seurat_celltype.Rdata')
seurat.all.filter <- readRDS(file.merge_3)

# surface marker
file.surface <- '/home/disk/drizzle/wgk/public_data/surface_marker_wlab.txt'
df.surface <- read.delim2(file.surface)
surface.genes <- unique(toupper(df.surface$ENTREZ.gene.symbol))

# chon and stroma lineage marker genes
seurat.CS <- subset(seurat.all.filter, subset = 
                        celltype %in% c('Stromal cell1', 'Stromal cell2', 'Stromal stem cell',
                                        'Chondral stem cell', 'Transitional chondrocyte', 
                                        'Chondrocyte1', 'Chondrocyte2'))
file.merge_3.CS <- paste0(path.data, 'seurat_celltype_CS.Rdata')
saveRDS(seurat.CS, file.merge_3.CS)
DimPlot(seurat.CS, group.by = "celltype", label = T)

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
# genes.chon <- c("CD83", "TSPAN6", 'ENPP1', 'ITM2B', 'CD99', 'TSPAN4', 
#                 'SCARA3', 'BOC', 'PTPRZ1', 'MXRA8', 'FGFR3', 'AQP1',
#                 'CNTFR', 'CADM1', 'SLC29A1', 'LRP1', 'CD46', 'A2M', 'ITGA10')
genes.chon <- c("CD83", "TSPAN6", 'ENPP1', 'ITM2B', 'CD99', 'TSPAN4', 
                'SCARA3', 'BOC', 'PTPRZ1', 'MXRA8', 'FGFR3', 'AQP1',
                'CNTFR', 'SLC29A1', 'LRP1', 'CD46', 'A2M', 'ITGA10')
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
       height = 16, width = 35, units = 'cm')

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
sel.cells <- c('Stromal cell1', 'Stromal cell2', 'Stromal stem cell', 
               'Chondral stem cell', 'Transitional chondrocyte',
               'Chondrocyte1', 'Chondrocyte2','Perivascular cell',
               'Endothelial cell', 'Immune cell')
# chondral stem cell
df.chon <- df.bubble[df.bubble$Cell == 'Chondral stem cell', ]
df.bubble$Gene <- factor(df.bubble$Gene, 
                         levels = df.chon[order(df.chon$ScaleExp, decreasing = T), 'Gene'])
df.bubble.chon <- df.bubble[df.bubble$Gene %in% genes.chon & df.bubble$Cell %in% sel.cells, ]
plot.bubble <- 
    ggplot(data = df.bubble.chon, 
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
       filename = 'marker_surface_chon.png',
       height = 18, width = 20, units = 'cm')

df.stroma <- df.bubble[df.bubble$Cell == 'Stromal stem cell', ]
df.bubble$Gene <- factor(df.bubble$Gene, 
                         levels = df.stroma[order(df.stroma$ScaleExp, decreasing = T), 'Gene'])
df.bubble.stroma <- df.bubble[df.bubble$Gene %in% genes.stroma & df.bubble$Cell %in% sel.cells, ]
plot.bubble <- 
    ggplot(data = df.bubble.stroma, 
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
       filename = 'marker_surface_stroma.png',
       height = 18, width = 30, units = 'cm')
