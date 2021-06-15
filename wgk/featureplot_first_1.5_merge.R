setwd('/home/zy/my_git/bioinformatics/wgk')
.libPaths('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)

path.data <- '/home/disk/drizzle/wgk/data/marker_1.5_merge/'
file.merge_1.5 <- paste0(path.data, 'seurat_first.Rdata')
seurat.first <- readRDS(file.merge_1.5)

# stem marker
file.bone <- '/home/disk/drizzle/wgk/bone_genes.txt'
df.genes.bone <- read.delim(file.bone, row.names = 1)
bone.genes <- intersect(rownames(seurat.first@assays$RNA@scale.data), rownames(df.genes.bone))
clusters <- as.character(seurat.first$vec.cluster)
sel.genes <- c()
for (cluster in as.character(c(1, 3, 5, 7, 4, 9))) {
    sel.genes <- c(sel.genes, intersect(bone.genes, rownames(list.marker[[cluster]])))
}
sel.genes <- unique(sel.genes)

clusters <- as.character(seurat.first$vec.cluster)
mat_exp <- seurat.first@assays$RNA@scale.data[sel.genes, clusters == '1']
annotation_c <- data.frame(cluster = rep('1', ncol(mat_exp)))
for (cluster in as.character(c(3, 5, 7, 4, 9))) {
    mat_exp <- cbind(mat_exp, seurat.first@assays$RNA@scale.data[sel.genes, clusters == cluster])
    annotation_c <- rbind(annotation_c, 
                          data.frame(cluster = rep(cluster, ncol(seurat.first@assays$RNA@scale.data[sel.genes, clusters == cluster]))))
}
row.names(annotation_c) <- colnames(mat_exp)
seurat.plot <- seurat.first
seurat.plot@assays$RNA@scale.data <- mat_exp
seurat.plot@meta.data <- seurat.plot@meta.data[colnames(mat_exp),]
seurat.plot$heatmap <- factor(annotation_c$cluster, as.character(c(1, 3, 5, 7, 4, 9)))
DoHeatmap(seurat.plot, features = sel.genes, cells = colnames(mat_exp), group.by = 'heatmap')


# feature plot 
# 1/3/5
# SC
FeaturePlot(seurat.first, features = c('COL1A1', 'COL1A2', 'COL3A1', 'OGN', 'DCN', 'LUM'), ncol = 3)
# FB
FeaturePlot(seurat.first, features = c('COL1A1', 'COL1A2', 'COL3A1', 'FBLN1', 'FBLN2', 'S100A4'), ncol = 3)
# 135
FeaturePlot(seurat.first, features = c('COL1A1', 'COL1A2', 'COL3A1', 
                                       'FBLN1', 'FBLN2', "DCN", 
                                       'OGN', 'VCAN', 'LUM'), ncol = 3)
FeaturePlot(seurat.first, features = c('CXCL14', 'MEDAG', 'FGF7', 
                                       'SERPINF1', "NAMPT", 'IL6'), ncol = 3)
FeaturePlot(seurat.first, features = c('FGF7', 'CXCL14', 'GSN', 
                                       'SERPINF1', "FGL2", "NAMPT", 
                                       'BASP1', 'RND1', 'FGF7', 
                                       'CXCL14', 'ABI3BP', 'VEGFA', 
                                       'TSHZ2', 'PTGES'), ncol = 4)


# 5
FeaturePlot(seurat.first, features = c('ID3', 'SKIL', 'MYC', 
                                       'CLK1', "MAP3K2", "MAP3K8"), ncol = 3)
FeaturePlot(seurat.first, features = c('GPNMB', "CBX4", "ATF3", 
                                       'FOS', 'FOSB', "JUN", "JUNB", 'ASPN'), ncol = 4)
# quiesence
FeaturePlot(seurat.first, features = c('ID3', 'NR4A1', 'HES1'))



# 1
# FeaturePlot(seurat.first, features = row.names(list.marker$`1`)[1:20], ncol = 5)
FeaturePlot(seurat.first, features = c('CFD', 'APOE', 'CSF3', 
                                       'MEDAG', "VCAN", "IL33", 
                                       'TYMP', 'VEGFA', 'PTGES'), ncol = 3)


# 3
FeaturePlot(seurat.first, features = c('MMP10', 'COCH', 'IGFBP6', 
                                       'COMP', "OGN", "AEBP1"), ncol = 3)


# 4/9
FeaturePlot(seurat.first, features = c('SERPINA1', 'TNFRSF11B', 
                                       'FN1', "C2orf82", 'LAMB3', 'SERPINE2'), ncol = 3)
FeaturePlot(seurat.first, features = c('COL9A2', "COL9A3", 'CSF1', 
                                       'CHI3L2', 'CHI3L1', "LOX"), ncol = 3)
FeaturePlot(seurat.first, features = c('ACAN', "BMP2", 'SOD3', 
                                       'LEPREL1', 'COL11A2', "WIF1"), ncol = 3)

# 9
FeaturePlot(seurat.first, features = c('STC1'))

# 4
FeaturePlot(seurat.first, features = c('FRZB', 'CTGF', 
                                       'SERPINA1', "SCRG1", 'COL9A3', 'FGFBP2'), ncol = 3)
FeaturePlot(seurat.first, features = c('PLA2G2A', 'KCNMA1', 
                                       'FGFBP2', "COL9A3", 'S100A1', 'SHDA'), ncol = 3)



# 7
FeaturePlot(seurat.first, features = c('CYTL1', "DUSP2", 'FOS', 
                                       'JUN', 'SERPINA1', "C2orf82"), ncol = 3)

# 11
FeaturePlot(seurat.first, features = c('IFIT1', "OAS3", "OAS2", "OAS1", 'XAF1', 
                                       'IFI44L', 'MX1', "RSAD2"), ncol = 4)

# 2
FeaturePlot(seurat.first, features = c('ACTA2', "MYH11", 'CNN1', 
                                       'SYNPO2', 'MCAM', "PDE3A"), ncol = 3)


# 2
FeaturePlot(seurat.first, features = c('ACTA2', "MYH11", 'CNN1', 
                                       'SYNPO2', 'MCAM', "PDE3A"), ncol = 3)


# 6
FeaturePlot(seurat.first, features = c('CDH5', "EMCN", 'PLVAP', 
                                       'PECAM1', 'ESAM', "SELE"), ncol = 3)


# 10
FeaturePlot(seurat.first, features = c('IL1B', "CCL3", 'CCL4', 
                                       'FCER1G', 'IL24', "IL10"), ncol = 3)


# cell type
vec.cluster <- seurat.first$vec.cluster
celltypes[vec.cluster == '1'] <- 'Stromal cell1'
celltypes[vec.cluster == '3'] <- 'Stromal cell2'
celltypes[vec.cluster == '5'] <- 'CSPC'
celltypes[vec.cluster == '7'] <- 'Pre Chondrocyte'
celltypes[vec.cluster == '4'] <- 'Chondrocyte1'
celltypes[vec.cluster == '9'] <- 'Chondrocyte2'
celltypes[vec.cluster == '11'] <- 'Unknown'
celltypes[vec.cluster == '2'] <- 'Smooth muscle cell'
celltypes[vec.cluster == '6'] <- 'Endothelial cell'
celltypes[vec.cluster == '8'] <- 'Neural cell'
celltypes[vec.cluster == '10'] <- 'Immune cell'
seurat.first$celltype <- celltypes
DimPlot(seurat.first, group.by = "celltype", label = T)
file.merge_1.5 <- paste0(path.data, 'seurat_first.Rdata')
saveRDS(seurat.first, file.merge_1.5)

# violin plot
library(ggplot2)
marker.genes <- c('CFD', 'APOE', 'VCAN', 'MMP10', 'COCH', 
                  'ASPN', 'OGN', 
                  'COL1A1', 'COL3A1', 'FBLN2', 
                  'FOS', 'CLK1', 'ID3', 'SKIL',
                  'CYTL1', 'COL9A2', 'COL9A3', 'ACAN', 'SERPINA1', 'C2orf82',
                  'ACTA2', 'MCAM', 'CDH5', 'ESAM', 'IL1B', 'CCL3', 
                  'ANK3', 'NTM', 'IFIT1', 'OAS3')
# marker.genes <- c('CFD', 'APOE', 'VCAN', 'MMP10', 'COCH', 
#                   'ASPN', 'OGN')
df.gene <- data.frame(stringsAsFactors = F)
for (gene in marker.genes) {
    df.sub <- data.frame(expvalue = seurat.first@assays$RNA@data[gene,],
                         gene = rep(gene, ncol(seurat.first@assays$RNA@data)),
                         celltype = seurat.first$celltype)
    df.gene <- rbind(df.gene, df.sub)
}
df.plot <- df.gene
df.plot$gene <- factor(df.gene$gene, levels = marker.genes)
df.plot$celltype <- factor(df.gene$celltype, 
                           levels = c('Stromal cell1', 'Stromal cell2', 'CSPC', 'Pre Chondrocyte',
                                      'Chondrocyte1', 'Chondrocyte2','Smooth muscle cell',
                                      'Endothelial cell', 'Immune cell', 'Neural cell', 'Unknown'))
ggplot(data = df.plot, aes(x = gene, y = expvalue, color = gene, fill = gene)) + 
    geom_violin(trim = T, scale = 'width') + 
    # geom_boxplot() + 
    facet_grid( ~ celltype) + 
    theme_classic() + coord_flip() +
    stat_summary(fun= mean, geom = "point",
                 shape = 23, size = 2, color = "black") + 
    # ylim(1,8) + 
    labs(x = 'Gene', y = 'Expression Level')
# VlnPlot(seurat.first, features = c('EDNRB', 'APOLD1', 'PDGFRB', 'SYNPO2', 'OLFM2', 'PDE3A'), 
#         group.by = 'celltype', pt.size = 0, ncol = 3)



# cell research
FeaturePlot(seurat.first, features = c('MSX1', 'HOXA9', 'PRRX1', 'HOXC6'))
FeaturePlot(seurat.first, features = c('TWIST2', 'PDGFRA', 'RUNX2', 'OSR2'))
FeaturePlot(seurat.first, features = c('NOV', 'SFRP2', 'SOX9', 'ACAN'))
FeaturePlot(seurat.first, features = c('SIX1', 'MYOG', 'CDH5', 'CD68'))
FeaturePlot(seurat.first, features = c('GYPA', 'FGF8', 'EPCAM', 'SOX10'))

FeaturePlot(seurat.first, features = c('PITX1', 'HOXA10', 'CRABP1', 'CD24'))
FeaturePlot(seurat.first, features = c('GAS1', 'GAS2', 'SOX4', 'SFRP2'))
FeaturePlot(seurat.first, features = c('CNMD', 'EPYC', 'COL9A2', 'COL11A2'))
FeaturePlot(seurat.first, features = c('DLX5', 'CDH11', 'OGN', 'COL1A1', 'COL1A2'), ncol = 3)

FeaturePlot(seurat.first, features = c('NR2F2', 'SHOX2', 'EBF1'))
FeaturePlot(seurat.first, features = c('MSX1', 'NR2F1', 'HOXA5'))
FeaturePlot(seurat.first, features = c('MECOM', 'FOXP1', 'FOXP2'))
FeaturePlot(seurat.first, features = c('DLX5', 'TBX2', 'ATF1'))
FeaturePlot(seurat.first, features = c('CEBPG', 'EGR1', 'ATF3'))
FeaturePlot(seurat.first, features = c('RARG', 'FOXP4', 'SOX6'))
FeaturePlot(seurat.first, features = c('SOX8', 'SOX9', 'MEF2C'))

FeaturePlot(seurat.first, features = c('NGFR', 'NES', 'BMP4', 'FOXC2'))
FeaturePlot(seurat.first, features = c('GJA1', 'SLC6A13', 'PTGDS'))
FeaturePlot(seurat.first, features = c('RUNX2', 'DLX5', 'CLEC11A'))
FeaturePlot(seurat.first, features = c('OSR2', 'POSTN', 'SOX9', 'COL9A2'))
FeaturePlot(seurat.first, features = c('PDGFRB', 'MCAM', 'ACTA2', 'MYF5'))


####################### microtia
FeaturePlot(seurat.first, features = c('PDGFRB', 'MCAM', 'ACTA2', 'MYF5'))





