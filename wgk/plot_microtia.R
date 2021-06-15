.libPaths('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)
library(ggplot2)
library(reshape2)

path.data <- '/home/disk/drizzle/wgk/data/AllSample_2_merge/'
file.all <- paste0(path.data, 'seurat_celltype.Rdata')
seurat.all <- readRDS(file.all)
# M1 M2 M3
seurat.chon <- subset(seurat.all, subset = celltype %in% 
                        c('Chondrocyte1', 'Chondrocyte2', 'Chondral stem cell', 
                          'Transitional chondrocyte') & batch %in% c('C4', 'C6', 'M1', 'M2', 'M3'))

# prop
table.prop <- table(seurat.chon$batch, seurat.chon$celltype)/as.vector(table(seurat.chon$batch))
df.prop <- melt(table.prop)
colnames(df.prop) <- c('Sample', 'CellType', 'Proportion')
df.prop$CellType <- factor(df.prop$CellType, 
                           levels = c('Chondral stem cell', 'Transitional chondrocyte', 
                                      'Chondrocyte1', 'Chondrocyte2'))
status <- rep('_', nrow(df.prop))
status[df.prop$Sample %in% c('C4', 'C6')] <- 'Normal'
status[df.prop$Sample %in% c('M1', 'M2', 'M3')] <- 'Microtia'
df.prop$Status <- status
plot.prop <- 
    ggplot(df.prop) + 
    geom_boxplot(aes(x = CellType, y = Proportion, color = Status)) + 
    geom_jitter(aes(x = CellType, y = Proportion, color = Status), 
                position = position_jitterdodge()) +
    theme_bw() + 
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(plot.prop, path = '/home/disk/drizzle/wgk/data/AllSample_2_merge', 
       filename = 'Prop.png', units = 'cm',
       width = 15, height = 10)

fc.cutoff <- 0.5
path.M123 <- '/home/disk/drizzle/wgk/microtia_child_M1M2M3/'
path.cutoff <- paste0(path.M123, 'cutoff_', fc.cutoff, '/')
# path.M12 <- '/home/disk/drizzle/wgk/microtia_child_M1M2/'
# path.cutoff <- paste0(path.M12, 'cutoff_', fc.cutoff, '/')
file.marker.go <- paste0(path.cutoff, 'marker_go.Rdata')
list.marker.go <- readRDS(file.marker.go)


file.go.BP <- paste0(path.cutoff, 'GO_BP.Rdata')
list.go.BP <- readRDS(file.go.BP)
file.go.MF <- paste0(path.cutoff, 'GO_MF.Rdata')
list.go.MF <- readRDS(file.go.MF)
file.go.CC <- paste0(path.cutoff, 'GO_CC.Rdata')
list.go.CC <- readRDS(file.go.CC)
file.kegg <- paste0(path.cutoff, 'kegg.Rdata')
list.kegg <- readRDS(file.kegg)

# select GO
df.GO <- data.frame(stringsAsFactors = F)
# Chondral stem cell
GO.BP.CSC.M <- c('response to oxidative stress', 'response to unfolded protein', 
                 'response to tumor necrosis factor', 'RNA splicing',
                 'RNA localization', 'positive regulation of defense response')
GO.BP.CSC.N <- c('ribosome biogenesis', 'response to copper ion', 'oxidative phosphorylation',
                 'extracellular matrix organization', 'skeletal system development',
                 'cell aggregation', 'cellular zinc ion homeostasis')
sel.GO.BP <- c('response to oxidative stress', 
               'response to unfolded protein', 
               'positive regulation of defense response', 
               # 'regulation of inflammatory response',
               'p38MAPK cascade', 'ERK1 and ERK2 cascade',
               # 'regulation of ERK1 and ERK2 cascade', 
               'intrinsic apoptotic signaling pathway', 
               'cell cycle arrest',
               'negative regulation of stem cell differentiation',
               'negative regulation of cell growth',
               'RNA splicing', 'RNA localization', 
               'vascular endothelial growth factor production', 'angiogenesis', 
               'positive regulation of vasculature development',
               'positive regulation of cell migration',
               'negative regulation of cell adhesion',
               'translational initiation', 'ribosome biogenesis', 
               'cartilage condensation', 
               'extracellular matrix organization', 
               # 'skeletal system development',
               'connective tissue development',
               'zinc ion homeostasis')
sel.GO.MF <- c('extracellular matrix structural constituent', 
               'extracellular matrix binding', 'S100 protein binding')

terms <- c("Chondral stem cell_Microtia_increase",
           "Chondral stem cell_Microtia_decrease",
           "Transitional chondrocyte_Microtia_increase",
           "Transitional chondrocyte_Microtia_decrease",
           "Chondrocyte1_Microtia_increase",
           "Chondrocyte1_Microtia_decrease",
           "Chondrocyte2_Microtia_increase",
           "Chondrocyte2_Microtia_decrease")

df.plot <- data.frame(stringsAsFactors = F)
for (term in terms) {
    cell <- strsplit(term, split = '_')[[1]][1]
    status <- strsplit(term, split = '_')[[1]][3]
    sub.BP <- list.go.BP[[term]]
    rownames(sub.BP) <- sub.BP$Description
    sub.BP <- sub.BP[sub.BP$p.adjust < 0.1,]
    sel.BP <- sub.BP[sel.GO.BP, c('Description', 'pvalue', 'geneID')]
    sel.BP$Description <- sel.GO.BP
    sub.MF <- list.go.MF[[term]]
    rownames(sub.MF) <- sub.MF$Description
    sub.MF <- sub.MF[sub.MF$p.adjust < 0.1,]
    sel.MF <- sub.MF[sel.GO.MF, c('Description', 'pvalue', 'geneID')]
    sel.MF$Description <- sel.GO.MF
    sub.plot <- rbind(sel.BP, sel.MF)
    sub.plot$pvalue[is.na(sub.plot$pvalue)] <- 1
    sub.plot$CellType <- rep(cell, nrow(sub.plot))
    if (status == 'increase') {
        sub.plot$Status <- rep('Microtia', nrow(sub.plot))
        sub.plot$Coeff <- rep(1, nrow(sub.plot))
    } else {
        sub.plot$Status <- rep('Normal', nrow(sub.plot))
        sub.plot$Coeff <- rep(-1, nrow(sub.plot))
    }
    df.plot <- rbind(df.plot, sub.plot)
}
df.plot$log10Pval <- -log10(df.plot$pvalue)
df.plot$log10Pval[abs(df.plot$log10Pval) > 10] = 10
df.plot$log10Pval <- df.plot$log10Pval * df.plot$Coeff
df.plot$abs_log10Pval <- abs(df.plot$log10Pval)
df.plot$Description <- factor(df.plot$Description, levels = c(sel.GO.BP, sel.GO.MF))
df.plot$CellType <- factor(df.plot$CellType, 
                           levels = c('Chondral stem cell', 'Transitional chondrocyte',
                                      'Chondrocyte1', 'Chondrocyte2'))

df.plot$col_name <- paste(df.plot$CellType, df.plot$Status, sep = '_')
mat.plot <- reshape2::dcast(df.plot, Description ~ col_name, value.var = 'log10Pval')
row.names(mat.plot) <- mat.plot$Description
mat.plot$Description <- NULL

# col annotation
annotation_col = data.frame(
    CellType = factor(c(rep('Chondral stem cell', 2),
                        rep('Chondrocyte1', 2),
                        rep('Chondrocyte2', 2),
                        rep('Transitional chondrocyte', 2)), 
                      levels = c('Chondral stem cell', 'Transitional chondrocyte',
                                 'Chondrocyte1', 'Chondrocyte2')), 
    Status = factor(rep(c('Microtia', 'Normal'), 4), levels = c('Normal', 'Microtia')),
    row.names = colnames(mat.plot)
)

cols <- c("Chondral stem cell_Microtia", "Transitional chondrocyte_Microtia", 
          "Chondrocyte1_Microtia", "Chondrocyte2_Microtia",
          "Chondral stem cell_Normal", "Transitional chondrocyte_Normal",
          "Chondrocyte1_Normal", "Chondrocyte2_Normal")
mat.plot <- mat.plot[rev(rownames(mat.plot)), cols]
annotation_col <- annotation_col[cols,]

pheatmap::pheatmap(mat.plot,
                   color = colorRampPalette(c('blue', 'white', 'red'))(100),
                   cluster_rows = F, cluster_cols = F, scale = "none",
                   display_numbers = F,
                   annotation_col = annotation_col ,
                   show_rownames = T, show_colnames = F, legend = T, 
                   # fontsize_row = 18, fontsize_col = 15,
                   gaps_col = c(4), 
                   filename = paste0(path.cutoff, 'heatmap_GO_enrich.png'), 
                   width = 8.5, height = 7
)


# plot.bubble <- 
#     ggplot(data = df.plot, 
#            aes(x = CellType, y = Description, size = abs_log10Pval, color = log10Pval)) + 
#     geom_point(fill = 'cornsilk') + 
#     scale_colour_gradient2(low = 'blue', mid = "white", high = 'red') + 
#     facet_grid( ~ Status) + 
#     labs(x = 'Cell Type', y = 'GO term', color = '-log10(P value)', 
#          size = 'abs(-log10(P value))') + 
#     theme(panel.background = element_rect(color = 'gray',
#                                           fill = 'transparent'),
#           axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
# 
# ggsave(plot = plot.bubble, path = path.cutoff, 
#        filename = 'GO_enrich.png',
#        height = 20, width = 20, units = 'cm')

# genes from GO
genes <- na.omit(unique(unlist(strsplit(df.plot$geneID, split = '/'))))
sel.genes <- c('SOD2', 'NFE2L2', 'PTGES3', 'TCP1', 'HSPA9', 'PMAIP1',
               'MAP2K3', 'MAP2K2', 'DUSP1', 'TNIP1', 
               'NFKB1', 'NFKB2', 'RELB', 'NFKBIA', 'NFKBIZ',
               'PPP2CA', 'CDKN1A', 'TSPYL2',
               'PDPN', 'CSF1', 'ICAM1', 'CLEC3A',
               'COL2A1', 'COL9A3', 'COL11A1', 'SPARC', 'BGN', 'SPINT2', 'A2M',
               'MGP', 'ELN', 'CD44', 'SERPINE2', 'LAMB3', 'LAMC2', 'FN1',
               'MMP1', 'MMP3', 'S100A1', 'S100B',
               'CYTL1', 'FRZB', 'CTGF', 'CHAD', 
               'ESM1', 'ZC3H12A')
cells <- c('Chondral stem cell', 'Transitional chondrocyte',
           'Chondrocyte1', 'Chondrocyte2')
file.gsea.marker <- paste0(path.M123, 'marker_GSEA.Rdata')
list.marker.gsea <- readRDS(file.gsea.marker)

df.plot.gene <- data.frame(stringsAsFactors = F)
for (cell in cells) {
    sub.marker <- list.marker.gsea[[cell]][sel.genes, c('p_val', 'avg_logFC')]
    sub.marker$CellType <- rep(cell, length(sel.genes))
    sub.marker$Gene <- sel.genes
    df.plot.gene <- rbind(df.plot.gene, sub.marker)
}
df.plot.gene$avg_logFC[is.na(df.plot.gene$avg_logFC)] <- 0
df.plot.gene$avg_logFC[df.plot.gene$avg_logFC < -1.5] <- -1.5
df.plot.gene$avg_logFC[df.plot.gene$avg_logFC > 1.5] <- 1.5
df.plot.gene$p_val[is.na(df.plot.gene$p_val)] <- 1
df.plot.gene$log_pval <- -log10(df.plot.gene$p_val)
df.plot.gene$log_pval[df.plot.gene$log_pval > 300] <- 300
df.plot.gene$log_log_pval <- log(df.plot.gene$log_pval + 1)
df.plot.gene$Gene <- factor(df.plot.gene$Gene, levels = sel.genes)
df.plot.gene$CellType <- factor(df.plot.gene$CellType, levels = cells)

plot.bubble <- 
    ggplot(data = df.plot.gene, 
           aes(x = CellType, y = Gene, size = log_pval, color = avg_logFC)) + 
    geom_point(fill = 'cornsilk') + 
    scale_colour_gradient2(low = 'blue', mid = "white", high = 'red',
                           breaks = c(-1, 0, 1)) + 
    scale_x_discrete(breaks = cells,
                     labels = c('CSC', 'TC', 'Chondrocyte1', 'Chondrocyte2')) +
    # facet_grid( ~ Status) + 
    labs(x = '', y = '', color = 'logFC', 
         size = '-log10(P value)') + 
    theme(panel.background = element_rect(color = 'gray',
                                          fill = 'transparent'),
          axis.text.y = element_text(
              size = 13, color = "black", face = 'bold.italic'), 
          axis.text.x = element_text(size = 15, color = "black", 
                                     angle = 45, vjust = 1, hjust = 1))

# ggsave(plot = plot.bubble, path = path.cutoff, 
#        filename = 'marker_genes.png',
#        height = 30, width = 13, units = 'cm')
ggsave(plot = plot.bubble, path = path.cutoff, 
       filename = 'marker_genes.pdf',
       height = 30, width = 13, units = 'cm')


# gene
library(ComplexHeatmap)
require("RColorBrewer")
genes <- c('CXCL1', 'CXCL2', 'CXCL3', 'IL8', 'CSF1', 
           'TNFRSF11B', 'ICAM1', 'NFKBIA', 'NFKBIZ', 'NFKB1', 'NFKB2', 
           'SOD2', 'NOS2', 'HSPA6', 'HSPA1B', 
           'MMP1', 'MMP3', 'COL1A1',
           'KLF2', 'KLF4',
           'COL2A1', 'COL11A1', 'COL9A3', 'MGP', 'ACAN', 'LAMB3', 'ELN', 
           'C2orf82', 'CYTL1', 'CTGF', 'FRZB', 
           'MT1X', 'MT1E', 'MT1M', 'S100B', 'S100A1',
           'RPS26', 'RPS28', 'RPS4X')

terms <- c("Chondral stem cell_Microtia_increase",
           "Transitional chondrocyte_Microtia_increase",
           "Chondrocyte1_Microtia_increase",
           "Chondrocyte2_Microtia_increase",
           "Chondral stem cell_Microtia_decrease",
           "Transitional chondrocyte_Microtia_decrease",
           "Chondrocyte1_Microtia_decrease",
           "Chondrocyte2_Microtia_decrease")
mat.chon <- seurat.chon@assays$RNA@data
genes <- c()
cells <- c('Chondral stem cell', 'Transitional chondrocyte', 
           'Chondrocyte1', 'Chondrocyte2')
for (cell in cells) {
    sub.markers <- list.marker.go[[as.character(cell)]]
    sub.markers <- sub.markers[order(sub.markers$avg_logFC, decreasing = T),]
    genes <- c(genes, row.names(sub.markers[1:40,]))
}
for (cell in cells) {
    sub.markers <- list.marker.go[[as.character(cell)]]
    sub.markers <- sub.markers[order(sub.markers$avg_logFC),]
    genes <- c(genes, row.names(sub.markers[1:40,]))
}
sel.genes <- unique(genes)

cell_ids.sort <- c()
df_info <- c()
types <- c('Microtia', 'Control')
for (cell in cells) {
    for (status in types) {
        sub.seurat <- subset(seurat.chon, subset = celltype == cell & type == status)
        cell_ids.sort <- c(cell_ids.sort, colnames(sub.seurat@assays$RNA@data))
        df_info <- rbind(df_info, sub.seurat@meta.data)
    }
}

mat.chon.sort <- t(scale(t(mat.chon[sel.genes, cell_ids.sort])))
mat.chon.sort[mat.chon.sort > 3] = 3
mat.chon.sort[mat.chon.sort < -3] = -3

cell_info <- factor(df_info$celltype, levels = cells)
color <- brewer.pal(4,"Set2")
names(color) <- levels(cell_info)
top_anno <- HeatmapAnnotation(
    cluster = anno_block(gp = gpar(fill = color), # 设置填充色
                         labels = levels(cell_info), 
                         labels_gp = gpar(cex = 0.5, col = "white"))) # 设置字体
Heatmap(as.matrix(mat.chon.sort),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = F,
        column_split = cell_info,
        column_title = NULL,
        top_annotation = top_anno)

seurat.chon$heatmap <- factor(seurat.chon$celltype, levels = cells)
plot.heatmap <- 
    DoHeatmap(seurat.chon, features = sel.genes, cells = cell_ids.sort, 
              group.by = 'heatmap', size = 6, draw.lines = 1) + guides(color = F)
ggsave(plot = plot.heatmap, path = path.cutoff, 
       filename = 'heatmap.png', units = 'cm',
       height = 50, width = 50)


# marker gene
VlnPlot(seurat.chon, features = c('IL8', 'CXCL2', 'CXCL1', 'MMP3', 'CSF1', 'SOD2'), 
        split.by = 'batch', pt.size = 0, group.by = 'celltype', ncol = 2)
VlnPlot(seurat.all, features = c('IL8', 'CXCL2', 'CXCL1', 'MMP3', 'CSF1', 'SOD2'), 
        split.by = 'batch', pt.size = 0, group.by = 'celltype', ncol = 2)

seurat.chon.allsample <- subset(seurat.all, subset = celltype %in% 
                          c('Chondrocyte1', 'Chondrocyte2', 'Chondral stem cell', 
                            'Transitional chondrocyte'))

VlnPlot(seurat.all, features = c('CD44'), 
        split.by = 'batch', pt.size = 0, group.by = 'celltype', ncol = 1)
VlnPlot(seurat.all, features = c('ICAM1'), 
        split.by = 'batch', pt.size = 0, group.by = 'celltype', ncol = 1)

VlnPlot(seurat.all, features = c('FRZB', 'CYTL1', 'SCRG1'), 
        split.by = 'batch', pt.size = 0, group.by = 'celltype', ncol = 1)
VlnPlot(seurat.all, features = c('RPS10', 'SSR3', 'RPS4X'), 
        split.by = 'batch', pt.size = 0, group.by = 'celltype', ncol = 1)
