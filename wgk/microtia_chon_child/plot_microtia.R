.libPaths('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)
library(ggplot2)
library(reshape2)

path.data <- '/home/disk/drizzle/wgk/data/AllSample_2_merge/'
path.lineage <- paste0(path.data, 'chon_lineage/')
file.chon <- paste0(path.lineage, 'seurat_celltype.Rdata')
seurat.chon <- readRDS(file.chon)
# M1 M2 M3
seurat.child <- subset(seurat.chon, subset = batch %in% c('C4', 'C6', 'M1', 'M2', 'M3'))

# prop
table.prop <- table(seurat.child$batch, seurat.child$celltype)/as.vector(table(seurat.child$batch))
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

fc.cutoff <- 0.4
path.M123 <- '/home/disk/drizzle/wgk/microtia_chon_child_M1M2M3/'
path.cutoff <- paste0(path.M123, 'cutoff_', fc.cutoff, '/')
# path.M12 <- '/home/disk/drizzle/wgk/microtia_child_M1M2/'
# path.cutoff <- paste0(path.M12, 'cutoff_', fc.cutoff, '/')
file.marker.go <- paste0(path.cutoff, 'marker_go.Rdata')
list.marker.go <- readRDS(file.marker.go)


file.go.BP <- paste0(path.cutoff, 'GO_BP.Rdata')
list.go.BP <- readRDS(file.go.BP)
file.go.MF <- paste0(path.cutoff, 'GO_MF.Rdata')
list.go.MF <- readRDS(file.go.MF)
# file.go.CC <- paste0(path.cutoff, 'GO_CC.Rdata')
# list.go.CC <- readRDS(file.go.CC)
# file.kegg <- paste0(path.cutoff, 'kegg.Rdata')
# list.kegg <- readRDS(file.kegg)

file.go.BP <- paste0(path.cutoff, 'GO_BP_all.Rdata')
list.go.BP <- readRDS(file.go.BP)
file.go.MF <- paste0(path.cutoff, 'GO_MF_all.Rdata')
list.go.MF <- readRDS(file.go.MF)

# select GO
# df.GO <- data.frame(stringsAsFactors = F)
# Chondral stem cell
# GO.BP.CSC.M <- c('response to oxidative stress', 'response to unfolded protein', 
#                  'response to tumor necrosis factor', 'RNA splicing',
#                  'RNA localization', 'positive regulation of defense response')
# GO.BP.CSC.N <- c('ribosome biogenesis', 'response to copper ion', 'oxidative phosphorylation',
#                  'extracellular matrix organization', 'skeletal system development',
#                  'cell aggregation', 'cellular zinc ion homeostasis')
sel.GO.BP.N <- c('protein localization to endoplasmic reticulum', 
                 'translational initiation', 'ribosome biogenesis', 
                 'oxidative phosphorylation', 'electron transport chain', 
                 'mitochondrion organization', 'cellular oxidant detoxification', 
                 'cartilage development', 'chondrocyte differentiation',
                 # 'cartilage condensation',
                 'cellular response to zinc ion', 'cellular response to copper ion',
                 'cell aggregation', 'extracellular matrix organization')
sel.GO.BP.M <- c('response to oxidative stress', 
                 'cell death in response to oxidative stress',
                 'nitric oxide biosynthetic process',
                 'I-kappaB kinase/NF-kappaB signaling',
                 'p38MAPK cascade', 'ERK1 and ERK2 cascade',
                 'intrinsic apoptotic signaling pathway', 
                 'extrinsic apoptotic signaling pathway',
                 'response to unfolded protein', 
                 'regulation of RNA stability',
                 'positive regulation of defense response', 
                 'cellular response to tumor necrosis factor',
                 'cellular response to interferon-gamma',
                 'cellular response to interleukin-1',
                 'regulation of inflammatory response',
                 'cell cycle arrest',
                 'negative regulation of stem cell differentiation',
                 'negative regulation of cell growth',
                 # 'RNA splicing', 'RNA localization', 
                 'angiogenesis', 
                 # 'vascular endothelial growth factor production', 
                 # 'positive regulation of vasculature development',
                 'positive regulation of cell migration',
                 'negative regulation of cell adhesion', 
                 'extracellular matrix organization')
# sel.GO.MF.N <- c('oxidoreductase activity, acting on a heme group of donors',
#                  # 'structural constituent of ribosome', 
#                  'extracellular matrix structural constituent', 
#                'extracellular matrix binding', 'S100 protein binding')
sel.GO.MF.N <- c('extracellular matrix structural constituent', 
                 'S100 protein binding')
sel.GO.MF.M <- c('extracellular matrix structural constituent')

terms <- c("CSC_Microtia_increase",
           "CSC_Microtia_decrease",
           "TC_Microtia_increase",
           "TC_Microtia_decrease",
           "C1_Microtia_increase",
           "C1_Microtia_decrease",
           "C2_Microtia_increase",
           "C2_Microtia_decrease")

df.plot <- data.frame(stringsAsFactors = F)
for (term in terms) {
    cell <- strsplit(term, split = '_')[[1]][1]
    status <- strsplit(term, split = '_')[[1]][3]
    sub.BP <- list.go.BP[[term]]
    rownames(sub.BP) <- sub.BP$Description
    sub.BP <- sub.BP[sub.BP$pvalue < 0.018,]
    if (status == 'increase') {
        sel.BP <- sub.BP[sel.GO.BP.M, c('Description', 'pvalue', 'geneID')]
        sel.BP$Description <- sel.GO.BP.M
    } else {
        sel.BP <- sub.BP[sel.GO.BP.N, c('Description', 'pvalue', 'geneID')]
        sel.BP$Description <- sel.GO.BP.N
    }
    sub.MF <- list.go.MF[[term]]
    rownames(sub.MF) <- sub.MF$Description
    sub.MF <- sub.MF[sub.MF$pvalue < 0.018,]
    if (status == 'increase') {
        sel.MF <- sub.MF[sel.GO.MF.M, c('Description', 'pvalue', 'geneID')]
        sel.MF$Description <- sel.GO.MF.M
    } else {
        sel.MF <- sub.MF[sel.GO.MF.N, c('Description', 'pvalue', 'geneID')]
        sel.MF$Description <- sel.GO.MF.N
    }
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
# df.plot$Description <- factor(df.plot$Description, 
#                               levels = unique(c(rev(sel.GO.MF.N), 
#                                                 rev(sel.GO.BP.N), 
#                                                 rev(sel.GO.BP.M))))
df.plot$Description <- factor(df.plot$Description, 
                              levels = rev(unique(c(rev(sel.GO.MF.N), 
                                                rev(sel.GO.BP.N), 
                                                rev(sel.GO.BP.M)))))
df.plot$CellType <- factor(df.plot$CellType, 
                           levels = c('CSC', 'TC', 'C1', 'C2'))

df.plot$col_name <- paste(df.plot$CellType, df.plot$Status, sep = '_')
mat.plot <- reshape2::dcast(df.plot, Description ~ col_name, value.var = 'log10Pval')
row.names(mat.plot) <- mat.plot$Description
mat.plot$Description <- NULL
mat.plot[is.na(mat.plot)] <- 0
# write
write.table(df.plot, file = paste0(path.cutoff, 'GO.txt'), sep = '\t')

# col annotation
annotation_col = data.frame(
    CellType = factor(c(rep('CSC', 2), rep('C1', 2), rep('C2', 2), rep('TC', 2)), 
                      levels = c('CSC', 'TC', 'C1', 'C2')), 
    Status = factor(rep(c('Microtia', 'Normal'), 4), levels = c('Normal', 'Microtia')),
    row.names = colnames(mat.plot)
)

cols <- c("CSC_Microtia", "TC_Microtia",
          "C1_Microtia", "C2_Microtia",
          "CSC_Normal", "TC_Normal",
          "C1_Normal", "C2_Normal")
mat.plot <- mat.plot[rev(rownames(mat.plot)), cols]
annotation_col <- annotation_col[cols,]

# ann_colors = list(
#     CellType = rep(c("#33A02C","#B2DF8A","#1F78B4","#A6CEE3"), 2),
#     Status = c(Microtia = rep("#F8766D",4), Normal = rep("#00BFC4",4)))
ann_colors = list(
    CellType = c(CSC="#EE9572",TC="#B2DF8A",C1="#A6CEE3",C2="#9999FF"),
    Status = c(Microtia = "#637FBF", Normal = "#6C6C6C")
)
cols.sort <- c("CSC_Normal", "TC_Normal",
          "C1_Normal", "C2_Normal",
          "CSC_Microtia", "TC_Microtia", 
          "C1_Microtia", "C2_Microtia")
mat.plot <- mat.plot[, cols.sort]

bk <- c(seq(-10,-0.1,by=0.01),seq(0,10,by=0.01))
labels_row <- gsub('oxidoreductase activity, acting on a heme group of donors', 
                   'oxidoreductase activity', rownames(mat.plot))
pheatmap::pheatmap(mat.plot,
                   color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
                   cluster_rows = F, cluster_cols = F, scale = "none",
                   display_numbers = F,
                   annotation_col = annotation_col, annotation_colors = ann_colors,
                   show_rownames = T, show_colnames = F, legend = T, 
                   fontsize_row = 10, labels_row = labels_row,
                   gaps_col = c(4), 
                   filename = paste0(path.cutoff, 'heatmap_GO_enrich.pdf'), 
                   width = 6.5, height = 6.5
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
# genes <- na.omit(unique(unlist(strsplit(df.plot$geneID, split = '/'))))
# sel.genes <- c('SOD2', 'NFE2L2', 'KDM6B',  # 氧化应激
#                'NFKB1', 'NFKB2', 'RELB', 'NFKBIA', 'NFKBIZ', # NF kappa
#                'MAP2K3', 'MAP2K2', 'DUSP1', 'TNIP1', # MAPK ERK
#                'CYLD', 'PMAIP1', 'FOXO3', 'FNIP2', 'RELA',  # 凋亡
#                'PTGES3', 'TCP1', 'HSPA9', 'HNRNPU', # 维持蛋白质/RNA稳定
#                'TNFAIP3', 'CD44', 'HLA-C', 'IFNGR2', # 免疫
#                'KMT2E', 'ZFP36L2', 'PPP2CA', 'CDKN1A', 'TSPYL2', # 抑制干细胞
#                'PDPN', 'DLC1', 'SDCBP', 'CSF1', 'CLEC3A', # 细胞迁移与粘连
#                'TNFAIP2', 'ESM1', 'ZC3H12A', # 血管生成
#                'MMP1', 'MMP3', # 基质解组装
#                'LAMB3', 'LAMC2', 'FN1', # extracellular matrix
#                'COL2A1', 'COL9A3', 'COL11A1', 'SPARC', 'BGN', 'SPINT2', 'A2M', 'ACAN',
#                'MGP', 'ELN', 
#                'CYTL1', 'FRZB', 'CTGF', 'CHAD', 'VIT', 'SFRP2', # 软骨发育
#                'S100A1', 'S100B', # S100 protein
#                'MT1X', 'MT1E', 'MT1G',  # 离子
#                'NDUFA1', 'SDHA', 'COX7A1', # 线粒体/呼吸链
#                'PRDX1', 'SOD1', 'SOD3', # 抗氧化
#                'RPS28', 'EIF3M') # 核糖体生成与翻译起始
sel.genes <- c('SOD2', 'NFE2L2', 'KDM6B',  # 氧化应激
               'PRDX1', 'SOD1', 'SOD3', # 抗氧化
               'NDUFA1', 'SDHA', 'COX7A1', # 线粒体/呼吸链
               'NFKB1', 'NFKB2', 'RELB', 'NFKBIA', 'NFKBIZ', # NF kappa
               'MAP2K3', 'MAP2K2', 'DUSP1', 'TNIP1', # MAPK ERK
               'CYLD', 'PMAIP1', 'FOXO3', 'FNIP2', 'RELA',  # 凋亡
               'PTGES3', 'TCP1', 'HSPA9', 'HNRNPU', # 维持蛋白质/RNA稳定
               'TNFAIP3', 'CD44', 'HLA-C', 'IFNGR2', # 免疫
               'KMT2E', 'ZFP36L2', 'PPP2CA', 'CDKN1A', 'TSPYL2', # 抑制干细胞
               'PDPN', 'DLC1', 'SDCBP', 'CSF1', 'CLEC3A', # 细胞迁移与粘连
               'TNFAIP2', 'ESM1', 'ZC3H12A', # 血管生成
               'MMP1', 'MMP3', # 基质解组装
               'LAMB3', 'LAMC2', 'FN1', # extracellular matrix
               'COL2A1', 'COL9A3', 'COL11A1', 'SPARC', 'BGN', 'SPINT2', 'A2M', 'ACAN',
               'MGP', 'ELN', 
               'CYTL1', 'FRZB', 'CTGF', 'CHAD', 'VIT', 'SFRP2', # 软骨发育
               'S100A1', 'S100B', # S100 protein
               'MT1X', 'MT1E', 'MT1G',  # 离子
               'RPS28', 'EIF3M') # 核糖体生成与翻译起始
sel.genes <- rev(sel.genes)
cells <- c('CSC', 'TC', 'C1', 'C2')
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
df.plot.gene$log_pval[df.plot.gene$log_pval > 200] <- 200
df.plot.gene$log_log_pval <- log(df.plot.gene$log_pval + 1)
df.plot.gene$Gene <- factor(df.plot.gene$Gene, levels = sel.genes)
df.plot.gene$CellType <- factor(df.plot.gene$CellType, levels = cells)

plot.bubble <- 
    ggplot(data = df.plot.gene, 
           aes(x = CellType, y = Gene, size = log_pval, color = avg_logFC)) + 
    geom_point(fill = 'cornsilk') + 
    scale_size_continuous(range = c(2,5)) +
    scale_colour_gradient2(low = 'blue', mid = "white", high = 'red',
                           breaks = c(-1, 0, 1)) + 
    scale_y_discrete(position = 'right') + 
    scale_x_discrete(breaks = cells,
                     labels = c('CSC', 'TC', 'C1', 'C2')) +
    # facet_grid( ~ Status) + 
    labs(x = '', y = '', color = 'logFC', 
         size = expression(paste("-log"[10], "(", italic("P"), "-value)"))) + 
    theme(panel.background = element_rect(color = 'gray',
                                          fill = 'transparent'),
          axis.text.y = element_text(
              size = 13, color = "black", face = 'bold.italic'), 
          axis.text.x = element_text(size = 15, color = "black", 
                                     angle = 45, vjust = 1, hjust = 1))

ggsave(plot = plot.bubble, path = path.cutoff,
       filename = 'marker_genes.png',
       height = 40, width = 11, units = 'cm')
ggsave(plot = plot.bubble, path = path.cutoff,
       filename = 'marker_genes.pdf',
       height = 40, width = 11, units = 'cm')


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
mat.chon <- seurat.child@assays$RNA@data
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
        sub.seurat <- subset(seurat.child, subset = celltype == cell & type == status)
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

seurat.child$heatmap <- factor(seurat.child$celltype, levels = cells)
plot.heatmap <- 
    DoHeatmap(seurat.child, features = sel.genes, cells = cell_ids.sort, 
              group.by = 'heatmap', size = 6, draw.lines = 1) + guides(color = F)
ggsave(plot = plot.heatmap, path = path.cutoff, 
       filename = 'heatmap.png', units = 'cm',
       height = 50, width = 50)


# marker gene
VlnPlot(seurat.child, features = c('IL8', 'CXCL2', 'CXCL1', 'MMP3', 'CSF1', 'SOD2'), 
        split.by = 'batch', pt.size = 0, group.by = 'celltype', ncol = 2)
VlnPlot(seurat.all, features = c('IL8', 'CXCL2', 'CXCL1', 'MMP3', 'CSF1', 'SOD2'), 
        split.by = 'batch', pt.size = 0, group.by = 'celltype', ncol = 2)

seurat.child.allsample <- subset(seurat.all, subset = celltype %in% 
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


# GSEA





