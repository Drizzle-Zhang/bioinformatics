setwd('/home/zy/my_git/bioinformatics/wgk')
.libPaths('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)
library(harmony)
library(ggplot2)
require("RColorBrewer")
library(gridExtra)

path.data <- '/home/disk/drizzle/wgk/data/AllSample_2_merge/'
file.merge_2 <- paste0(path.data, 'seurat_celltype.Rdata')
seurat.all.filter <- readRDS(file.merge_2)

path.lineage <- paste0(path.data, 'stroma_lineage/')
if (!file.exists(path.lineage)) {
    dir.create(path.lineage)
}

seurat.stroma <- subset(seurat.all.filter, subset = celltype %in% 
                          c('Stromal cell1', 'Stromal cell2', 'Stromal stem cell'))
seurat.stroma <- NormalizeData(seurat.stroma)
seurat.stroma <- FindVariableFeatures(seurat.stroma, nfeatures = 3000)
seurat.stroma <- ScaleData(seurat.stroma, split.by = "batch")
seurat.stroma <- RunPCA(seurat.stroma, verbose = F, npcs = 50)
seurat.stroma <- RunHarmony(seurat.stroma, "batch", reduction.save = "harmony2")
seurat.stroma <- RunUMAP(seurat.stroma, reduction = "harmony2", 
                       dims = 1:50, n_neighbors = 30)
DimPlot(seurat.stroma, group.by = "celltype.abbr", label = T, reduction = 'umap', pt.size = 1)

seurat.stroma <- FindNeighbors(seurat.stroma, reduction = "harmony2", dims = 1:50)
seurat.stroma <- FindClusters(seurat.stroma, resolution = 0.15)
DimPlot(seurat.stroma, group.by = "RNA_snn_res.0.15", label = T, reduction = 'umap', pt.size = 1)

# seurat.stroma <- subset(seurat.stroma, 
#                             subset = RNA_snn_res.0.2 %in% 
#                                 setdiff(unique(seurat.stroma$RNA_snn_res.0.3), c(3)))

# marker cluster
clusters <- unique(seurat.stroma$RNA_snn_res.0.3)
list.marker.cluster <- list()
for (cluster in clusters) {
    sub.markers <- FindMarkers(seurat.stroma, ident.1 = cluster, group.by = 'RNA_snn_res.0.3',
                               logfc.threshold = 0.3, min.diff.pct = 0.05, only.pos = T)
    list.marker.cluster[[cluster]] <- sub.markers
}

cluster_0.3 <- seurat.stroma$RNA_snn_res.0.15
celltypes <- rep('_', length(cluster_0.3))
celltypes[cluster_0.3 %in% c(2)] <- 'SSC'
celltypes[cluster_0.3 %in% c(1)] <- 'SC1'
celltypes[cluster_0.3 %in% c(0)] <- 'SC2'
seurat.stroma$celltype <- celltypes

file.stroma <- paste0(path.lineage, 'seurat_celltype.Rdata')
saveRDS(seurat.stroma, file.stroma)
seurat.stroma <- readRDS(file.stroma)


colors <- c("#BC80BD", "#80B1D3", "#F4A460")
names(colors) <- c('SSC', 'SC1', 'SC2')

plot.umap <-
    DimPlot(seurat.stroma, reduction = 'umap', 
            group.by = "celltype", label = T,
            pt.size = 0.8, label.size = 5.5,
            cols = colors) + 
    # ylim(-12, 15) + 
    labs(x = 'UMAP1', y = 'UMAP2') +
    theme(panel.background=element_rect(fill='transparent', color='black',size = 1),
          legend.key=element_rect(fill='transparent', color='transparent'), 
          legend.text = element_text(size=16), 
          axis.title = element_text(size = 16),
          axis.text = element_blank(), 
          axis.ticks = element_blank()) + 
    guides(colour = guide_legend(override.aes = list(size=10)))
ggsave(plot.umap, path = path.lineage, 
       filename = 'stroma_lineage_umap.pdf',
       height = 10, width = 15, units = 'cm')


plot.cluster <-
    DimPlot(seurat.stroma, reduction = 'umap', 
            group.by = "RNA_snn_res.0.15", label = T,
            pt.size = 0.8, label.size = 5.5) + 
    # ylim(-12, 15) + 
    labs(x = 'UMAP1', y = 'UMAP2') +
    theme(panel.background=element_rect(fill='transparent', color='black',size = 1),
          legend.key=element_rect(fill='transparent', color='transparent'), 
          legend.text = element_text(size=16), 
          axis.title = element_text(size = 16),
          axis.text = element_blank(), 
          axis.ticks = element_blank()) + 
    guides(colour = guide_legend(override.aes = list(size=10)))
ggsave(plot.cluster, path = path.lineage, 
       filename = 'stroma_lineage_umap_cluster.pdf',
       height = 10, width = 15, units = 'cm')

plot.lineage <- 
    DimPlot(seurat.stroma, reduction = 'harmony', dims = c(2, 1),
            group.by = "celltype", label = T,
            cols = colors,
            pt.size = 0.8, label.size = 5.5) + 
    labs(x = 'PC2 (15%)', y = 'PC1 (17%)') + 
    theme(axis.title = element_text(size = 16),
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          axis.line = element_line(arrow = arrow(length = unit(0.5, 'cm'))),
          legend.position = 'none')
ggsave(plot.lineage, path = path.lineage, 
       filename = 'stroma_lineage_harmony.pdf',
       height = 10, width = 15, units = 'cm')

FeaturePlot(seurat.stroma, reduction = 'umap', features = 'IL6')

# genes
# PC2 <- seurat.stroma.bak@reductions$pca@cell.embeddings[, 'PC_2']
Harmony2 <- seurat.stroma@reductions$harmony@cell.embeddings[, 'harmony_2']
mat.gene <- seurat.stroma@assays$RNA@data
df.pc.gene <- data.frame(t(as.matrix(mat.gene)))
df.pc.gene$Harmony2 <- Harmony2
df.pc.gene$celltype <- seurat.stroma$celltype
df.pc.gene$status <- seurat.stroma$type
df.pc.gene <- df.pc.gene[order(Harmony2, decreasing = F),]
df.pc.gene$idx <- 1:nrow(df.pc.gene)

# 
genes <- c('EGR1', 'COL1A1', 'IL6', 
           'HES1', 'ASPN', 'IL8', 
           'FBLN1', 'VCAN', 'MMP3')
list.plot <- list()
for (gene in genes) {
    # gene <- 'FBLN1'
    df.plot <- df.pc.gene[, c('idx', 'Harmony2', 'celltype', gene)]
    names(df.plot) <- c('idx', 'Harmony2', 'celltype', 'gene')
    p.gene <-
        ggplot(data = df.plot, aes(x = idx, y = gene)) + 
        geom_point(aes(color = celltype), size = 0.3) + 
        scale_color_manual(labels = c('SSC', 'SC1', 'SC2'),
                           values = colors) + 
        # xlim(-25, 8) +
        geom_smooth(color = '#696969') + 
        labs(x = '', y = '') + 
        theme(panel.background=element_rect(fill='transparent', color='black',size = 1),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.position = 'none') +
        annotate('text', label = gene, x = 5800, y = max(df.plot$gene), 
                 hjust = 1, vjust = 1, size = 7)
    ggsave(plot = p.gene, path = path.lineage, 
           filename = paste0(gene, '.png'),
           height = 6, width = 9, units = 'cm')
    list.plot[[gene]] <- p.gene
}
p.merge <- marrangeGrob(list.plot, ncol = 3, nrow = 3, top = NULL)
ggsave(plot = p.merge, path = path.lineage, 
       filename = 'stroma_genes.pdf',
       height = 20, width = 30, units = 'cm')


# cell marker
# clusters <- unique(seurat.stroma$celltype)
# list.marker.all <- list()
# for (cluster in clusters) {
#     sub.markers <- FindMarkers(seurat.stroma, ident.1 = cluster, group.by = 'celltype',
#                                logfc.threshold = 0.3, min.diff.pct = 0.05, only.pos = T)
#     list.marker.all[[cluster]] <- sub.markers
# }
# file.marker <- paste0(path.lineage, 'marker_genes.Rdata')
# saveRDS(list.marker.all, file.marker)
# list.marker.all <- readRDS(file.marker)

# heatmap
sort.cells <- c('SSC', 'SC1', 'SC2')

sel.genes <- c()
for (cell in sort.cells) {
    sub.markers <- list.marker.all[[cell]]
    sub.markers$diff.pct <- sub.markers$pct.1 - sub.markers$pct.2
    sub.markers <- sub.markers[sub.markers$p_val_adj < 0.01 & sub.markers$avg_logFC > 0.5 & 
                                   sub.markers$diff.pct > 0.05, ]
    sub.markers <- sub.markers[order(sub.markers$avg_logFC, decreasing = T),]
    sel.genes <- c(sel.genes, rownames(sub.markers)[1:min(nrow(sub.markers), 50)])
}
sel.genes <- intersect(unique(sel.genes), rownames(seurat.stroma@assays$RNA@scale.data))

seurat.stroma$celltype <- factor(seurat.stroma$celltype, levels = sort.cells)
seurat.stroma.sample <- subset(seurat.stroma,
                             cells = sample(rownames(seurat.stroma@meta.data), 3000))
plot.heatmap <- 
    DoHeatmap(seurat.stroma.sample, features = sel.genes,
              draw.lines = T, lines.width = 20,
              group.by = 'celltype',
              label = T, group.colors = colors, size = 6.5,
              slot = 'scale.data', disp.min = -1.5, disp.max = 1.5) +
    guides(color = F) + 
    labs(fill = 'Scaled Expression') + 
    theme(axis.text = element_blank(),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18),
          legend.position = 'bottom', legend.direction = 'horizontal',
          plot.margin = unit(c(0.5,0,0,0), "cm"))+
    scale_fill_gradientn(colors = c("navy", "white", "firebrick3"), 
                         breaks = c(-1, 0, 1)) 
ggsave(plot = plot.heatmap, path = path.lineage, 
       filename = 'DoHeatmap_markers_sample.png',
       height = 20, width = 20, units = 'cm')


# GO
fc.cutoff <- 0.4
pct.cutoff <- 0.05
clusters <- unique(seurat.stroma$celltype)
list.marker.all <- list()
for (cluster in clusters) {
    sub.markers <- FindMarkers(seurat.stroma, ident.1 = cluster, group.by = 'celltype',
                               logfc.threshold = fc.cutoff, min.diff.pct = pct.cutoff, only.pos = T)
    sub.markers <- sub.markers[sub.markers$p_val_adj < 0.01,]
    list.marker.all[[cluster]] <- sub.markers
}

path.cutoff <- paste0(path.lineage, 'cutoff_', fc.cutoff, '_', pct.cutoff, '/')
if (!file.exists(path.cutoff)) {
    dir.create(path.cutoff)
}
file.marker.all <- paste0(path.cutoff, 'marker_go.Rdata')
saveRDS(list.marker.all, file = file.marker.all)
list.marker.all <- readRDS(file.marker.all)

# enrich GO
library(dplyr)
df.gene_id <- read.delim('/home/disk/drizzle/wgk/ncbi/Homo_sapiens.gene_info')
df.gene_id <- df.gene_id %>% distinct(Symbol, .keep_all = T)
df.symbol2geneid <- data.frame(Gene_id = df.gene_id$GeneID, row.names = df.gene_id$Symbol)
df.geneid2symbol <- data.frame(Symbol = df.gene_id$Symbol, row.names = df.gene_id$GeneID)
library(clusterProfiler)
library(org.Hs.eg.db)
all.genes <- rownames(seurat.stroma@assays$RNA@data)
use.genes <- intersect(keys(org.Hs.eg.db, keytype = "SYMBOL"), all.genes)
list.go.BP <- list()
list.go.MF <- list()
list.go.CC <- list()
list.kegg <- list()
for (cell in clusters) {
    sub.markers <- list.marker.all[[as.character(cell)]]
    genes.input <- intersect(row.names(sub.markers), use.genes)
    egmt <- enrichGO(gene = genes.input, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", 
                     universe = use.genes, pvalueCutoff = 0.1, ont = 'BP')
    # res.egmt <- simplify(egmt)@result
    res.egmt <- egmt@result
    list.go.BP[[as.character(cell)]] <- res.egmt
    # egmt <- enrichGO(gene = genes.input, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", 
    #                  universe = use.genes, pvalueCutoff = 0.5, ont = 'MF')
    # # res.egmt <- simplify(egmt)@result
    # res.egmt <- egmt@result
    # list.go.MF[[as.character(cell)]] <- res.egmt
    # egmt <- enrichGO(gene = genes.input, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", 
    #                  universe = use.genes, pvalueCutoff = 0.5, ont = 'CC')
    # # res.egmt <- simplify(egmt)@result
    # res.egmt <- egmt@result
    # list.go.CC[[as.character(cell)]] <- res.egmt
    # res.kegg <- enrichKEGG(gene = df.symbol2geneid[genes.input, 'Gene_id'], 
    #                        universe = as.character(df.symbol2geneid[use.genes, 'Gene_id']), 
    #                        pvalueCutoff = 0.5, keyType = 'kegg')
    # list.kegg[[as.character(cell)]] <- res.kegg@result
}

file.go.BP <- paste0(path.cutoff, 'GO_BP.Rdata')
saveRDS(list.go.BP, file = file.go.BP)
# file.go.MF <- paste0(path.cutoff, 'GO_MF.Rdata')
# saveRDS(list.go.MF, file = file.go.MF)
# file.go.CC <- paste0(path.cutoff, 'GO_CC.Rdata')
# saveRDS(list.go.CC, file = file.go.CC)
# file.kegg <- paste0(path.cutoff, 'kegg.Rdata')
# saveRDS(list.kegg, file = file.kegg)

list.go.BP <- readRDS(file.go.BP)

list.sel.GO <- list()
list.sel.GO$SSC <- c('transforming growth factor beta receptor signaling pathway', 
                     'extracellular matrix organization', 
                     'connective tissue development',
                     'BMP signaling pathway', 'cell cycle arrest')
list.sel.GO$SC2 <- c('extracellular matrix organization', 
                     'extracellular matrix disassembly',
                     'regulation of immunoglobulin secretion',
                     'positive regulation of cell migration')
list.sel.GO$SC1 <- c('extracellular matrix organization', 
                     'regulation of inflammatory response', 
                     'negative regulation of cell migration', 
                     'chemokine-mediated signaling pathway',
                     'regulation of angiogenesis', 
                     'positive regulation of insulin-like growth factor receptor signaling pathway',
                     'regulation of smooth muscle cell proliferation')
# barplot
sort.cells <- c('SSC', 'SC1', 'SC2')
colors <- c("#BC80BD", "#80B1D3", "#F4A460")
df.plot <- data.frame()
i = 0
for (cell in sort.cells) {
    i = i + 1
    sub.go <- list.go.BP[[cell]]
    sel.go.term <- list.sel.GO[[cell]]
    sel.go <- sub.go[sub.go$Description %in% sel.go.term, 
                     c('Description', 'pvalue')]
    sel.go$log10Pval <- -log10(sel.go$pvalue)
    sel.go$celltype <- rep(cell, nrow(sel.go))
    # sel.go$Description <- factor(sel.go$Description, levels = rev(sel.go.term))
    df.plot <- rbind(df.plot, sel.go)
}
col_name <- paste(df.plot$celltype, df.plot$Description, sep = '_')
df.plot$col_name <- factor(col_name, levels = rev(col_name))
df.plot$celltype <- factor(df.plot$celltype, levels = sort.cells)

p <- ggplot(df.plot, aes(x = celltype, y = col_name, 
                         color = celltype, size = log10Pval)) + 
    geom_point(fill = 'cornsilk') +
    scale_color_manual(breaks = c('SSC', 'SC1', 'SC2'),
                       values = c("#BC80BD", "#80B1D3", "#F4A460")) + 
    scale_size_continuous(range = c(3,5)) +
    scale_y_discrete(breaks = col_name,
                     labels = c(list.sel.GO$SSC, list.sel.GO$SC1, list.sel.GO$SC2)) +
    labs(x = '', y = 'GO term', color = 'Cell type',
         size = expression(paste("-log"[10], "(", italic("P"), "-value)"))) +
    theme(panel.background = element_rect(color = 'black',
                                          fill = 'transparent'), 
          panel.grid.major = element_line(colour = 'gray', size = 0.2, linetype = 5),
          axis.title = element_text(size = 14, face = 'bold', color = 'black'), 
          axis.text.y = element_text(size = 12, face = 'bold', color = 'black'), 
          axis.text.x = element_text(size = 12, face = 'bold', color = 'black'),
          legend.text = element_text(size = 12, color = 'black'),
          legend.title = element_text(size = 14, face = 'bold', color = 'black'),
          legend.key = element_blank()) + 
    guides(colour = guide_legend(override.aes = list(size=4)))
ggsave(plot = p, path = path.cutoff, 
       filename = paste0('GO.pdf'),
       height = 12, width = 24, units = 'cm')


