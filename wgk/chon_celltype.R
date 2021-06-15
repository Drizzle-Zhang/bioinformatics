setwd('/home/zy/my_git/bioinformatics/wgk')
.libPaths('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)
library(harmony)
library(ggplot2)
require("RColorBrewer")

path.data <- '/home/disk/drizzle/wgk/data/AllSample_2_merge/'
file.merge_2 <- paste0(path.data, 'seurat_celltype.Rdata')
seurat.all.filter <- readRDS(file.merge_2)

path.lineage <- paste0(path.data, 'chon_lineage/')

# seurat.chon.bak <- subset(seurat.all.filter, subset = celltype %in% 
#                           c('Chondrocyte1', 'Chondrocyte2', 'Chondral stem cell', 
#                             'Transitional chondrocyte'))
seurat.chon <- subset(seurat.all.filter, subset = celltype %in% 
                          c('Chondrocyte1', 'Chondrocyte2', 'Chondral stem cell', 
                            'Transitional chondrocyte'))
seurat.chon <- NormalizeData(seurat.chon)
seurat.chon <- FindVariableFeatures(seurat.chon, nfeatures = 3000)
seurat.chon <- ScaleData(seurat.chon, split.by = "batch")
seurat.chon <- RunPCA(seurat.chon, verbose = F, npcs = 100)
seurat.chon <- RunHarmony(seurat.chon, "batch", reduction.save = "harmony2")
seurat.chon <- RunUMAP(seurat.chon, reduction = "harmony2", 
                       dims = 1:80, n_neighbors = 50, min.dist = 1)
DimPlot(seurat.chon, group.by = "celltype.abbr", label = T, reduction = 'umap', pt.size = 1)

seurat.chon <- FindNeighbors(seurat.chon, reduction = "harmony", dims = 1:80)
seurat.chon <- FindClusters(seurat.chon, resolution = 0.3)
DimPlot(seurat.chon, group.by = "RNA_snn_res.0.3", label = T, reduction = 'umap', pt.size = 1)


# color <- c('#531F6F', '#5483C7', '#8FBA3B', '#A83B54')
cluster_0.3 <- seurat.chon$RNA_snn_res.0.3
celltypes <- rep('_', length(cluster_0.3))
celltypes[cluster_0.3 %in% c(3)] <- 'CSC'
celltypes[cluster_0.3 %in% c(2)] <- 'TC'
celltypes[cluster_0.3 %in% c(1)] <- 'C1'
celltypes[cluster_0.3 %in% c(0)] <- 'C2'
seurat.chon$celltype <- celltypes
DimPlot(seurat.chon, group.by = "celltype", label = T, pt.size = 1)

seurat.chon$celltype <- factor(seurat.chon$celltype, 
                                    levels = c('CSC', 'TC', 'C1', 'C2'))
# colors <- c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
colors <- c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
names(colors) <- c('CSC', 'TC', 'C1', 'C2')

plot.umap <-
    DimPlot(seurat.chon, reduction = 'umap', 
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
       filename = 'chon_lineage_umap.pdf',
       height = 10, width = 15, units = 'cm')


plot.cluster <-
    DimPlot(seurat.chon, reduction = 'umap', 
            group.by = "RNA_snn_res.0.3", label = T,
            pt.size = 0.8, label.size = 5.5,
            cols = rev(colors)) + 
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
       filename = 'chon_lineage_umap_cluster.pdf',
       height = 10, width = 15, units = 'cm')


file.chon <- paste0(path.lineage, 'seurat_celltype.Rdata')
saveRDS(seurat.chon, file.chon)
seurat.chon <- readRDS(file.chon)

pc.stdev <- seurat.chon@reductions$harmony2@stdev[1:10]
pc.stdev[1]/sum(pc.stdev)
pc.stdev[2]/sum(pc.stdev)

# seurat.chon.bak$celltype <- seurat.chon$celltype

plot.lineage <- 
    DimPlot(seurat.chon, reduction = 'harmony', dims = c(2, 1),
            group.by = "celltype", label = T,
            pt.size = 0.8, label.size = 5.5,
            cols = colors) + 
    labs(x = 'PC2 (15%)', y = 'PC1 (17%)') + 
    theme(axis.title = element_text(size = 16),
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          axis.line = element_line(arrow = arrow(length = unit(0.5, 'cm'))),
          legend.position = 'none')
ggsave(plot.lineage, path = path.lineage, 
       filename = 'chon_lineage_harmony.pdf',
       height = 10, width = 15, units = 'cm')

# genes
# PC2 <- seurat.chon.bak@reductions$pca@cell.embeddings[, 'PC_2']
Harmony2 <- seurat.chon@reductions$harmony@cell.embeddings[, 'harmony_2']
mat.gene <- seurat.chon@assays$RNA@data
df.pc.gene <- data.frame(t(as.matrix(mat.gene)))
df.pc.gene$Harmony2 <- Harmony2
df.pc.gene$celltype <- seurat.chon$celltype
df.pc.gene$status <- seurat.chon$type

genes <- c('EGR1', 'HES1', 'VIT', 
           # 'FRZB', 'CTGF', 
           'CYTL1', 'COL2A1', 
           'COL9A2', 'ACAN', 
           'IL8', 'MMP3')
genes <- c('EGR1', 'FRZB', 'ACAN', 
           'HES1', 'COL2A1', 'IL8', 
           'CTGF', 'COL9A2', 'MMP3')
list.plot <- list()
for (gene in genes) {
    df.plot <- df.pc.gene[, c('Harmony2', 'celltype', gene)]
    names(df.plot) <- c('Harmony2', 'celltype', 'gene')
    p.gene <- 
        ggplot(data = df.plot, aes(x = Harmony2, y = gene)) + 
        geom_point(aes(color = celltype), size = 0.3) + 
        scale_color_manual(labels = c('CSC', 'TC', 'C1', 'C2'),
                           values = colors) + 
        xlim(-30, 10) + 
        geom_smooth(color = '#696969') + 
        labs(x = '', y = '') + 
        theme(panel.background=element_rect(fill='transparent', color='black',size = 1),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.position = 'none') +
        annotate('text', label = gene, x = 10, y = max(df.plot$gene), 
                 hjust = 1, vjust = 1, size = 7)
    list.plot[[gene]] <- p.gene
}

library(gridExtra)
p.merge <- marrangeGrob(list.plot, ncol = 3, nrow = 3, top = NULL)
ggsave(plot = p.merge, path = path.lineage, 
       filename = 'chon_genes.pdf',
       height = 20, width = 30, units = 'cm')

# supp
genes <- c('ING1', 'CYTL1','ICAM1', 
           'CCNL1', 'SPARC', 'IL6', 
           'VIT', 'COL11A1', 'MMP1',
           'SCRG1', 'COL9A3', 'CXCL3',
           'ELN', 'MT1X', 'CXCL2')
list.plot <- list()
for (gene in genes) {
    df.plot <- df.pc.gene[, c('Harmony2', 'celltype', gene)]
    names(df.plot) <- c('Harmony2', 'celltype', 'gene')
    p.gene <- 
        ggplot(data = df.plot, aes(x = Harmony2, y = gene)) + 
        geom_point(aes(color = celltype), size = 0.15) + 
        scale_color_manual(labels = c('CSC', 'TC', 'C1', 'C2'),
                           values = colors) + 
        xlim(-30, 10) + 
        geom_smooth(color = '#696969') + 
        labs(x = '', y = '') + 
        theme_bw() + 
        theme(panel.background=element_rect(fill='transparent', color='black',size = 1),
              panel.grid  = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.position = 'none') +
        annotate('text', label = gene, x = 10, y = max(df.plot$gene), 
                 hjust = 1, vjust = 1, size = 5, fontface = 'italic')
    list.plot[[gene]] <- p.gene
}

library(gridExtra)
p.merge <- marrangeGrob(list.plot, ncol = 5, nrow = 3, top = NULL)
ggsave(plot = p.merge, path = path.lineage, 
       filename = 'chon_genes_supp.png',
       height = 12, width = 30, units = 'cm')
ggsave(plot = p.merge, path = path.lineage, 
       filename = 'chon_genes_supp.pdf',
       height = 12, width = 30, units = 'cm')

# cell marker
clusters <- unique(seurat.chon$celltype)
list.marker.all <- list()
for (cluster in clusters) {
    sub.markers <- FindMarkers(seurat.chon, ident.1 = cluster, group.by = 'celltype',
                               logfc.threshold = 0.3, min.diff.pct = 0.05, only.pos = T)
    list.marker.all[[cluster]] <- sub.markers
}
file.marker <- paste0(path.lineage, 'marker_genes.Rdata')
saveRDS(list.marker.all, file.marker)
list.marker.all <- readRDS(file.marker)

# heatmap
sort.cells <- c('CSC', 'TC', 'C1', 'C2')

sel.genes <- c()
for (cell in sort.cells) {
    sub.markers <- list.marker.all[[cell]]
    sub.markers$diff.pct <- sub.markers$pct.1 - sub.markers$pct.2
    sub.markers <- sub.markers[sub.markers$p_val_adj < 0.01 & sub.markers$avg_logFC > 0.5 & 
                                   sub.markers$diff.pct > 0.05, ]
    sub.markers <- sub.markers[order(sub.markers$avg_logFC, decreasing = T),]
    sel.genes <- c(sel.genes, rownames(sub.markers)[1:min(nrow(sub.markers), 50)])
}
sel.genes <- intersect(unique(sel.genes), rownames(seurat.chon@assays$RNA@scale.data))

seurat.chon$celltype <- factor(seurat.chon$celltype, levels = sort.cells)
seurat.chon.sample <- subset(seurat.chon, 
                             cells = sample(rownames(seurat.chon@meta.data), 5000))
plot.heatmap <- 
    DoHeatmap(seurat.chon.sample, features = sel.genes,
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
       filename = 'DoHeatmap_markers_sample.pdf',
       height = 20, width = 20, units = 'cm')

# GO
fc.cutoff <- 0.3
pct.cutoff <- 0.05
clusters <- unique(seurat.chon$celltype)
list.marker.all <- list()
for (cluster in clusters) {
    sub.markers <- FindMarkers(seurat.chon, ident.1 = cluster, group.by = 'celltype',
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
all.genes <- rownames(seurat.chon@assays$RNA@data)
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
list.sel.GO$CSC <- c('cartilage development', 'chondrocyte differentiation',
                     'extracellular matrix organization', 
                     'BMP signaling pathway', 'cell cycle arrest', 
                     'Wnt signaling pathway', 'Notch signaling pathway')
list.sel.GO$TC <- c('cellular response to zinc ion', 'extracellular matrix organization',
                    'cellular response to copper ion', 'cartilage development',
                    'chondrocyte differentiation')
list.sel.GO$C1 <- c('extracellular matrix organization', 'nitric oxide biosynthetic process',
                    'cellular response to interferon-gamma', 'lymphocyte migration',
                    'cellular response to tumor necrosis factor', 
                    'regulation of cell shape')
list.sel.GO$C2 <- c('interleukin-1-mediated signaling pathway', 'neutrophil migration', 
                    'response to oxidative stress', 'extracellular matrix organization',
                    'nitric oxide biosynthetic process', 'extracellular matrix disassembly')
# barplot
sort.cells <- c('CSC', 'TC', 'C1', 'C2')
colors <- c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
i = 0
for (cell in sort.cells) {
    i = i + 1
    sub.go <- list.go.BP[[cell]]
    sel.go.term <- list.sel.GO[[cell]]
    sel.go <- sub.go[sub.go$Description %in% sel.go.term, 
                     c('Description', 'pvalue')]
    sel.go$log10Pval <- -log10(sel.go$pvalue)
    sel.go$Description <- factor(sel.go$Description, levels = rev(sel.go.term))
    p <- ggplot(sel.go, aes(x = Description, y = log10Pval)) + 
        geom_bar(stat = 'identity', color = colors[i], fill = colors[i]) + 
        theme_classic() + coord_flip() +
        labs(x = 'GO term', y = expression(paste("-log"[10], "(", italic("P"), "-value)"))) +
        theme(axis.title = element_text(size = 16, face = 'bold'), 
              axis.text.y = element_text(size = 16, face = 'bold'), 
              axis.text.x = element_text(size = 16, face = 'bold'))
    ggsave(plot = p, path = path.cutoff, 
           filename = paste0(cell, '_bar.pdf'),
           height = 8, width = 18, units = 'cm')
}













