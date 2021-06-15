setwd('/home/zy/my_git/bioinformatics/wgk')
.libPaths('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)
library(ggplot2)
require("RColorBrewer")

# path.data <- '/home/disk/drizzle/wgk/data/AllSample_2_merge/'
# file.merge_2 <- paste0(path.data, 'seurat_celltype.Rdata')
# seurat.all.filter <- readRDS(file.merge_2)
# 
# path.out <- paste0(path.data, 'GO_celltype/')
path.out <- path.fig

# cell marker
fc.cutoff <- 0.5
pct.cutoff <- 0.1
clusters <- unique(seurat.all.filter$celltype.abbr)
list.marker.all <- list()
for (cluster in clusters) {
    sub.markers <- FindMarkers(seurat.all.filter, ident.1 = cluster, group.by = 'celltype.abbr',
                               logfc.threshold = fc.cutoff, min.diff.pct = pct.cutoff, only.pos = T)
    sub.markers <- sub.markers[sub.markers$p_val_adj < 0.01,]
    list.marker.all[[cluster]] <- sub.markers
}

path.cutoff <- paste0(path.out, 'cutoff_', fc.cutoff, '_', pct.cutoff, '/')
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
all.genes <- rownames(seurat.all.filter@assays$RNA@data)
use.genes <- intersect(keys(org.Hs.eg.db, keytype = "SYMBOL"), all.genes)
list.go.BP <- list()
list.go.MF <- list()
list.go.CC <- list()
list.kegg <- list()
for (cell in clusters) {
        sub.markers <- list.marker.all[[as.character(cell)]]
        genes.input <- intersect(row.names(sub.markers), use.genes)
        egmt <- enrichGO(gene = genes.input, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", 
                         universe = use.genes, pvalueCutoff = 0.5, ont = 'BP')
        # res.egmt <- simplify(egmt)@result
        res.egmt <- egmt@result
        list.go.BP[[as.character(cell)]] <- res.egmt
        egmt <- enrichGO(gene = genes.input, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", 
                         universe = use.genes, pvalueCutoff = 0.5, ont = 'MF')
        # res.egmt <- simplify(egmt)@result
        res.egmt <- egmt@result
        list.go.MF[[as.character(cell)]] <- res.egmt
        egmt <- enrichGO(gene = genes.input, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", 
                         universe = use.genes, pvalueCutoff = 0.5, ont = 'CC')
        # res.egmt <- simplify(egmt)@result
        res.egmt <- egmt@result
        list.go.CC[[as.character(cell)]] <- res.egmt
        res.kegg <- enrichKEGG(gene = df.symbol2geneid[genes.input, 'Gene_id'], 
                               universe = as.character(df.symbol2geneid[use.genes, 'Gene_id']), 
                               pvalueCutoff = 0.5, keyType = 'kegg')
        list.kegg[[as.character(cell)]] <- res.kegg@result
}

file.go.BP <- paste0(path.cutoff, 'GO_BP.Rdata')
saveRDS(list.go.BP, file = file.go.BP)
file.go.MF <- paste0(path.cutoff, 'GO_MF.Rdata')
saveRDS(list.go.MF, file = file.go.MF)
file.go.CC <- paste0(path.cutoff, 'GO_CC.Rdata')
saveRDS(list.go.CC, file = file.go.CC)
file.kegg <- paste0(path.cutoff, 'kegg.Rdata')
saveRDS(list.kegg, file = file.kegg)

list.go.BP <- readRDS(file.go.BP)


# plot
# main figure
GO.BP <- c('granulocyte activation',
           'leukocyte migration',
           'antigen processing and presentation',
           'angiogenesis',
           'endothelium development',
           'muscle system process',
           'angiogenesis',
           'extracellular matrix organization',
           'regulation of smooth muscle cell proliferation',
           'regulation of inflammatory response',
           'regulation of angiogenesis',
           'extracellular matrix organization',
           'cellular response to transforming growth factor beta stimulus',
           'connective tissue development',
           'BMP signaling pathway',
           'cell cycle arrest',
           'extracellular matrix organization',
           'chondrocyte differentiation',
           'response to BMP',
           'cell cycle arrest',
           'cellular response to transforming growth factor beta stimulus',
           'response to metal ion',
           'extracellular structure organization',
           'positive regulation of cell-substrate adhesion')
GO.cell <- c('IC', 'IC', 'IC', 'EC', 'EC', 'PVC', 'PVC', 
             'SC', 'SC', 'SC', 'SC', 
             'SSC', 'SSC', 'SSC', 'SSC', 'SSC', 
             'CSC', 'CSC', 'CSC', 'CSC', 'CSC',  
             'C', 'C', 'C')

sort.cells <- c('IC', 'EC', 'PVC', 'SC', 'SSC', 'CSC', 'C')
color.cell <- c("#E31A1C","#6A3D9A","#FF7F00" ,"#33A02C","#CAB2D6" ,"#A6CEE3","#1F78B4")

df.plot <- data.frame(stringsAsFactors = F)
for (i in 1:length(GO.BP)) {
    GO.term <- GO.BP[i]
    cell <- GO.cell[i]
    sub.BP <- list.go.BP[[cell]]
    df.plot[i, 'GO'] <- GO.term
    df.plot[i, 'Cell'] <- cell
    df.plot[i, 'pvalue'] <- sub.BP[
        sub.BP$Description == GO.term, 'pvalue']
    df.plot[i, 'name'] <- paste0(cell, '_', GO.term)
}
df.plot$log10Pval <- -log10(df.plot$pvalue)
df.plot$Cell <- factor(df.plot$Cell, levels = sort.cells)
df.plot$name <- factor(df.plot$name, levels = rev(paste(GO.cell, GO.BP, sep = '_')))
df.plot$log10Pval[df.plot$log10Pval > 20] <- 20

plot.bar <- 
    ggplot(data = df.plot, aes(x = name, y = log10Pval, color = Cell, fill = Cell)) + 
    geom_bar(stat = 'identity') + 
    theme_classic() + coord_flip() +
    scale_x_discrete(breaks = paste(GO.cell, GO.BP, sep = '_'), 
                     labels = GO.BP,
                     position = "top") + 
    scale_color_manual(values = color.cell) +
    scale_fill_manual(values = color.cell) +
    labs(x = 'GO term', y = expression(paste("-log"[10], "(", italic("P"), "-value)"))) +
    theme(axis.title = element_text(size = 23, face = 'bold', color = 'black'), 
          axis.text.y = element_text(size = 20, face = 'bold', color = 'black'), 
          axis.text.x = element_text(size = 20, face = 'bold', color = 'black'), 
          legend.position = 'none')
ggsave(plot = plot.bar, path = path.cutoff, 
       filename = 'heatmap_bar.png',
       height = 30, width = 30, units = 'cm')
ggsave(plot = plot.bar, path = path.cutoff, 
       filename = 'heatmap_bar.pdf',
       height = 30, width = 30, units = 'cm')


# supp figure
sel.GO.BP <- c('extracellular matrix organization',
               'extracellular structure organization',
               'cellular response to zinc ion',
               'cellular response to copper ion',
               'collagen fibril organization',
               'cartilage condensation',
               'cartilage development',
               'chondrocyte differentiation',
               'connective tissue development',
               'response to BMP',
               'BMP signaling pathway',
               'regulation of cyclin-dependent protein kinase activity',
               'response to transforming growth factor beta',
               'Wnt signaling pathway',
               'Notch signaling pathway',
               'ear development',
               'cell cycle arrest',
               'regulation of angiogenesis',
               'regulation of inflammatory response',
               'regulation of smooth muscle cell proliferation',
               'regulation of epithelial cell proliferation',
               'extracellular matrix disassembly',
               'regulation of epithelial cell migration',
               'regulation of endothelial cell migration',
               'regulation of smooth muscle cell migration',
               'angiogenesis', 'vasculogenesis',
               'endothelium development',
               'endothelial cell differentiation',
               'endothelial cell migration',
               'muscle system process',
               'muscle cell proliferation',
               'muscle contraction',
               'neutrophil activation',
               'leukocyte chemotaxis',
               'neutrophil chemotaxis',
               'antigen processing and presentation',
               'T cell activation')

celltypes <- setdiff(names(list.go.BP), 'Chondrocyte2')
df.plot <- data.frame(stringsAsFactors = F)
for (cell in celltypes) {
    sub.BP <- list.go.BP[[cell]]
    rownames(sub.BP) <- sub.BP$Description
    sub.BP <- sub.BP[sub.BP$p.adjust < 0.1,]
    sel.BP <- sub.BP[sel.GO.BP, c('Description', 'pvalue')]
    sel.BP$Description <- sel.GO.BP
    sel.BP$pvalue[is.na(sel.BP$pvalue)] <- 1
    sel.BP$CellType <- rep(cell, nrow(sel.BP))
    df.plot <- rbind(df.plot, sel.BP)
}
df.plot$log10Pval <- -log10(df.plot$pvalue)
mat.plot <- reshape2::dcast(df.plot, Description ~ CellType, value.var = 'log10Pval')
# write mat
file.mat <- paste0(path.cutoff, 'GO_mat.txt')
write.table(mat.plot, file.mat, sep = '\t')
mat.plot <- read.delim(file.mat, sep = '\t', check.names = F)

df.plot.2 <- reshape2::melt(mat.plot)
names(df.plot.2) <- c('Description', 'CellType', 'log10Pval')
df.plot.2$Description <- factor(df.plot.2$Description,
                                levels = sel.GO.BP)
df.plot.2$CellType <- factor(df.plot.2$CellType,
                             levels = c("Chondrocyte1", "Transitional chondrocyte", "Chondral stem cell",
                                        "Stromal stem cell", "Stromal cell1", "Stromal cell2",
                                        "Endothelial cell", "Perivascular cell", "Immune cell"))

plot.bar <- 
    ggplot(data = df.plot.2, aes(x = Description, y = log10Pval, color = CellType, fill = CellType)) + 
    geom_bar(stat = 'identity') + 
    facet_grid( ~ CellType, scales = "free_x") + 
    theme_classic() + coord_flip() +
    labs(x = 'GO term', y = expression(paste("-log"[10], "(P value)"))) +
    theme(legend.position = 'none')

ggsave(plot = plot.bar, path = path.cutoff, 
       filename = 'bar.png',
       height = 20, width = 40, units = 'cm')


