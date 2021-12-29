setwd('/home/zy/my_git/bioinformatics/wgk')
.libPaths(c('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0'))
library(Seurat)

path.data <- '/home/disk/drizzle/wgk/data/AllSample_2_merge/'
path.lineage <- paste0(path.data, 'chon_lineage/')
file.chon <- paste0(path.lineage, 'seurat_celltype.Rdata')
seurat.chon <- readRDS(file.chon)

DimPlot(seurat.chon, group.by = "batch")
DimPlot(seurat.chon, group.by = "type")
DimPlot(seurat.chon, group.by = "celltype", label = T)

# M1 M2 M3
seurat.child <- subset(seurat.chon, subset = batch %in% c('C4', 'C6', 'M1', 'M2', 'M3'))
# old.cells <- seurat.child$celltype
# new.cells <- as.character(old.cells)
# new.cells[old.cells %in% c('Chondrocyte1', 'Chondrocyte2')] <- 'Chondrocyte'
# seurat.child$celltype <- new.cells
DimPlot(seurat.child, group.by = "celltype", label = T)
path.M123 <- '/home/disk/drizzle/wgk/microtia_chon_child_M1M2M3/'
if (!file.exists(path.M123)) {
    dir.create(path.M123)
}

# marker gene
fc.cutoff <- 0.4
pct.cutoff <- 0
celltypes <- unique(seurat.child$celltype)
list.marker.go <- list()
for (cell in celltypes) {
    sub.seurat <- subset(seurat.child, subset = celltype == cell)
    sub.markers <- FindMarkers(sub.seurat, ident.1 = 'Microtia', group.by = 'type',
                               logfc.threshold = fc.cutoff, min.diff.pct = pct.cutoff)
    sub.markers <- sub.markers[sub.markers$p_val_adj < 0.01,]
    list.marker.go[[as.character(cell)]] <- sub.markers
}
# FeaturePlot(sub.seurat, features = c('FRZB'))

path.cutoff <- paste0(path.M123, 'cutoff_', fc.cutoff, '/')
if (!file.exists(path.cutoff)) {
    dir.create(path.cutoff)
}
file.marker.go <- paste0(path.cutoff, 'marker_go.Rdata')
saveRDS(list.marker.go, file = file.marker.go)
list.marker.go <- readRDS(file.marker.go)

# enrich GO
library(dplyr)
df.gene_id <- read.delim('/home/disk/drizzle/wgk/ncbi/Homo_sapiens.gene_info')
df.gene_id <- df.gene_id %>% distinct(Symbol, .keep_all = T)
df.symbol2geneid <- data.frame(Gene_id = df.gene_id$GeneID, row.names = df.gene_id$Symbol)
df.geneid2symbol <- data.frame(Symbol = df.gene_id$Symbol, row.names = df.gene_id$GeneID)
library(clusterProfiler)
library(org.Hs.eg.db)
all.genes <- rownames(seurat.child@assays$RNA@data)
use.genes <- intersect(keys(org.Hs.eg.db, keytype = "SYMBOL"), all.genes)
status <- c('Microtia_increase', 'Microtia_decrease')
list.go.BP <- list()
list.go.MF <- list()
list.go.CC <- list()
list.kegg <- list()
list.go.BP.all <- list()
list.go.MF.all <- list()
list.go.CC.all <- list()
celltypes <- names(list.marker.go)
for (cell in celltypes) {
    for (st in status) {
        type <- paste0(as.character(cell), '_', st)
        sub.markers <- list.marker.go[[as.character(cell)]]
        if (st == 'Microtia_increase') {
            genes <- row.names(sub.markers[sub.markers$avg_logFC > 0,])
        } else {
            genes <- row.names(sub.markers[sub.markers$avg_logFC < 0,])
        }
        genes.input <- intersect(genes, use.genes)
        egmt <- enrichGO(gene = genes.input, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", 
                         universe = use.genes, pvalueCutoff = 0.1, ont = 'BP')
        res.egmt <- egmt@result
        list.go.BP.all[[as.character(type)]] <- res.egmt
        res.egmt <- simplify(egmt)@result
        list.go.BP[[as.character(type)]] <- res.egmt
        egmt <- enrichGO(gene = genes.input, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                         universe = use.genes, pvalueCutoff = 0.5, ont = 'MF')
        res.egmt <- egmt@result
        list.go.MF.all[[as.character(type)]] <- res.egmt
        res.egmt <- simplify(egmt)@result
        list.go.MF[[as.character(type)]] <- res.egmt
        egmt <- enrichGO(gene = genes.input, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                         universe = use.genes, pvalueCutoff = 0.5, ont = 'CC')
        list.go.CC.all[[as.character(type)]] <- egmt@result
        res.egmt <- simplify(egmt)@result
        list.go.CC[[as.character(type)]] <- res.egmt
        res.kegg <- enrichKEGG(gene = df.symbol2geneid[genes.input, 'Gene_id'],
                               universe = as.character(df.symbol2geneid[use.genes, 'Gene_id']),
                               pvalueCutoff = 0.5, keyType = 'kegg')
        list.kegg[[as.character(type)]] <- res.kegg@result
    }
}

file.go.BP.all <- paste0(path.cutoff, 'GO_BP_all.Rdata')
saveRDS(list.go.BP.all, file = file.go.BP.all)
file.go.MF.all <- paste0(path.cutoff, 'GO_MF_all.Rdata')
saveRDS(list.go.MF.all, file = file.go.MF.all)
file.go.CC.all <- paste0(path.cutoff, 'GO_CC_all.Rdata')
saveRDS(list.go.CC.all, file = file.go.CC.all)

file.go.BP <- paste0(path.cutoff, 'GO_BP.Rdata')
saveRDS(list.go.BP, file = file.go.BP)
file.go.MF <- paste0(path.cutoff, 'GO_MF.Rdata')
saveRDS(list.go.MF, file = file.go.MF)
file.go.CC <- paste0(path.cutoff, 'GO_CC.Rdata')
saveRDS(list.go.CC, file = file.go.CC)
file.kegg <- paste0(path.cutoff, 'kegg.Rdata')
saveRDS(list.kegg, file = file.kegg)


# GSEA marker
celltypes <- unique(seurat.child$celltype)
list.marker.gsea <- list()
for (cell in celltypes) {
    sub.seurat <- subset(seurat.child, subset = celltype == cell)
    sub.markers <- FindMarkers(sub.seurat, ident.1 = 'Microtia', group.by = 'type',
                               logfc.threshold = 0, min.diff.pct = 0)
    list.marker.gsea[[as.character(cell)]] <- sub.markers
}

file.gsea.marker <- paste0(path.M123, 'marker_GSEA.Rdata')
saveRDS(list.marker.gsea, file = file.gsea.marker)
list.marker.gsea <- readRDS(file.gsea.marker)

library(clusterProfiler)
library(org.Hs.eg.db)
df.gene_id <- read.delim('/home/disk/drizzle/wgk/ncbi/Homo_sapiens.gene_info')
df.gene_id <- df.gene_id %>% distinct(Symbol, .keep_all = T)
df.symbol2geneid <- data.frame(Gene_id = df.gene_id$GeneID, row.names = df.gene_id$Symbol)
df.geneid2symbol <- data.frame(Symbol = df.gene_id$Symbol, row.names = df.gene_id$GeneID)
# use.genes <- intersect(keys(org.Hs.eg.db, keytype = "SYMBOL"), all.genes)
list.gsea <- list()
list.gsea.kegg <- list()
celltypes <- unique(seurat.child$celltype)
for (cell in celltypes) {
    sub.markers <- list.marker.gsea[[as.character(cell)]]
    geneList <- sub.markers[, 'avg_logFC']
    names(geneList) <- rownames(sub.markers)
    geneList <- geneList[order(geneList, decreasing = T)]
    # GO
    egmt <- gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", pvalueCutoff = 0.05)
    res.egmt <- egmt@result
    list.gsea[[as.character(cell)]] <- res.egmt
    # KEGG
    names(geneList) <- df.symbol2geneid[rownames(sub.markers), 'Gene_id']
    res.kegg <- gseKEGG(geneList = geneList, pvalueCutoff = 0.05)
    list.gsea.kegg[[cell]] <- res.kegg@result
}

file.gsea.go <- paste0(path.M123, 'GSEA_GO.Rdata')
saveRDS(list.gsea, file = file.gsea.go)
file.gsea.kegg <- paste0(path.M123, 'GSEA_KEGG.Rdata')
saveRDS(list.gsea.kegg, file = file.gsea.kegg)





