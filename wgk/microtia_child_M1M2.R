.libPaths('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)

path.data <- '/home/disk/drizzle/wgk/data/AllSample_2_merge/'
file.all <- paste0(path.data, 'seurat_celltype.Rdata')
seurat.all <- readRDS(file.all)

DimPlot(seurat.all, group.by = "batch")
DimPlot(seurat.all, group.by = "type")
DimPlot(seurat.all, group.by = "celltype", label = T)

# M1 M2
seurat.cs <- subset(seurat.all, subset = celltype %in% 
                          c('Chondrocyte1', 'Chondrocyte2', 'Chondral stem cell', 
                            'Transitional chondrocyte', 'Stromal stem cell',
                            'Stromal cell1', 'Stromal cell2') & batch %in% c('C4', 'C6', 'M1', 'M2'))
DimPlot(seurat.cs, group.by = "celltype", label = T)
path.M12 <- '/home/disk/drizzle/wgk/microtia_child_M1M2/'
if (!file.exists(path.M12)) {
    dir.create(path.M12)
}

# marker gene
fc.cutoff <- 0.7
pct.cutoof <- 0
celltypes <- unique(seurat.cs$celltype)
list.marker.go <- list()
for (cell in celltypes) {
    sub.seurat <- subset(seurat.cs, subset = celltype == cell)
    sub.markers <- FindMarkers(sub.seurat, ident.1 = 'Microtia', group.by = 'type',
                               logfc.threshold = fc.cutoff, min.diff.pct = pct.cutoof)
    sub.markers <- sub.markers[sub.markers$p_val_adj < 0.01,]
    list.marker.go[[as.character(cell)]] <- sub.markers
}
# FeaturePlot(sub.seurat, features = c('FRZB'))

path.cutoff <- paste0(path.M12, 'cutoff_', fc.cutoff, '/')
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
all.genes <- rownames(seurat.cs@assays$RNA@data)
use.genes <- intersect(keys(org.Hs.eg.db, keytype = "SYMBOL"), all.genes)
status <- c('Microtia_increase', 'Microtia_decrease')
list.go.BP <- list()
list.go.MF <- list()
list.go.CC <- list()
list.kegg <- list()
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
                         universe = use.genes, pvalueCutoff = 0.5, ont = 'BP')
        res.egmt <- simplify(egmt)@result
        list.go.BP[[as.character(type)]] <- res.egmt
        egmt <- enrichGO(gene = genes.input, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", 
                         universe = use.genes, pvalueCutoff = 0.5, ont = 'MF')
        res.egmt <- simplify(egmt)@result
        list.go.MF[[as.character(type)]] <- res.egmt
        egmt <- enrichGO(gene = genes.input, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", 
                         universe = use.genes, pvalueCutoff = 0.5, ont = 'CC')
        res.egmt <- simplify(egmt)@result
        list.go.CC[[as.character(type)]] <- res.egmt
        res.kegg <- enrichKEGG(gene = df.symbol2geneid[genes.input, 'Gene_id'], 
                               universe = as.character(df.symbol2geneid[use.genes, 'Gene_id']), 
                               pvalueCutoff = 0.5, keyType = 'kegg')
        list.kegg[[as.character(type)]] <- res.kegg@result
    }
}

file.go.BP <- paste0(path.cutoff, 'GO_BP.Rdata')
saveRDS(list.go.BP, file = file.go.BP)
file.go.MF <- paste0(path.cutoff, 'GO_MF.Rdata')
saveRDS(list.go.MF, file = file.go.MF)
file.go.CC <- paste0(path.cutoff, 'GO_CC.Rdata')
saveRDS(list.go.CC, file = file.go.CC)
file.kegg <- paste0(path.cutoff, 'kegg.Rdata')
saveRDS(list.kegg, file = file.kegg)



