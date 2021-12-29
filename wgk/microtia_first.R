library(Seurat)

path.data <- '/home/disk/drizzle/wgk/data/marker_1.5_merge/'
file.merge_1.5 <- paste0(path.data, 'seurat_first.Rdata')
seurat.first <- readRDS(file.merge_1.5)

DimPlot(seurat.first, group.by = "status")

DimPlot(seurat.first, group.by = "celltype", label = T)

seurat.chon <- subset(seurat.first, subset = celltype %in% 
                          c('Chondrocyte1', 'Chondrocyte2', 'CSPC', 
                            'Stromal cell1', 'Stromal cell2', 'Pre Chondrocyte'))
DimPlot(seurat.chon, group.by = "celltype", label = T)

# marker gene
fc.cutoff <- 0.5
pct.cutoof <- 0
celltypes <- unique(seurat.chon$celltype)
list.marker.go <- list()
for (cell in celltypes) {
    sub.seurat <- subset(seurat.chon, subset = celltype == cell)
    sub.markers <- FindMarkers(sub.seurat, ident.1 = 'Microtia', group.by = 'status',
                               logfc.threshold = fc.cutoff, min.diff.pct = pct.cutoof)
    list.marker.go[[as.character(cell)]] <- sub.markers
}
# FeaturePlot(sub.seurat, features = c('FRZB'))

file.marker.go <- '/home/disk/drizzle/wgk/microtia_first/marker_go.Rdata'
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
all.genes <- rownames(seurat.chon@assays$RNA@data)
use.genes <- intersect(keys(org.Hs.eg.db, keytype = "SYMBOL"), all.genes)
status <- c('Microtia_increase', 'Microtia_decrease')
list.go <- list()
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
                         universe = use.genes, pvalueCutoff = 0.05)
        res.egmt <- egmt@result
        list.go[[as.character(type)]] <- res.egmt
        res.kegg <- enrichKEGG(gene = df.symbol2geneid[genes.input, 'Gene_id'], 
                               universe = as.character(df.symbol2geneid[use.genes, 'Gene_id']), 
                               pvalueCutoff = 0.05, keyType = 'kegg')
        list.kegg[[as.character(type)]] <- res.kegg@result
    }
}

file.go <- '/home/disk/drizzle/wgk/microtia_first/GO.Rdata'
saveRDS(list.go, file = file.go)
file.kegg <- '/home/disk/drizzle/wgk/microtia_first/kegg.Rdata'
saveRDS(list.kegg, file = file.kegg)

df.geneid2symbol[unlist(strsplit((list.kegg$`Chondrocyte1_Microtia_decrease`)['hsa04978', 'geneID'], '/')),]
FeaturePlot(seurat.chon, features = c('MT1X'))


