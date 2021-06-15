library(Seurat)

path.data <- '/home/disk/drizzle/wgk/data/AllSample_2_merge/'
file.all <- paste0(path.data, 'seurat_celltype.Rdata')
seurat.all <- readRDS(file.all)

DimPlot(seurat.all, group.by = "batch")
DimPlot(seurat.all, group.by = "type")
DimPlot(seurat.all, group.by = "celltype", label = T)
table(seurat.all$batch, seurat.all$celltype)/as.vector(table(seurat.all$batch))

# chond
seurat.chon <- subset(seurat.all, subset = celltype %in% 
                          c('Chondrocyte1', 'Chondrocyte2', 'Chondral stem cell', 
                            'Transitional chondrocyte'))
DimPlot(seurat.chon, group.by = "celltype", label = T)
table(seurat.chon$batch, seurat.chon$celltype)/as.vector(table(seurat.chon$batch))


# chond & stromal
seurat.cs <- subset(seurat.all, subset = celltype %in% 
                          c('Chondrocyte1', 'Chondrocyte2', 'Chondral stem cell', 
                            'Transitional chondrocyte', 'Stromal stem cell',
                            'Stromal cell1', 'Stromal cell2'))
DimPlot(seurat.cs, group.by = "celltype", label = T)

# marker gene
fc.cutoff <- 0.5
pct.cutoof <- 0.1
celltypes <- unique(seurat.cs$celltype)
list.marker.go <- list()
for (cell in celltypes) {
    sub.seurat <- subset(seurat.cs, subset = celltype == cell)
    sub.markers <- FindMarkers(sub.seurat, ident.1 = 'Microtia', group.by = 'type',
                               logfc.threshold = fc.cutoff, min.diff.pct = pct.cutoof)
    list.marker.go[[as.character(cell)]] <- sub.markers
}
# FeaturePlot(sub.seurat, features = c('FRZB'))

file.marker.go <- paste0(path.M3, 'marker_go.Rdata')
saveRDS(list.marker.go, file = file.marker.go)
list.marker.go <- readRDS(file.marker.go)

# enrich GO
library(clusterProfiler)
library(org.Hs.eg.db)
all.genes <- rownames(seurat.cs@assays$RNA@data)
use.genes <- intersect(keys(org.Hs.eg.db, keytype = "SYMBOL"), all.genes)
status <- c('Microtia_increase', 'Microtia_decrease')
list.go <- list()
for (cell in celltypes) {
    for (st in status) {
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
        type <- paste0(as.character(cell), '_', st)
        list.go[[as.character(type)]] <- res.egmt
    }
}


file.go <- '/home/disk/drizzle/wgk/microtia_first/GO.Rdata'
saveRDS(list.go, file = file.go)
file.kegg <- '/home/disk/drizzle/wgk/microtia_first/kegg.Rdata'
saveRDS(list.kegg, file = file.kegg)

df.geneid2symbol[unlist(strsplit((list.kegg$`Chondrocyte1_Microtia_decrease`)['hsa04978', 'geneID'], '/')),]
FeaturePlot(seurat.chon, features = c('MT1X'))


