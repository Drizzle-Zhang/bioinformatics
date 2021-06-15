.libPaths("/home/zy/tools/R-4.0.0/library")
setwd('/home/zy/my_git/bioinformatics/wgk')
library(Seurat)

file.seurat <- '/home/disk/drizzle/wgk/data/ear.seurat.first.Rdata'
seurat.first <- readRDS(file.seurat)
seurat.first <- subset(seurat.first, subset = sample %in% c('C1', 'C2', 'C3'))

DimPlot(seurat.first, group.by = "sample")
DimPlot(seurat.first, group.by = "RNA_snn_res.0.8", label = T)

table(seurat.first$sample, seurat.first$RNA_snn_res.0.8)/as.vector(table(seurat.first$sample))

clusters <- unique(seurat.first$RNA_snn_res.0.8)
list.marker <- list()
for (cluster in clusters) {
    sub.markers <- FindMarkers(seurat.first, ident.1 = cluster, group.by = 'RNA_snn_res.0.8',
                               logfc.threshold = 0.3, min.diff.pct = 0.1, only.pos = T)
    list.marker[[cluster]] <- sub.markers
}

file.marker.first <- '/home/disk/drizzle/wgk/data/seurat_marker_first_normal_0.8.Rdata'
saveRDS(list.marker, file = file.marker.first)
list.marker <- readRDS(file.marker.first)

# file.marker.first <- '/home/disk/drizzle/wgk/data/seurat_first_marker.Rdata'
# list.marker <- readRDS(file.marker.first)

# enrich GO
library(dplyr)
df.gene_id <- read.delim('/home/disk/drizzle/wgk/ncbi/Homo_sapiens.gene_info')
df.gene_id <- df.gene_id %>% distinct(Symbol, .keep_all = T)
df.symbol2geneid <- data.frame(Gene_id = df.gene_id$GeneID, row.names = df.gene_id$Symbol)
df.geneid2symbol <- data.frame(Symbol = df.gene_id$Symbol, row.names = df.gene_id$GeneID)
library(clusterProfiler)
library(org.Hs.eg.db)
list.go <- list()
list.kegg <- list()
clusters <- unique(seurat.first$RNA_snn_res.0.8)
for (cluster in clusters) {
    sub.markers <- list.marker[[cluster]]
    genes <- row.names(sub.markers)
    back.genes <- rownames(seurat.first@assays$RNA@counts)
    egmt <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                     universe = rownames(seurat.first@assays$RNA@counts), pvalueCutoff = 0.05)
    res.egmt <- egmt@result
    list.go[[cluster]] <- res.egmt
    res.kegg <- enrichKEGG(gene = df.symbol2geneid[genes, 'Gene_id'], 
                           universe = df.symbol2geneid[back.genes, 'Gene_id'], 
                           pvalueCutoff = 0.05, keyType = 'kegg')
    list.kegg[[cluster]] <- res.kegg@result
}

file.marker.first.go <- '/home/disk/drizzle/wgk/data/seurat_first_marker_go_0.8.Rdata'
saveRDS(list.go, file = file.marker.first.go)
list.go <- readRDS(file.marker.first.go)
file.marker.first.kegg <- '/home/disk/drizzle/wgk/data/seurat_first_marker_kegg_0.8.Rdata'
saveRDS(list.kegg, file = file.marker.first.kegg)
list.kegg <- readRDS(file.marker.first.kegg)


view.gene <- unlist(strsplit((list.kegg$`1`)['hsa05022', 'geneID'], '/'))
df.geneid2symbol[unlist(strsplit((list.kegg$`9`)['hsa04350', 'geneID'], '/')),]

# GSEA
library(dplyr)
df.gene_id <- read.delim('/home/disk/drizzle/wgk/ncbi/Homo_sapiens.gene_info')
df.gene_id <- df.gene_id %>% distinct(Symbol, .keep_all = T)
df.symbol2geneid <- data.frame(Gene_id = df.gene_id$GeneID, row.names = df.gene_id$Symbol)
df.geneid2symbol <- data.frame(Symbol = df.gene_id$Symbol, row.names = df.gene_id$GeneID)

list.marker.gsea <- list()
clusters <- unique(seurat.first$RNA_snn_res.0.8)
for (cluster in clusters) {
    sub.markers <- FindMarkers(seurat.first, ident.1 = cluster, group.by = 'RNA_snn_res.0.8',
                               logfc.threshold = 0, min.pct = 0)
    list.marker.gsea[[cluster]] <- sub.markers
}

file.marker.first.gsea <- '/home/disk/drizzle/wgk/data/seurat_first_marker_GSEA_0.8.Rdata'
saveRDS(list.marker.gsea, file = file.marker.first.gsea)
list.marker.gsea <- readRDS(file.marker.first.gsea)

library(clusterProfiler)
library(org.Hs.eg.db)
list.gsea <- list()
list.gsea.kegg <- list()
clusters <- unique(seurat.first$RNA_snn_res.0.8)
for (cluster in clusters) {
    sub.markers <- list.marker.gsea[[cluster]]
    geneList <- sub.markers[, 'avg_logFC']
    names(geneList) <- row.names(sub.markers)
    geneList <- geneList[order(geneList, decreasing = T)]
    egmt <- gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", pvalueCutoff = 0.7)
    res.egmt <- egmt@result
    list.gsea[[cluster]] <- res.egmt
    # KEGG
    names(geneList) <- df.symbol2geneid[row.names(sub.markers), 'Gene_id']
    res.kegg <- gseKEGG(geneList = geneList, 
                        pvalueCutoff = 0.7)
    list.gsea.kegg[[cluster]] <- res.kegg@result
}

file.gsea.first.go <- '/home/disk/drizzle/wgk/data/GSEA_first_go_0.8.Rdata'
saveRDS(list.go, file = file.gsea.first.go)
list.go <- readRDS(file.gsea.first.go)
file.gsea.first.kegg <- '/home/disk/drizzle/wgk/data/GSEA_first_kegg_0.8.Rdata'
saveRDS(list.kegg, file = file.gsea.first.kegg)
list.kegg <- readRDS(file.gsea.first.kegg)

view.gene <- df.geneid2symbol[unlist(strsplit((list.gsea.kegg$`2`)['hsa04380', 'core_enrichment'], '/')),]


View(list.gsea$`1`[list.gsea$`1`$enrichmentScore > 0,])

# feature plot 
# bone morphogenesis
FeaturePlot(seurat.first, features = c('SFRP2', 'COL1A1', 'COCH', 'SERPINH1'))
# 
FeaturePlot(seurat.first, features = c('HSPA1A', 'HSPA1B', 'CRYAB', 'HSP90AA1'))
# response to temperature stimulus
FeaturePlot(seurat.first, features = c('HSPA6', 'HSPA1A', 'HSPA1B', 'DNAJB1'))
# skeletal system development
FeaturePlot(seurat.first, features = c('SFRP2', 'COL1A1', 'COCH', 'CYTL1'))

FeaturePlot(seurat.first, features = c('GATA6'))

# regulation of stem cell differentiation
FeaturePlot(seurat.first, features = c('ANK3', 'HSP90AA1', 'CTNNA1', 'MYOT'))

FeaturePlot(seurat.first, 
            features = c('CFD', 'MMP10', 'COCH', 'BASP1', 'SFRP4', 'FGL2', 'FBLN1', 'DPT'), 
            ncol = 3)
FeaturePlot(seurat.first, 
            features = c('HSPA6', 'COL1A1', 'IER5L', 'GADD45G', 'IGFBP5', 'OGN', 
                         'ASPN', 'PI16', 'CRABP2', 'FBLN2', 'CLIC2', 'GPNMB'), 
            ncol = 3)
FeaturePlot(seurat.first, 
            features = c('FST', 'CTGF', 'FRZB', 'FGFBP2', 'PLA2G2A', 'S100A1', 
                         'MT1M', 'S100B', 'PRELP'), 
            ncol = 3)
FeaturePlot(seurat.first, 
            features = c('LAMB3', 'NOS2', 'LOX', 'ELN', 'COL2A1', 'CX3CL1', 
                         'COL11A2', 'SCUBE3', 'COLGALT2'), 
            ncol = 3)
FeaturePlot(seurat.first, 
            features = c('CYTL1', 'SNHG12', 'HSPB1', 'DLGAP1-AS2', 'NDRG2', 
                         'RASL11B', 'DUSP2', 'HSPA6', 'MDFI', 'GPRC5C'), 
            ncol = 3)

FeaturePlot(seurat.first, 
            features = c('ACTA2', 'ACTG2'), 
            ncol = 2)

##### extracellular matrix
# Fibulin
FeaturePlot(seurat.first, features = c('FBLN1', 'FBLN2', 'FBLN3', 'FBLN4', 'FBLN5', 'FBLN6', 'FBLN7'))
# Fibronectin

# Laminin
FeaturePlot(seurat.first, features = c('LAMA1', 'LAMA2', 'LAMA3', 'LAMA4', 'LAMA5', 
                                       'LAMB1', 'LAMB2', 'LAMB3', 'LAMB4',
                                       'LAMC1', 'LAMC2', 'LAMC3'))


# Thrombospondin
FeaturePlot(seurat.first, features = c('THBS1', 'THBS2', 'THBS3', 'THBS4'))


##### growth factor binding
# IGF
FeaturePlot(seurat.first, features = c('IGFBP1', 'IGFBP2', 'IGFBP3', 'IGFBP4', 'IGFBP5', 'IGFBP6', 'IGFBP7'))


##### genes in cluster
# 9
unlist(strsplit((list.kegg$`1`)['hsa05022', 'geneID'], '/'))

# findmarker <- function(seurat.first, cluster) {
#     sub.markers <- FindMarkers(seurat.first, ident.1 = cluster, group.by = 'vec.cluster',
#                                logfc.threshold = 0.3, min.diff.pct = 0.2)
#     return(sub.markers)
# }
# 
# library(foreach)
# library(doParallel)
# registerDoParallel(10)
# clusters <- unique(vec.cluster)
# list.marker <- foreach(cluster = clusters) %do% 
#     findmarker(seurat.first, cluster)

