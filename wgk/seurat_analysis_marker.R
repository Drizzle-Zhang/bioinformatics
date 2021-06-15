setwd('/home/zy/my_git/bioinformatics/wgk')
library(Seurat)

file.seurat <- '/home/disk/drizzle/wgk/data/ear.seurat.first.Rdata'
seurat.first <- readRDS(file.seurat)

cluster_0.8 <- seurat.first$RNA_snn_res.0.8
vec.cluster <- rep('_', length(cluster_0.8))
vec.cluster[cluster_0.8 %in% c(4, 5, 10)] <- '1'
vec.cluster[cluster_0.8 %in% c(9)] <- '2'
vec.cluster[cluster_0.8 %in% c(6)] <- '3'
vec.cluster[cluster_0.8 %in% c(0, 1)] <- '4'
vec.cluster[cluster_0.8 %in% c(2, 3)] <- '5'
vec.cluster[cluster_0.8 %in% c(14)] <- '6'
vec.cluster[cluster_0.8 %in% c(7, 11, 12)] <- '7'
vec.cluster[cluster_0.8 %in% c(8, 15)] <- '8'
vec.cluster[cluster_0.8 %in% c(13, 17)] <- '9'
vec.cluster[cluster_0.8 %in% c(16)] <- '10'
table(seurat.first$sample, vec.cluster)/as.vector(table(seurat.first$sample))

seurat.first$vec.cluster <- vec.cluster

clusters <- unique(seurat.first$RNA_snn_res.0.8)
list.marker <- list()
for (cluster in clusters) {
    sub.markers <- FindMarkers(seurat.first, ident.1 = cluster, group.by = 'RNA_snn_res.0.8',
                               logfc.threshold = 0.3, min.diff.pct = 0.1, only.pos = T)
    list.marker[[cluster]] <- sub.markers
}

file.marker.first <- '/home/disk/drizzle/wgk/data/seurat_first_marker.Rdata'
saveRDS(list.marker, file = file.marker.first)
list.marker <- readRDS(file.marker.first)


file.marker.first <- '/home/disk/drizzle/wgk/data/seurat_first_marker.Rdata'
list.marker <- readRDS(file.marker.first)



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
clusters <- unique(vec.cluster)
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

file.marker.first.go <- '/home/disk/drizzle/wgk/data/seurat_first_marker_go.Rdata'
saveRDS(list.go, file = file.marker.first.go)
list.go <- readRDS(file.marker.first.go)
file.marker.first.kegg <- '/home/disk/drizzle/wgk/data/seurat_first_marker_kegg.Rdata'
saveRDS(list.kegg, file = file.marker.first.kegg)
list.kegg <- readRDS(file.marker.first.kegg)


view.gene <- unlist(strsplit((list.kegg$`1`)['hsa05022', 'geneID'], '/'))


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

file.marker.first.gsea <- '/home/disk/drizzle/wgk/data/seurat_first_marker_GSEA.Rdata'
saveRDS(list.marker.gsea, file = file.marker.first.gsea)
list.marker.gsea <- readRDS(file.marker.first.gsea)

library(clusterProfiler)
library(org.Hs.eg.db)
list.gsea <- list()
list.gsea.kegg <- list()
clusters <- unique(seurat.first$vec.cluster)
for (cluster in clusters) {
    sub.markers <- list.marker.gsea[[cluster]]
    geneList <- sub.markers[, 'avg_logFC']
    names(geneList) <- row.names(sub.markers)
    geneList <- geneList[order(geneList, decreasing = T)]
    # egmt <- gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", pvalueCutoff = 0.7)
    # res.egmt <- egmt@result
    # list.gsea[[cluster]] <- res.egmt
    # KEGG
    names(geneList) <- df.symbol2geneid[row.names(sub.markers), 'Gene_id']
    res.kegg <- gseKEGG(geneList = geneList, 
                        pvalueCutoff = 0.7)
    list.gsea.kegg[[cluster]] <- res.kegg@result
}

file.res.gsea <- '/home/disk/drizzle/wgk/data/GSEA_first.Rdata'
saveRDS(list.gsea, file = file.res.gsea)
list.gsea <- readRDS(file.res.gsea)

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

