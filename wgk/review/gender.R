library(Seurat)
library(ggplot2)

file.chon <- '/mdshare/node9/zy/wgk/RDS/seurat_celltype.Rdata'
seurat.chon <- readRDS(file.chon)

VlnPlot(seurat.chon, features = 'SOX8', group.by = 'batch', split.by = 'celltype')
VlnPlot(seurat.chon, features = 'EGR1', group.by = 'batch', split.by = 'celltype')

# C2 vs C3
seurat.C2.C3 <- subset(seurat.chon, subset = batch %in% c('C2', 'C3'))
# marker gene
fc.cutoff <- 0.4
pct.cutoof <- 0
celltypes <- unique(seurat.C2.C3$celltype)
list.marker.gender <- list()
for (cell in celltypes) {
    sub.seurat <- subset(seurat.C2.C3, subset = celltype == cell)
    sub.markers <- FindMarkers(sub.seurat, ident.1 = 'C2', group.by = 'batch',
                               logfc.threshold = fc.cutoff, min.diff.pct = pct.cutoof)
    sub.markers <- sub.markers[sub.markers$p_val_adj < 0.05,]
    list.marker.gender[[as.character(cell)]] <- sub.markers
}

# microtia vs NC
seurat.child <- subset(seurat.chon, subset = batch %in% c('C4', 'C6', 'M1', 'M2', 'M3'))
list.marker.microtia <- list()
for (cell in celltypes) {
    sub.seurat <- subset(seurat.child, subset = celltype == cell)
    sub.markers <- FindMarkers(sub.seurat, ident.1 = 'Microtia', group.by = 'type',
                               logfc.threshold = fc.cutoff, min.diff.pct = pct.cutoof)
    sub.markers <- sub.markers[sub.markers$p_val_adj < 0.05,]
    list.marker.microtia[[as.character(cell)]] <- sub.markers
}

# test
df_test <- data.frame(stringsAsFactors = F)
for (cell in celltypes) {
    sub.seurat.C2.C3 <- subset(seurat.C2.C3, subset = celltype == cell)
    gene.gender <- rownames(sub.seurat.C2.C3@assays$RNA@counts)[rowSums(as.matrix(sub.seurat.C2.C3@assays$RNA@counts)) > 0]
    sub.seurat.child <- subset(seurat.child, subset = celltype == cell)
    gene.microtia <- rownames(sub.seurat.child@assays$RNA@counts)[rowSums(as.matrix(sub.seurat.child@assays$RNA@counts)) > 0]
    gene.background <- union(gene.gender, gene.microtia)
    marker.gender <- list.marker.gender[[cell]]
    merker.microtia <- list.marker.microtia[[cell]]
    # up
    marker.gender.up <- rownames(marker.gender[marker.gender$avg_logFC > 0,])
    merker.microtia.up <- rownames(merker.microtia[merker.microtia$avg_logFC > 0,])
    df_up <- data.frame(gender = (gene.background %in% marker.gender.up),
                        microtia = (gene.background %in% merker.microtia.up),
                        row.names = gene.background)
    table_up <- xtabs(~ gender + microtia, data = df_up)
    fisher.test(table_up)
}


path_C2M <- '/mdshare/node9/zy/wgk/C2M/'
# C2 vs M1M2M3
seurat.C2.M <- subset(seurat.chon, subset = batch %in% c('C2', 'M1', 'M2', 'M3'))
# marker gene
fc.cutoff <- 0
pct.cutoof <- 0
celltypes <- unique(seurat.C2.M$celltype)
list.marker.go.C2.M <- list()
for (cell in celltypes) {
    sub.seurat <- subset(seurat.C2.M, subset = celltype == cell)
    sub.markers <- FindMarkers(sub.seurat, ident.1 = 'Microtia', group.by = 'type',
                               logfc.threshold = fc.cutoff, min.diff.pct = pct.cutoof)
    list.marker.go.C2.M[[as.character(cell)]] <- sub.markers
}

list.all.diff <- list.marker.go
file.all.diff <- paste0(path_C2M, 'diff_genes.Rdata')
saveRDS(list.all.diff, file = file.all.diff)

# C3 vs M1M2M3
seurat.C3.M <- subset(seurat.chon, subset = batch %in% c('C3', 'M1', 'M2', 'M3'))



