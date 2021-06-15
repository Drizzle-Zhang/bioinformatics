library(Seurat)

file.seurat <- '/home/disk/drizzle/wgk/data/ear.seurat.first.Rdata'
seurat.first <- readRDS(file.seurat)

file.bone <- '/home/disk/drizzle/wgk/bone_genes.txt'
df.genes.bone <- read.delim(file.bone, row.names = 1)
bone.genes <- intersect(rownames(seurat.first@assays$RNA@scale.data), rownames(df.genes.bone))
clusters <- as.character(seurat.first$RNA_snn_res.0.8)
sel.genes <- c()
for (cluster in as.character(c(4, 5, 10, 9, 6, 0, 1, 2, 3, 14))) {
    sel.genes <- c(sel.genes, intersect(bone.genes, rownames(list.marker[[cluster]])))
}
sel.genes <- unique(sel.genes)
# DoHeatmap(seurat.first, features = unique(sel.genes), group.by = 'RNA_snn_res.0.8')


clusters <- as.character(seurat.first$RNA_snn_res.0.8)
mat_exp <- seurat.first@assays$RNA@scale.data[sel.genes, clusters == '4']
annotation_c <- data.frame(cluster = rep('4', ncol(mat_exp)))
for (cluster in as.character(c(5, 10, 9, 6, 0, 1, 2, 3, 14))) {
    mat_exp <- cbind(mat_exp, seurat.first@assays$RNA@scale.data[sel.genes, clusters == cluster])
    annotation_c <- rbind(annotation_c, 
                          data.frame(cluster = rep(cluster, ncol(seurat.first@assays$RNA@scale.data[sel.genes, clusters == cluster]))))
}
row.names(annotation_c) <- colnames(mat_exp)
# pheatmap::pheatmap(mat_exp, 
#                    cluster_rows = T, cluster_cols = F, scale = "row",
#                    annotation_col = annotation_c,
#                    annotation_legend = T,
#                    show_rownames = T, show_colnames = F,
#                    color = colorRampPalette(c("#FF00FF", "#000000","#FFFF00"))(100)
#             )

seurat.plot <- seurat.first
seurat.plot@assays$RNA@scale.data <- mat_exp
seurat.plot@meta.data <- seurat.plot@meta.data[colnames(mat_exp),]
seurat.plot$heatmap <- factor(annotation_c$cluster, as.character(c(4, 5, 10, 9, 6, 0, 1, 2, 3, 14)))
DoHeatmap(seurat.plot, features = sel.genes, cells = colnames(mat_exp), group.by = 'heatmap')

