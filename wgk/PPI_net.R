
path.out <- '/home/disk/drizzle/wgk/microtia_chon_child_M1M2M3/PPI_net/'

path.scenic <- '/home/yzj/JingMA_NEW/res/SCENIC_main/'
RegulonInfo <- readRDS(paste0(path.scenic, 'int/2.5_regulonTargetsInfo.Rds'))
RegulonInfo_HighConf <- RegulonInfo[RegulonInfo$highConfAnnot == TRUE,]
genes.scenic <- unique(RegulonInfo_HighConf$gene)

path.M123 <- '/home/disk/drizzle/wgk/microtia_chon_child_M1M2M3/'
file.gsea.marker <- paste0(path.M123, 'marker_GSEA.Rdata')
list.marker.gsea <- readRDS(file.gsea.marker)

cell <- 'CSC'
sub.marker <- list.marker.gsea[[cell]]
sub.marker.1 <- sub.marker[sub.marker$p_val_adj < 0.01 & 
                               sub.marker$avg_logFC < -0.1,]
sub.marker.2 <- sub.marker[sub.marker$p_val_adj < 0.01 & 
                               sub.marker$avg_logFC > 0.4,]
CSC.diff.genes <- intersect(c(rownames(sub.marker.1), rownames(sub.marker.2)), genes.scenic)

df.CSC <- data.frame(CSC.diff.genes)
write.table(df.CSC, file = paste0(path.out, 'diff_gene_CSC.txt'), 
            row.names = F, col.names = F, quote = F)
