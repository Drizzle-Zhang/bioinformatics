setwd('/home/zy/my_git/bioinformatics/wgk')
.libPaths('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(igraph)
library(ggraph)
library(tidygraph)

path.TF.net <- '/home/disk/drizzle/wgk/TF_net/'

path.scenic <- '/home/yzj/JingMA_NEW/res/SCENIC_All_2/'
path.env <- paste0(path.scenic, 'int/')

setwd(path.scenic)
RegulonInfo <- readRDS('int/2.5_regulonTargetsInfo.Rds')
Regulons <- readRDS('int/3.1_regulons_forAUCell.Rds')

RegulonInfo_HighConf <- RegulonInfo[RegulonInfo$highConfAnnot == TRUE,]

mat_weight <- reshape2::dcast(RegulonInfo_HighConf, TF ~ gene, value.var = 'Genie3Weight')
rownames(mat_weight) <- mat_weight$TF
mat_weight$TF <- NULL
mat_weight <- mat_weight[colnames(mat_weight), colnames(mat_weight)]
rownames(mat_weight) <- colnames(mat_weight)
mat_weight[is.na(mat_weight)] <- 0

igraph <- 
    graph_from_adjacency_matrix(as.matrix(mat_weight), 
                                mode = 'directed', 
                                weighted = NULL, diag = T)

vcount(igraph)


# use diff genes
path.M123 <- '/home/disk/drizzle/wgk/microtia_chon_child_M1M2M3/'
file.gsea.marker <- paste0(path.M123, 'marker_GSEA.Rdata')
list.marker.gsea <- readRDS(file.gsea.marker)

#######################################
# CSC
marker.CSC <- list.marker.gsea$CSC
down.marker.CSC <- marker.CSC[marker.CSC$p_val_adj < 0.01 & marker.CSC$avg_logFC < -0.1,]
down.gene.CSC <- rownames(down.marker.CSC)
RegulonInfo_CSC_down <- 
    RegulonInfo_HighConf[RegulonInfo_HighConf$TF %in% down.gene.CSC & 
                             RegulonInfo_HighConf$gene %in% down.gene.CSC,]

# down
mat_weight_down <- reshape2::dcast(RegulonInfo_CSC_down, 
                                   TF ~ gene, value.var = 'Genie3Weight')
rownames(mat_weight_down) <- mat_weight_down$TF
mat_weight_down$TF <- NULL
row_add <- setdiff(colnames(mat_weight_down), rownames(mat_weight_down))
df.row_add <- data.frame(matrix(rep(NA, length(row_add)*ncol(mat_weight_down)), 
                                nrow = length(row_add), ncol = ncol(mat_weight_down)),
                         row.names = row_add)
colnames(df.row_add) <- colnames(mat_weight_down)
mat_weight_down <- rbind(mat_weight_down, df.row_add)
col_add <- setdiff(rownames(mat_weight_down), colnames(mat_weight_down))
df.col_add <- data.frame(matrix(rep(NA, length(col_add)*nrow(mat_weight_down)), 
                                nrow = nrow(mat_weight_down), ncol = length(col_add)),
                         row.names = rownames(mat_weight_down))
colnames(df.col_add) <- col_add
mat_weight_down <- cbind(mat_weight_down, df.col_add)
mat_weight_down[is.na(mat_weight_down)] <- 0
TFgene <- sort(colnames(mat_weight_down))
mat_weight_down <- mat_weight_down[TFgene, TFgene]

igraph_down <- 
    graph_from_adjacency_matrix(t(as.matrix(mat_weight_down)), 
                                mode = 'directed', 
                                weighted = T, diag = T)
V(igraph_down)$degree <- degree(igraph_down, normalized = T)
V(igraph_down)$weight_degree <- strength(igraph_down)
V(igraph_down)$closeness_centrality <- closeness(igraph_down, normalized = T)
V(igraph_down)$betweenness_centrality <- betweenness(igraph_down, normalized = T)
V(igraph_down)$eigenvector_centrality <- evcent(igraph_down)$vector
V(igraph_down)$page_rank <- page_rank(igraph_down)$vector
V(igraph_down)$betweenness_centrality_directed <- betweenness(igraph_down, directed = T, normalized = T)
V(igraph_down)$eigenvector_centrality_directed <- evcent(igraph_down, directed = T)$vector

node_list_down <- data.frame(
    node_id = V(igraph_down)$name,
    degree = V(igraph_down)$degree,
    weight_degree = V(igraph_down)$weight_degree,
    closeness_centrality = V(igraph_down)$closeness_centrality,
    betweenness_centrality = V(igraph_down)$betweenness_centrality,
    eigenvector_centrality = V(igraph_down)$eigenvector_centrality,
    page_rank = V(igraph_down)$page_rank,
    betweenness_centrality_directed = V(igraph_down)$betweenness_centrality_directed,
    eigenvector_centrality_directed = V(igraph_down)$eigenvector_centrality_directed)

igraph_down_2 <- 
    graph_from_adjacency_matrix((as.matrix(mat_weight_down)), 
                                mode = 'directed', 
                                weighted = T, diag = T)
V(igraph_down_2)$page_rank <- page_rank(igraph_down)$vector
FCs <- abs(down.marker.CSC[V(igraph_down)$name, 'avg_logFC'])
norm_FCs <- FCs
norm_FCs[FCs > 0.6] <- 0.6
# norm_FCs <- (FCs - min(FCs))/(max(FCs) - min(FCs))
# names(norm_FCs) <- V(igraph_down)$name
V(igraph_down_2)$FC <- norm_FCs

# node group
node_type <- V(igraph_down)$name
node_group <- rep('3', length(node_type))
node_group[node_type %in% c('SOX5', 'SOX8', 'DBP', 'ARID5B')] <- '1'
node_group[node_type %in% c('A2M', 'ACAN', 'CHAD', 'COL11A1', 'COL2A1', 'COL9A3',
                            "CTGF", "ELN", "FMOD", "MT1G", "MT2A", "S100A1", 
                            "S100B", "SCRG1", "SPARC", "VIT", "WWP2", "SOD3", 'TXN')] <- '2'
V(igraph_down_2)$group <- node_group

ggraph_CSC_down <- as_tbl_graph(igraph_down_2)

# plot_CSC_down <- 
#     # ggraph(ggraph_CSC_down, layout = 'centrality',cent=page_rank) + 
#     ggraph(ggraph_CSC_down, layout = 'stress') + 
#     geom_edge_link(aes(edge_width=weight),color="gray", alpha = 0.7,
#                    arrow = arrow(length = unit(2, 'mm')), 
#                    end_cap = circle(1.5, 'mm'), linejoin = 'bevel', 
#                    start_cap = circle(3, 'mm')) +
#     scale_edge_width(range=c(0.3,1)) + 
#     scale_size_continuous(range = c(2,15)) + 
#     geom_node_point(aes(size = page_rank, fill = group, color = group, alpha = group),
#                     shape=21) + 
#     scale_color_manual(values = c('#B6CEF2', '#4B6BF2', 'dimgray')) + 
#     scale_fill_manual(values = c('#B6CEF2', '#4B6BF2', 'dimgray')) + 
#     # scale_alpha_manual(values = c(0.8, 1, 0.6)) + 
#     geom_node_text(aes(filter = group == 1,label=name),size=5) + 
#     geom_node_text(aes(filter = group == 2,label=name),size=4, repel = T) + 
#     theme_void() + theme(legend.position = 'none')
# ggsave(plot = plot_CSC_down, path = path.TF.net, 
#        filename = 'TF_CSC_down.pdf',
#        height = 16, width = 16, units = 'cm')

plot_CSC_down <-
    # ggraph(ggraph_CSC_down, layout = 'centrality',cent=page_rank) +
    ggraph(ggraph_CSC_down, layout = 'stress') +
    geom_edge_link(aes(edge_width=weight, alpha = weight),color="gray",
                   arrow = arrow(length = unit(2, 'mm')),
                   end_cap = circle(1.5, 'mm'), 
                   start_cap = circle(3, 'mm')) +
    scale_edge_width(range=c(0.5,1)) +
    scale_edge_alpha(range=c(0.2,1)) + 
    scale_size_continuous(range = c(2,15)) +
    geom_node_point(aes(size = page_rank, fill = group, alpha = FC),
                    shape=21, color = 'transparent') +
    # scale_color_manual(values = c('#0000CD', '#DC143C', 'dimgray')) +
    scale_fill_manual(values = c('#4169E1', '#FF4500', 'dimgray')) +
    # scale_alpha_manual(values = c(0.8, 1, 0.6)) +
    geom_node_text(aes(filter = group == 1,label=name),size=5) +
    geom_node_text(aes(filter = group == 2,label=name),size=4, repel = T) +
    theme_void() + theme(legend.position = 'none')
ggsave(plot = plot_CSC_down, path = path.TF.net,
       filename = 'TF_CSC_down.png',
       height = 16, width = 16, units = 'cm')


# up
up.marker.CSC <- marker.CSC[marker.CSC$p_val_adj < 0.01 & marker.CSC$avg_logFC > 0.4,]
up.gene.CSC <- rownames(up.marker.CSC)
RegulonInfo_CSC_up <- 
    RegulonInfo_HighConf[RegulonInfo_HighConf$TF %in% up.gene.CSC & 
                             RegulonInfo_HighConf$gene %in% up.gene.CSC,]

mat_weight_up <- reshape2::dcast(RegulonInfo_CSC_up, 
                                   TF ~ gene, value.var = 'Genie3Weight')
rownames(mat_weight_up) <- mat_weight_up$TF
mat_weight_up$TF <- NULL
row_add <- setdiff(colnames(mat_weight_up), rownames(mat_weight_up))
df.row_add <- data.frame(matrix(rep(NA, length(row_add)*ncol(mat_weight_up)), 
                                nrow = length(row_add), ncol = ncol(mat_weight_up)),
                         row.names = row_add)
colnames(df.row_add) <- colnames(mat_weight_up)
mat_weight_up <- rbind(mat_weight_up, df.row_add)
col_add <- setdiff(rownames(mat_weight_up), colnames(mat_weight_up))
df.col_add <- data.frame(matrix(rep(NA, length(col_add)*nrow(mat_weight_up)), 
                                nrow = nrow(mat_weight_up), ncol = length(col_add)),
                         row.names = rownames(mat_weight_up))
colnames(df.col_add) <- col_add
mat_weight_up <- cbind(mat_weight_up, df.col_add)
mat_weight_up[is.na(mat_weight_up)] <- 0
TFgene <- sort(colnames(mat_weight_up))
mat_weight_up <- mat_weight_up[TFgene, TFgene]

igraph_up <- 
    graph_from_adjacency_matrix(t(as.matrix(mat_weight_up)), 
                                mode = 'directed', 
                                weighted = T, diag = T)
V(igraph_up)$degree <- degree(igraph_up, normalized = T)
V(igraph_up)$weight_degree <- strength(igraph_up)
V(igraph_up)$closeness_centrality <- closeness(igraph_up, normalized = T)
V(igraph_up)$betweenness_centrality <- betweenness(igraph_up, normalized = T)
V(igraph_up)$eigenvector_centrality <- evcent(igraph_up)$vector
V(igraph_up)$page_rank <- page_rank(igraph_up)$vector
V(igraph_up)$betweenness_centrality_directed <- betweenness(igraph_up, directed = T, normalized = T)
V(igraph_up)$eigenvector_centrality_directed <- evcent(igraph_up, directed = T)$vector

node_list_up <- data.frame(
    node_id = V(igraph_up)$name,
    degree = V(igraph_up)$degree,
    weight_degree = V(igraph_up)$weight_degree,
    closeness_centrality = V(igraph_up)$closeness_centrality,
    betweenness_centrality = V(igraph_up)$betweenness_centrality,
    eigenvector_centrality = V(igraph_up)$eigenvector_centrality,
    page_rank = V(igraph_up)$page_rank,
    betweenness_centrality_directed = V(igraph_up)$betweenness_centrality_directed,
    eigenvector_centrality_directed = V(igraph_up)$eigenvector_centrality_directed)

igraph_up_2 <- 
    graph_from_adjacency_matrix((as.matrix(mat_weight_up)), 
                                mode = 'directed', 
                                weighted = T, diag = T)
V(igraph_up_2)$page_rank <- page_rank(igraph_up)$vector
FCs <- abs(up.marker.CSC[V(igraph_up)$name, 'avg_logFC'])
norm_FCs <- FCs
norm_FCs[FCs > 1] <- 1
# norm_FCs <- (FCs - min(FCs))/(max(FCs) - min(FCs))
# names(norm_FCs) <- V(igraph_down)$name
V(igraph_up_2)$FC <- norm_FCs

# node group
node_type <- V(igraph_up)$name
node_group <- rep('3', length(node_type))
node_group[node_type %in% c('EGR1', 'KLF10', 'JUNB', 'REL', 'EGR3',
                            'ATF3', 'HIVEP3', 'IRX2', 'EGR2', 'KLF2', 'BCL3',
                            'CEBPB', 'CEBPD', 'STAT5A')] <- '1'
node_group[node_type %in% c('MMP3', 'IL8', 'CXCL2', 'CRYAB', 'CXCL3',
                            "KDM6B", "CXCL1", "ICAM1", "NFKBIZ", 'NFKBIA', "ZC3H12A", "SOD2", 
                            "CLK1", "MAP2K3", "BMP2", "MMP2", "PMAIP1", 'CSF1',
                            'CDKN1A', 'DUSP1', 'TNIP1', 'PPP1R15A', 'NR4A2', 
                            # 'SKIL', 
                            'APOD', 'HSPA1A', 'NR1D1', 'ZFP36L1', 
                            'HSPB1', 'MCL1', 'NR4A3', 'ZFHX3', 'TNFAIP2',
                            'TNFAIP3', 'CYLD', '')] <- '2'
V(igraph_up_2)$group <- node_group

ggraph_CSC_up <- as_tbl_graph(igraph_up_2)
# plot_CSC_up <- 
#     # ggraph(ggraph_CSC_up, layout = 'linear', circular = TRUE) + 
#     # geom_edge_arc(aes(edge_width=weight),color="gray", alpha = 0.7,
#     #                              arrow = arrow(length = unit(2, 'mm')),
#     #                              end_cap = circle(1.5, 'mm'),
#     #                              start_cap = circle(3, 'mm')) +
#     # ggraph(ggraph_CSC_up, layout = 'centrality',cent=page_rank) +
#     ggraph(ggraph_CSC_up, layout = 'stress') +
#     geom_edge_link(aes(edge_width=weight),color="gray", alpha = 0.7,
#                    arrow = arrow(length = unit(2, 'mm')),
#                    end_cap = circle(1.5, 'mm'),
#                    start_cap = circle(3, 'mm')) +
#     scale_edge_width(range=c(0.3,1)) + 
#     scale_size_continuous(range = c(2,15)) + 
#     geom_node_point(aes(size = page_rank, fill = group, color = group, alpha = group),
#                     shape=21) + 
#     scale_color_manual(values = c('#6495ED', '#DC143C', 'dimgray')) + 
#     scale_fill_manual(values = c('#6495ED', '#DC143C', 'dimgray')) + 
#     scale_alpha_manual(values = c(0.8, 1, 0.6)) + 
#     geom_node_text(aes(filter = group == 1,label=name),size=4) + 
#     geom_node_text(aes(filter = group == 2,label=name),size=3, repel = T, 
#                    nudge_x = 0.08, nudge_y = 0.005) + 
#     theme_void() + theme(legend.position = 'none')
plot_CSC_up <- 
    ggraph(ggraph_CSC_up, layout = 'stress') +
    geom_edge_link(aes(edge_width=weight, alpha = weight),color="gray",
                   arrow = arrow(length = unit(2, 'mm')),
                   end_cap = circle(1.5, 'mm'), 
                   start_cap = circle(3, 'mm')) +
    scale_edge_width(range=c(0.5,1)) +
    scale_edge_alpha(range=c(0.2,1)) + 
    scale_size_continuous(range = c(2,15)) +
    geom_node_point(aes(size = page_rank, fill = group, alpha = FC),
                    shape=21, color = 'transparent') +
    # scale_color_manual(values = c('#0000CD', '#DC143C', 'dimgray')) +
    scale_fill_manual(values = c('#4169E1', '#FF4500', 'dimgray')) +
    # scale_alpha_manual(values = c(0.8, 1, 0.6)) + 
    geom_node_text(aes(filter = group == 1,label=name),size=4) + 
    geom_node_text(aes(filter = group == 2,label=name), 
                   # nudge_x = 0.08, nudge_y = 0.005,
                   size=3, repel = T) + 
    theme_void() + theme(legend.position = 'none')
ggsave(plot = plot_CSC_up, path = path.TF.net, 
       filename = 'TF_CSC_up.png',
       height = 16, width = 20, units = 'cm')

#################################################
# TC
marker.TC <- list.marker.gsea$TC
down.marker.TC <- marker.TC[marker.TC$p_val_adj < 0.01 & marker.TC$avg_logFC < -0.1,]
down.gene.TC <- rownames(down.marker.TC)
RegulonInfo_TC_down <- 
    RegulonInfo_HighConf[RegulonInfo_HighConf$TF %in% down.gene.TC & 
                             RegulonInfo_HighConf$gene %in% down.gene.TC,]

# down
mat_weight_down <- reshape2::dcast(RegulonInfo_TC_down, 
                                   TF ~ gene, value.var = 'Genie3Weight')
rownames(mat_weight_down) <- mat_weight_down$TF
mat_weight_down$TF <- NULL
TFgene <- unique(c(rownames(mat_weight_down), colnames(mat_weight_down)))
mat_weight_down <- as.matrix(mat_weight_down)[TFgene, TFgene]
rownames(mat_weight_down) <- TFgene
colnames(mat_weight_down) <- TFgene
mat_weight_down[is.na(mat_weight_down)] <- 0

igraph_down <- 
    graph_from_adjacency_matrix(t(as.matrix(mat_weight_down)), 
                                mode = 'directed', 
                                weighted = T, diag = T)
V(igraph_down)$degree <- degree(igraph_down, normalized = T)
V(igraph_down)$weight_degree <- strength(igraph_down)
V(igraph_down)$closeness_centrality <- closeness(igraph_down, normalized = T)
V(igraph_down)$betweenness_centrality <- betweenness(igraph_down, normalized = T)
V(igraph_down)$eigenvector_centrality <- evcent(igraph_down)$vector
V(igraph_down)$page_rank <- page_rank(igraph_down)$vector
V(igraph_down)$betweenness_centrality_directed <- betweenness(igraph_down, directed = T, normalized = T)
V(igraph_down)$eigenvector_centrality_directed <- evcent(igraph_down, directed = T)$vector

node_list_down <- data.frame(
    node_id = V(igraph_down)$name,
    degree = V(igraph_down)$degree,
    weight_degree = V(igraph_down)$weight_degree,
    closeness_centrality = V(igraph_down)$closeness_centrality,
    betweenness_centrality = V(igraph_down)$betweenness_centrality,
    eigenvector_centrality = V(igraph_down)$eigenvector_centrality,
    page_rank = V(igraph_down)$page_rank,
    betweenness_centrality_directed = V(igraph_down)$betweenness_centrality_directed,
    eigenvector_centrality_directed = V(igraph_down)$eigenvector_centrality_directed)

