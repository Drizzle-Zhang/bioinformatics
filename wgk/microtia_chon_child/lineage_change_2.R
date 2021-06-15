setwd('/home/zy/my_git/bioinformatics/wgk')
.libPaths('/home/zy/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)
.libPaths('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0')
library(ggplot2)
library(SCENIC)
require("RColorBrewer")
library(maSigPro)


path.data <- '/home/disk/drizzle/wgk/data/AllSample_2_merge/'
path.lineage <- paste0(path.data, 'chon_lineage/')
file.chon <- paste0(path.lineage, 'seurat_celltype.Rdata')
seurat.chon <- readRDS(file.chon)
seurat.child <- subset(seurat.chon, subset = batch %in% c('C4', 'C6', 'M1', 'M2', 'M3'))
# seurat.child <- NormalizeData(seurat.child)
# seurat.child <- FindVariableFeatures(seurat.child, nfeatures = 5000)
# highvar.genes <- VariableFeatures(seurat.child)
# seurat.child <- ScaleData(seurat.child, split.by = "batch", 
#                           features = rownames(seurat.child@assays$RNA@counts))

path.M123 <- '/home/disk/drizzle/wgk/microtia_chon_child_M1M2M3/'
path.change <- paste0(path.M123, 'chon_lineage_1/')
if (!file.exists(path.change)) {
    dir.create(path.change)
}


seurat.child.N <- subset(seurat.child, subset = type == 'Normal')
seurat.child.M <- subset(seurat.child, subset = type == 'Microtia')

library(foreach)
library(doParallel)
registerDoParallel(cores = 10)
corr.PC.exp <- function(mat.gene, PC2, gene) {
    vec.exp <- mat.gene[gene,]
    corr <- cor(PC2, vec.exp, method = 'spearman')
    return(data.frame(gene = gene, corr = corr))
}

# Normal
Harmony2_N <- seurat.child.N@reductions$harmony@cell.embeddings[, 'harmony_2']
mat.gene_N <- seurat.child.N@assays$RNA@data
df.corr_N <- 
    foreach(gene = highvar.genes, .combine = rbind) %dopar% 
    corr.PC.exp(mat.gene_N, Harmony2_N, gene)
rownames(df.corr_N) <- highvar.genes

# Microtia
Harmony2_M <- seurat.child.M@reductions$harmony@cell.embeddings[, 'harmony_2']
mat.gene_M <- seurat.child.M@assays$RNA@data
df.corr_M <- 
    foreach(gene = highvar.genes, .combine = rbind) %dopar% 
    corr.PC.exp(mat.gene_M, Harmony2_M, gene)
rownames(df.corr_M) <- highvar.genes

df.corr <- merge(df.corr_N, df.corr_M, by = 'gene')
df.corr[is.na(df.corr)] <- 0
df.corr$diff_M_N <- df.corr$corr.y - df.corr$corr.x


# plot single gene
Harmony2 <- seurat.child@reductions$harmony@cell.embeddings[, 'harmony_2']
mat.gene <- seurat.child@assays$RNA@data
# AUC
regulonAUC <- readRDS(file='/home/yzj/JingMA_NEW/res/SCENIC_main/int/3.4_regulonAUC.Rds')
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), colnames(mat.gene)]
mat.auc <- as.matrix(regulonAUC@assays@data@listData$AUC)
# mat.gene.TF <- rbind(as.matrix(mat.gene), mat.auc)
df.pc.gene <- data.frame(t(rbind(as.matrix(mat.gene), mat.auc)), check.names = F)
df.pc.gene$Harmony2 <- Harmony2
df.pc.gene$celltype <- seurat.child$celltype
df.pc.gene$status <- seurat.child$type
df.pc.gene <- df.pc.gene[order(Harmony2, decreasing = F),]
df.pc.gene$idx <- 1:nrow(df.pc.gene)
colors <- c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
names(colors) <- c('CSC', 'TC', 'C1', 'C2')

gene <- 'DDIT3'
# gene <- 'SOX5 (218g)'
df.plot <- df.pc.gene[, c('idx', 'Harmony2', 'celltype', 'status', gene)]
names(df.plot) <- c('idx', 'Harmony2', 'celltype', 'status', 'gene')
p.gene <-
    ggplot(data = df.plot, aes(x = idx, 
                               linetype = status, 
                               y = gene)) + 
    geom_point(aes(color = celltype), size = 0.3) + 
    scale_color_manual(labels = c('CSC', 'TC', 'C1', 'C2'),
                       values = colors) + 
    # xlim(-30, 10) + 
    geom_smooth(color = '#696969') +
    # geom_smooth(color = '#696969', formula = y~poly(x, 3), method = lm) + 
    labs(x = '', y = '') + 
    theme(panel.background=element_rect(fill='transparent', color='black',size = 1),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'none') +
    theme(panel.background=element_rect(fill='transparent', color='black',size = 1)) +
    annotate('text', label = gene, x = 22000, y = max(df.plot$gene), 
             hjust = 1, vjust = 1, size = 7)

# TF change
# gene <- 'SOX5 (218g)'
vec.TF <- c('SOX5 (218g)', 'DBP (45g)', 'SOX8 (158g)',
            'EGR1 (264g)', 'EGR3 (53g)', 'KLF10 (16g)',
            'ATF3 (83g)', 'REL (834g)', 'BCL3 (162g)',
            'JUNB (158g)', 'CEBPB (258g)', 'CEBPD (199g)')
vec.TF <- c('SOX5 (218g)', 'DBP (45g)', 'SOX8 (158g)',
            'EGR1 (264g)', 'EGR3 (53g)', 'KLF10 (16g)',
            'ATF3 (83g)')
for (i in 1:length(vec.TF)) {
    TF <- vec.TF[i]
    df.plot <- df.pc.gene[, c('idx', 'Harmony2', 'celltype', 'status', TF)]
    names(df.plot) <- c('idx', 'Harmony2', 'celltype', 'status', 'TF')
    p1 <- ggplot(data = df.plot, aes(x = Harmony2, 
                                     linetype = status, 
                                     y = TF)) + 
        geom_point(aes(color = celltype), size = 0.1) + 
        scale_color_manual(labels = c('CSC', 'TC', 'C1', 'C2'),
                           values = colors) + 
        xlim(-30, 10) +
        geom_smooth(color = '#696969') +
        # geom_smooth(color = '#696969', formula = y~poly(x, 3), method = lm) + 
        labs(x = '', y = '') + 
        theme(panel.background=element_rect(fill='transparent', color='black',size = 1),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.position = 'none') +
        theme(panel.background=element_rect(fill='transparent', color='black',size = 1)) +
        annotate('text', label = TF, x = 9, y = max(df.plot$TF), 
                 hjust = 1, vjust = 1, size = 7)
    if (i == 1) {
        p <- p1
    } else {
        p <- p / p1
    }
}
ggsave(plot = p, path = path.change, 
       filename = 'TF_AUC.png',
       height = 35, width = 10, units = 'cm')

vec.TF.exp <- c('SOX5', 'DBP', 'SOX8',
            'EGR1', 'EGR3', 'KLF10',
            'ATF3')
for (i in 1:length(vec.TF.exp)) {
    TF <- vec.TF.exp[i]
    df.plot <- df.pc.gene[, c('idx', 'Harmony2', 'celltype', 'status', TF)]
    names(df.plot) <- c('idx', 'Harmony2', 'celltype', 'status', 'TF')
    p1 <- ggplot(data = df.plot, aes(x = Harmony2, 
                                     linetype = status, 
                                     y = TF)) + 
        geom_point(aes(color = celltype), size = 0.1) + 
        scale_color_manual(labels = c('CSC', 'TC', 'C1', 'C2'),
                           values = colors) + 
        xlim(-30, 10) +
        geom_smooth(color = '#696969') +
        # geom_smooth(color = '#696969', formula = y~poly(x, 3), method = lm) + 
        labs(x = '', y = '') + 
        theme(panel.background=element_rect(fill='transparent', color='black',size = 1),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.position = 'none') +
        theme(panel.background=element_rect(fill='transparent', color='black',size = 1)) +
        annotate('text', label = TF, x = 9, y = max(df.plot$TF), 
                 hjust = 1, vjust = 1, size = 7)
    if (i == 1) {
        p <- p1
    } else {
        p <- p / p1
    }
}
ggsave(plot = p, path = path.change, 
       filename = 'TF_EXP.png',
       height = 35, width = 10, units = 'cm')

    
# time series
# four time
library(maSigPro)
mat.count <- seurat.child@assays$RNA@counts[highvar.genes,]
Time <- as.character(seurat.child$celltype)
Time[Time == 'CSC'] <- 1
Time[Time == 'TC'] <- 2
Time[Time == 'C1'] <- 3
Time[Time == 'C2'] <- 4

# 20 points
mat.count <- seurat.child@assays$RNA@counts
Time <- as.character(seurat.child$celltype)
Harmony2 <- seurat.child@reductions$harmony@cell.embeddings[, 'harmony_2']
quantile.H2 <- quantile(Harmony2, probs = seq(0, 1, 0.05))
for (i in 1:20) {
    Time[Harmony2 >= quantile.H2[i] & Harmony2 <= quantile.H2[i+1]] <- i
}

# time series analysis
Time <- as.numeric(Time)
Replicates <- as.numeric(as.factor(seurat.child$batch))
Normal <- rep(0, length(seurat.child$type))
Normal[seurat.child$type == 'Normal'] <- 1
Microtia <- rep(0, length(seurat.child$type))
Microtia[seurat.child$type == 'Microtia'] <- 1
edesign <- cbind(Time,Replicates,Normal,Microtia)
# rownames(edesign) = paste("Cell",c(1:nrow(edesign)), sep= "_")
edesign <- data.frame(Time,Replicates,Normal,Microtia,
                      row.names = colnames(mat.count))
d = make.design.matrix(edesign, degree = 1)

# filter genes
path.M123 <- '/home/disk/drizzle/wgk/microtia_chon_child_M1M2M3/'
file.gsea.marker <- paste0(path.M123, 'marker_GSEA.Rdata')
list.marker.gsea <- readRDS(file.gsea.marker)
diff.genes <- c()
for (cell in names(list.marker.gsea)) {
    sub.marker <- list.marker.gsea[[cell]]
    sub.marker.1 <- sub.marker[sub.marker$p_val_adj < 0.01 & 
                                   abs(sub.marker$avg_logFC) > 0.25,]
    diff.genes <- c(diff.genes, rownames(sub.marker.1))
}
diff.genes <- unique(diff.genes)

high.exp <- function(mat.gene, gene) {
    vec.single <- mat.gene[gene,]
    agge.group <- aggregate(vec.single, by = list(cond_times), FUN = mean)
    if (max(agge.group$x) > 0.25) {
        return(data.frame(gene = gene, high.exp = 1))
    } else {
        return(data.frame(gene = gene, high.exp = 0))
    }
}

library(foreach)
library(doParallel)
registerDoParallel(cores = 10)
cond_times <- paste(seurat.child$type, Time, sep = '_')
df.highexp <- 
    foreach(gene = diff.genes, .combine = rbind) %dopar% high.exp(mat.gene, gene)
genes.highexp <- df.highexp[df.highexp$high.exp == 1, 'gene']

# maSigPro
fit = p.vector(mat.count[genes.highexp,], d, Q=0.01, 
               MT.adjust = "bonferroni",counts = T)
file.fit <- paste0(path.change, 'fit.Rdata')
saveRDS(fit, file.fit)

fit <- readRDS(file.fit)
tstep = T.fit(fit,step.method="backward")
file.tstep <- paste0(path.change, 'tstep.Rdata')
saveRDS(tstep, file.tstep)

tstep <- readRDS(file.tstep)

sigs = get.siggenes(tstep, rsq = 0.25, vars = "groups")

sigs.new <- sigs$sig.genes$MicrotiavsNormal
new.genes <- intersect(rownames(sigs.new$sig.profiles), rownames(mat.gene))
length(new.genes)

# XY genes
df.genes.v19 <- read.delim(
    '/home/yzj/publicData/annotation/hg19/gencode_v19_gene_pos.txt',
    header = F)
XY.genes <- df.genes.v19$V1[df.genes.v19$V2 %in% c('chrY', 'chrX')]
new.genes <- setdiff(new.genes, XY.genes)
length(new.genes)


sigs.new$sig.profiles <- sigs.new$sig.profiles[new.genes,]
sigs.new$coefficients <- sigs.new$coefficients[new.genes,]
sigs.new$group.coeffs <- sigs.new$group.coeffs[new.genes,]
sigs.new$sig.pvalues <- sigs.new$sig.pvalues[new.genes,]

# df.rsq <- tstep$sol
df.pval <- sigs.new$sig.pvalues
df.coff <- sigs.new$coefficients

pdf(file = paste(path.change,"cluster_gene_1.pdf",sep=''), width = 20, height = 20);
res.cluster <- see.genes(
    sigs$sig.genes$MicrotiavsNormal,
    show.fit = T,
    dis =d$dis,
    cluster.method="hclust" ,
    cluster.data = 1,
    k = 12,
    # k.mclust = T,
    ylim = c(0, 20),
    newX11 = F)
dev.off()
file.cluster <- paste0(path.change, 'cluster.Rdata')
saveRDS(res.cluster, file.cluster)
res.cluster <- readRDS(file.cluster)

View(merge(res.cluster$cut, df.pval, by = 'row.names'))
View(merge(res.cluster$cut, df.coff, by = 'row.names'))


# plot by clust
path.clust <- paste0(path.change, 'clust_plot/')
if (!file.exists(path.clust)) {
    dir.create(path.clust)
}
df.pc.gene.sel <- df.pc.gene[, c('idx', 'Harmony2', 'celltype', 'status', new.genes)]
clusters <- names(table(res.cluster$cut))
sort.plot.clusters <- c('Module_3', 'Module_12', 'Module_11', 
                        'Module_7', 'Module_4', 'Module_5',
                        'Module_1', 'Module_2', 'Module_9', 
                        'Module_8', 'Module_10', 'Module_6')
sort.plot.Modules <- c('Module_3', 'Module_4', 'Module_9',
                        'Module_12', 'Module_5', 'Module_8',
                        'Module_11', 'Module_1', 'Module_10',
                        'Module_7', 'Module_2', 'Module_6')
colors <- c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
names(colors) <- c('CSC', 'TC', 'C1', 'C2')

list.plot <- list()
for (cluster in sort.plot.clusters) {
    cluster_id <- as.numeric(strsplit(cluster, split = '_')[[1]][2])
    sub.gene <- intersect(names(res.cluster$cut[res.cluster$cut == cluster_id]), new.genes)
    if (length(sub.gene) < 1) {next()}
    sub.df.gene <- df.pc.gene.sel[,sub.gene]
    sub.df <- data.frame(idx = df.pc.gene.sel$idx, status = df.pc.gene.sel$status, 
                         celltype = df.pc.gene.sel$celltype, 
                         mean.exp = rowMeans(sub.df.gene))
    p.mean <-
        ggplot(data = sub.df, aes(x = idx, y = mean.exp, linetype = status)) + 
        geom_point(aes(color = celltype), size = 0.3) + 
        scale_color_discrete(breaks = c('CSC', 'TC', 'C1', 'C2'),
                             labels = c('CSC', 'C0', 'C1', 'C2')) + 
        scale_color_manual(labels = c('CSC', 'TC', 'C1', 'C2'),
                           values = colors) + 
        # xlim(-30, 10) + 
        geom_smooth(color = '#696969') + 
        labs(x = '', y = '') + 
        theme(panel.background=element_rect(fill='transparent', color='black',size = 1),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.position = 'none') +
        theme(panel.background=element_rect(fill='transparent', color='black',size = 1)) +
        annotate('text', label = cluster, 
                 x = 22000, y = max(sub.df$mean.exp), 
                 hjust = 1, vjust = 1, size = 7, family = 'Arial')
    ggsave(p.mean, path = path.clust, 
           filename = paste0(cluster, '.png'),
           height = 6, width = 9, units = 'cm')
    list.plot[[cluster]] <- p.mean
}


library(gridExtra)
p.merge <- marrangeGrob(list.plot, ncol = 4, nrow = 3, top = NULL)
ggsave(plot = p.merge, path = path.change, 
       filename = 'clust_plot.png',
       height = 20, width = 30, units = 'cm')


# plot
Time <- as.character(seurat.child$celltype)
Harmony2 <- seurat.child@reductions$harmony@cell.embeddings[, 'harmony_2']
quantile.H2 <- quantile(Harmony2, probs = seq(0, 1, 1/101))
for (i in 0:100) {
    Time[Harmony2 >= quantile.H2[i+1] & Harmony2 <= quantile.H2[i+2]] <- i
}
Time <- as.numeric(Time)

# mat.sel <- seurat.child@assays$RNA@scale.data[new.genes,]
# length(new.genes)
# dim(mat.sel)
mat.gene <- seurat.child@assays$RNA@data
df.sig.genes <- data.frame(t(as.matrix(mat.gene[new.genes,])), check.names = F)
df.sig.genes$Harmony2 <- Harmony2
df.sig.genes$celltype <- seurat.child$celltype
df.sig.genes$status <- seurat.child$type
df.sig.genes$Time.point <- Time
df.sig.genes <- df.sig.genes[order(Harmony2, decreasing = F),]
df.sig.genes$idx <- 1:nrow(df.sig.genes)

df.time <- data.frame()
df.time.status <- data.frame()
for (cluster in clusters) {
    sub.gene <- intersect(names(res.cluster$cut[res.cluster$cut == cluster]), new.genes)
    vec.single <- rowMeans(df.sig.genes[,sub.gene])
    agge.group <- aggregate(vec.single, by = list(df.sig.genes$Time.point), FUN = mean)
    sub.df <- data.frame(Time.point = agge.group$Group.1, 
                         cluster = rep(paste0('Module_', cluster), nrow(agge.group)), 
                         mean.exp = agge.group$x)
    df.time <- rbind(df.time, sub.df)
    for (st in unique(df.sig.genes$status)) {
        vec.single.st <- vec.single[df.sig.genes$status == st]
        vec.time.st <- df.sig.genes$Time.point[df.sig.genes$status == st]
        agge.group.st <- aggregate(vec.single.st, by = list(vec.time.st), FUN = mean)
        sub.df.st <- data.frame(Time.point = agge.group.st$Group.1, 
                             cluster = rep(paste0('Module_', cluster), nrow(agge.group.st)), 
                             status = rep(st, nrow(agge.group.st)),
                             mean.exp = agge.group.st$x)
        df.time.status <- rbind(df.time.status, sub.df.st)
    }
}

mat.plot <- reshape2::dcast(df.time, cluster ~ Time.point, value.var = 'mean.exp')
row.names(mat.plot) <- mat.plot$cluster
mat.plot$cluster <- NULL
mat.plot <- t(scale(t(as.matrix(mat.plot))))
# no status / clust tree
bk <- c(seq(-10,-0.1,by=0.01),seq(0,10,by=0.01))
pheatmap::pheatmap(mat.plot,
                   color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
                   cluster_rows = T, cluster_cols = F, scale = "none",
                   display_numbers = F,
                   # annotation_col = annotation_col, annotation_colors = ann_colors,
                   show_rownames = T, show_colnames = F, legend = T, 
                   fontsize_row = 10, 
                   # labels_row = labels_row,
                   # gaps_col = c(4), 
                   filename = paste0(path.change, 'heatmap_time_change.png'), 
                   width = 10, height = 6
)

# add status
sort.clusters <- c('Module_3', 'Module_12', 'Module_11', 
                   'Module_7', 'Module_4', 'Module_5',
                   'Module_1', 'Module_2', 'Module_9', 
                   'Module_8', 'Module_10', 'Module_6')
sort.row <- c()
for (sub in rev(sort.clusters)) {
    sort.row <- c(sort.row, c(paste0('Microtia_', sub), paste0('Normal_', sub)))
}
df.time.status$row_name <- paste(df.time.status$Time.point, df.time.status$status, sep = '_')
mat.heatmap <- reshape2::dcast(df.time.status, row_name ~ cluster, value.var = 'mean.exp')
row.names(mat.heatmap) <- mat.heatmap$row_name
mat.heatmap$row_name <- NULL
mat.heatmap <- scale(as.matrix(mat.heatmap))
mat.heatmap[mat.heatmap > 2.95] <- 2.95
df.heatmap <- reshape2::melt(mat.heatmap)
names(df.heatmap) <- c('row_name', 'cluster', 'MeanExp')
df.heatmap <- merge(df.heatmap, df.time.status, by = c('cluster', 'row_name'))
df.heatmap$cluster <- factor(df.heatmap$cluster, levels = sort.clusters)
df.heatmap$status <- factor(df.heatmap$status, levels = c("Microtia", "Normal"))
df.heatmap$Time.point <- as.factor(df.heatmap$Time.point)
df.heatmap$st_cluster <- paste(df.heatmap$status, df.heatmap$cluster, sep = '_')
df.heatmap$st_cluster <- factor(df.heatmap$st_cluster, levels = sort.row)
plot.heatmap <- 
    ggplot(data = df.heatmap, aes(x = Time.point, y = st_cluster)) + 
    geom_tile(aes(fill = MeanExp)) + 
    facet_grid(cluster ~ ., scales = 'free', space = 'free') +
    scale_fill_gradientn(
        colors = colorRampPalette(rev(brewer.pal(n = 11, name =  "Spectral"))[1:9])(100)) +
    scale_y_discrete(position = 'left') + 
    scale_x_discrete(breaks = c(0, 25, 50, 75, 100)) + 
    coord_cartesian(xlim=c(0, 100)) + 
    labs(x = 'Pseudotime', y = '', fill = 'Scaled \nExpression') +
    theme(panel.background = element_rect(color = 'transparent',
                                          fill = 'transparent'), 
          axis.text.x = element_text(size = 12, color = "black", 
                                     face = "bold", family = 'Arial'),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.x = element_text(
              size = 12, color = "black", face = "bold"), 
          strip.background = element_rect(
              color = "transparent", fill = "transparent"),
          strip.text.y = element_text(angle = 0, size = 12, color = "black", 
                                      face = "bold", family = 'Arial'),
          panel.spacing = unit(0.1, "cm"),
          strip.placement = "inside",
          axis.title.x = element_text(
              size = 14, color = "black", face = "bold", family = 'Arial'),
          legend.title = element_text(
              size = 13, color = "black", face = "bold", family = 'Arial'),
          legend.text = element_text(
              size = 12, color = "black", face = "bold", family = 'Arial'), 
          legend.position = 'right',
          plot.margin = unit(c(0, 0.6, 0.6, -0.6),"cm"))
ggsave(filename = 'heatmap_time_change_status.png', 
       path = path.change, plot = plot.heatmap,
       units = 'cm', height = 12, width = 18)

# library(tidyverse)
# library(aplot)
library(patchwork)
df.status <- data.frame(status=rep(c("Microtia","Normal"), times=12),
                        p = rep('white', 24), row_name = sort.row, 
                        cluster = rep(rev(sort.clusters), c(rep(2, 12))))
df.status$row_name <- factor(df.status$row_name, levels = (sort.row))
df.status$cluster <- factor(df.status$cluster, levels = sort.clusters)
df.status$status <- factor(df.status$status, levels = c("Microtia", "Normal"))
plot.Status <- 
    ggplot(data = df.status, aes(x=p,y=row_name,fill=status))+
    geom_tile() + 
    facet_grid(cluster ~ ., scales = 'free', space = 'free') +
    scale_y_discrete(position="right") +
    scale_fill_manual(values = c("#637FBF", "#6C6C6C")) + 
    xlab(NULL) + ylab(NULL) +
    theme(panel.background = element_rect(color = 'transparent',
                                          fill = 'transparent'),
          strip.text.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.spacing = unit(0.1, "cm"),
          legend.title = element_text(
              size = 13, color = "black", face = "bold", family = 'Arial'),
          legend.text = element_text(
              size = 12, color = "black", face = "bold", family = 'Arial'),
          legend.position = 'left',
          plot.margin = unit(c(0, -1, 0.6, 0.6),"cm"))+
    labs(fill = "Status")
ggsave(filename = 'heatmap_bar.png', 
       path = path.change, plot = plot.Status,
       units = 'cm', height = 12, width = 4.2)

plot.final <- plot.Status + plot.heatmap + plot_layout(widths = c(1, 30),
                                         guides = 'collect')

ggsave(filename = 'heatmap_time_change_final.png', 
       path = path.change, plot = plot.final,
       units = 'cm', height = 12, width = 24)

colors <- c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
names(colors) <- c('CSC', 'TC', 'C1', 'C2')
plot.dens <- 
    ggplot(df.sig.genes, aes(x = Time.point, color = celltype, fill = celltype)) + 
    geom_density(alpha = 0.3) + 
    scale_color_manual(breaks = names(colors), values = colors) + 
    scale_color_discrete(breaks = c('CSC', 'TC', 'C1', 'C2'),
                         labels = c('CSC', 'C0', 'C1', 'C2')) + 
    scale_fill_discrete(breaks = c('CSC', 'TC', 'C1', 'C2'),
                         labels = c('CSC', 'C0', 'C1', 'C2')) + 
    labs(x = '', y = 'Probability density', 
         color = 'Cell Type', fill = 'Cell Type') + 
    theme(panel.background = element_rect(color = 'transparent',
                                          fill = 'transparent'),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(
              size = 10, color = "black", face = "bold", family = 'Arial',
              margin=margin(0,-20,0,0)),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_text(
              size = 13, color = "black", face = "bold", family = 'Arial'),
          legend.text = element_text(
              size = 12, color = "black", face = "bold", family = 'Arial'),
          plot.margin = unit(c(0.6, 0, 0, 0.6),"cm"))
ggsave(filename = 'heatmap_density.png', 
       path = path.change, plot = plot.dens,
       units = 'cm', height = 4, width = 20)

# percentage of TF
df.TF <- read.table('/home/disk/drizzle/DataBase/Human_TF/TF_names_v_1.01.txt')
length(intersect(new.genes, df.TF$V1))
percent.TF <- data.frame()
for (cluster in clusters) {
    sub.gene <- intersect(names(res.cluster$cut[res.cluster$cut == cluster]), new.genes)
    percent.TF <- 
        rbind(percent.TF, 
              data.frame(Cluster = paste0('Module_', cluster),
                         Percent = length(intersect(sub.gene, df.TF$V1))/length(sub.gene)))
}
percent.TF$Cluster <- factor(percent.TF$Cluster, levels = rev(sort.clusters))
plot.TF <- 
    ggplot(percent.TF, aes(x = Cluster, y = Percent)) + 
    geom_bar(stat = 'identity', color = 'transparent', fill = '#3CB371') + 
    labs(y = 'Percentage of TF') + 
    scale_y_continuous(labels = c(0, 0.2, 0.4, 0.6)) + 
    theme_classic() + 
    coord_flip() + 
    theme(panel.background = element_rect(color = 'transparent',
                                          fill = 'transparent'),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(
              size = 14, color = "black", face = "bold", family = 'Arial'),
          axis.text.x = element_text(
              size = 12, color = "black", face = "bold", family = 'Arial'),
          plot.margin = unit(c(0.6, 0.6, 0.6, -0.1),"cm"))
ggsave(filename = 'heatmap_percent_TF.png', 
       path = path.change, plot = plot.TF,
       units = 'cm', height = 12.5, width = 5)


# GO enrich per cluster
path.GO <- paste0(path.change, 'GO/')
if (!file.exists(path.GO)) {
    dir.create(path.GO)
}

library(org.Hs.eg.db)
library(clusterProfiler)
all.genes <- rownames(mat.count)
use.genes <- intersect(keys(org.Hs.eg.db, keytype = "SYMBOL"), all.genes)
list.go.BP <- list()
list.go.MF <- list()
list.go.CC <- list()
list.go.BP.all <- list()
list.go.MF.all <- list()
list.go.CC.all <- list()
for (cluster in clusters) {
    type <- paste0('Module_', cluster)
    sub.gene <- intersect(names(res.cluster$cut[res.cluster$cut == cluster]), new.genes)
    genes.input <- intersect(sub.gene, use.genes)
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
}

file.go.BP.all <- paste0(path.GO, 'GO_BP_all.Rdata')
saveRDS(list.go.BP.all, file = file.go.BP.all)
file.go.MF.all <- paste0(path.GO, 'GO_MF_all.Rdata')
saveRDS(list.go.MF.all, file = file.go.MF.all)
file.go.CC.all <- paste0(path.GO, 'GO_CC_all.Rdata')
saveRDS(list.go.CC.all, file = file.go.CC.all)
list.go.BP.all <- readRDS(file.go.BP.all)

file.go.BP <- paste0(path.GO, 'GO_BP.Rdata')
saveRDS(list.go.BP, file = file.go.BP)
file.go.MF <- paste0(path.GO, 'GO_MF.Rdata')
saveRDS(list.go.MF, file = file.go.MF)
file.go.CC <- paste0(path.GO, 'GO_CC.Rdata')
saveRDS(list.go.CC, file = file.go.CC)

# GO plot
list.sel.GO <- list()
list.sel.GO$Module_3 <- c('cellular response to reactive oxygen species')
list.sel.GO$Module_11 <- c('p38MAPK cascade')
list.sel.GO$Module_7 <- c('response to unfolded protein', 
                           'extrinsic apoptotic signaling pathway')
list.sel.GO$Module_4 <- c('negative regulation of protein phosphorylation')
list.sel.GO$Module_5 <- c('cellular response to zinc ion', 
                           'cellular response to copper ion')
list.sel.GO$Module_1 <- c('cartilage development')
list.sel.GO$Module_2 <- c('extracellular matrix organization', 
                           'chondrocyte differentiation',
                           'cartilage development')
list.sel.GO$Module_9 <- c('cellular transition metal ion homeostasis')
list.sel.GO$Module_8 <- c('chemokine-mediated signaling pathway')
list.sel.GO$Module_6 <- c('interferon-gamma-mediated signaling pathway', 
                           'regulation of inflammatory response',
                           'extracellular matrix disassembly')

# bubble plot
sort.cells <- names(list.sel.GO)
colors <- brewer.pal(10,"Paired")
df.plot <- data.frame()
i = 0
GOterms <- c()
for (cell in sort.cells) {
    i = i + 1
    sub.go <- list.go.BP.all[[cell]]
    sel.go.term <- list.sel.GO[[cell]]
    sel.go <- sub.go[sub.go$Description %in% sel.go.term, 
                     c('Description', 'pvalue')]
    sel.go$log10Pval <- -log10(sel.go$pvalue)
    sel.go$celltype <- rep(cell, nrow(sel.go))
    # sel.go$Description <- factor(sel.go$Description, levels = rev(sel.go.term))
    df.plot <- rbind(df.plot, sel.go)
    GOterms <- c(GOterms, sel.go.term)
}
col_name <- paste(df.plot$celltype, df.plot$Description, sep = '_')
df.plot$col_name <- factor(col_name, levels = rev(col_name))
df.plot$celltype <- factor(df.plot$celltype, levels = sort.cells)


p <- ggplot(df.plot, aes(x = celltype, y = col_name, 
                         color = celltype, size = log10Pval)) + 
    geom_point(fill = 'cornsilk') +
    scale_color_manual(breaks = sort.cells,
                       values = colors) + 
    scale_size_continuous(range = c(3,5)) +
    scale_y_discrete(breaks = col_name, labels = GOterms) +
    labs(x = '', y = 'GO term', color = 'Module',
         size = expression(paste("-log"[10], "(", italic("P"), "-value)"))) +
    theme(panel.background = element_rect(color = 'black',
                                          fill = 'transparent'), 
          panel.grid.major = element_line(colour = 'gray', size = 0.2, linetype = 5),
          axis.title = element_text(size = 14, face = 'bold', 
                                    color = 'black', family = 'Arial'), 
          axis.text.y = element_text(size = 12, face = 'bold', 
                                     color = 'black', family = 'Arial'), 
          axis.text.x = element_text(size = 12, face = 'bold', 
                                     color = 'black', family = 'Arial',
                                     angle = 45, hjust = 1, vjust = 1),
          legend.text = element_text(size = 12, color = 'black', family = 'Arial'),
          legend.title = element_text(size = 14, face = 'bold', 
                                      color = 'black', family = 'Arial'),
          legend.key = element_blank()) + 
    guides(colour = guide_legend(override.aes = list(size=4)))
ggsave(plot = p, path = path.change, 
       filename = paste0('GO.png'),
       height = 15, width = 24, units = 'cm')
