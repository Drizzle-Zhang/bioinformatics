setwd('/home/zy/my_git/bioinformatics/wgk')
.libPaths('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)
library(ggplot2)
require("RColorBrewer")
library(SCENIC)


file.chon <- '/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Control_Chond.Rdata'
seurat.normal <- readRDS(file.chon)
# seurat.normal <- NormalizeData(seurat.normal)
# seurat.normal <- FindVariableFeatures(seurat.normal, nfeatures = 3000)
# highvar.genes <- VariableFeatures(seurat.normal)

path.M123 <- '/home/disk/drizzle/wgk/age/'
path.change <- paste0(path.M123, 'chon_lineage_1/')
if (!file.exists(path.change)) {
    dir.create(path.change)
}


seurat.normal.N <- subset(seurat.normal, subset = type == 'Normal')
seurat.normal.M <- subset(seurat.normal, subset = type == 'Microtia')

library(foreach)
library(doParallel)
registerDoParallel(cores = 10)
corr.PC.exp <- function(mat.gene, PC2, gene) {
    vec.exp <- mat.gene[gene,]
    corr <- cor(PC2, vec.exp, method = 'spearman')
    return(data.frame(gene = gene, corr = corr))
}

# Normal
Harmony2_N <- seurat.normal.N@reductions$harmony@cell.embeddings[, 'harmony_2']
mat.gene_N <- seurat.normal.N@assays$RNA@data
df.corr_N <- 
    foreach(gene = highvar.genes, .combine = rbind) %dopar% 
    corr.PC.exp(mat.gene_N, Harmony2_N, gene)
rownames(df.corr_N) <- highvar.genes

# Microtia
Harmony2_M <- seurat.normal.M@reductions$harmony@cell.embeddings[, 'harmony_2']
mat.gene_M <- seurat.normal.M@assays$RNA@data
df.corr_M <- 
    foreach(gene = highvar.genes, .combine = rbind) %dopar% 
    corr.PC.exp(mat.gene_M, Harmony2_M, gene)
rownames(df.corr_M) <- highvar.genes

df.corr <- merge(df.corr_N, df.corr_M, by = 'gene')
df.corr[is.na(df.corr)] <- 0
df.corr$diff_M_N <- df.corr$corr.y - df.corr$corr.x


# plot single gene
Harmony2 <- seurat.normal@reductions$harmony@cell.embeddings[, 'harmony_2']
mat.gene <- seurat.normal@assays$RNA@data
df.pc.gene <- data.frame(t(as.matrix(mat.gene)))
df.pc.gene$Harmony2 <- Harmony2
df.pc.gene$celltype <- seurat.normal$celltype
df.pc.gene$status <- seurat.normal$Phase
df.pc.gene <- df.pc.gene[order(Harmony2, decreasing = F),]
df.pc.gene$idx <- 1:nrow(df.pc.gene)
colors <- c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
names(colors) <- c('CSC', 'TC', 'C1', 'C2')

gene <- 'JUNB'
df.plot <- df.pc.gene[, c('idx', 'Harmony2', 'celltype', 'status', gene)]
names(df.plot) <- c('idx', 'Harmony2', 'celltype', 'status', 'gene')
# p.gene <-
    ggplot(data = df.plot, aes(x = idx, y = gene, linetype = status)) + 
    geom_point(aes(color = celltype), size = 0.3) + 
    scale_color_manual(labels = c('CSC', 'TC', 'C1', 'C2'),
                       values = colors) + 
    # xlim(-30, 10) + 
    geom_smooth(color = '#696969') +
    # geom_smooth(color = '#696969', formula = y~poly(x, 3), method = lm) + 
    labs(x = '', y = '') + 
    # theme(panel.background=element_rect(fill='transparent', color='black',size = 1),
    #       axis.text = element_blank(),
    #       axis.ticks = element_blank(),
    #       legend.position = 'none') +
    theme(panel.background=element_rect(fill='transparent', color='black',size = 1)) +
    annotate('text', label = gene, x = 22000, y = max(df.plot$gene), 
             hjust = 1, vjust = 1, size = 7)
    
    
# time series
# four time
library(maSigPro)
mat.count <- seurat.normal@assays$RNA@counts[highvar.genes,]
Time <- as.character(seurat.normal$celltype)
Time[Time == 'CSC'] <- 1
Time[Time == 'TC'] <- 2
Time[Time == 'C1'] <- 3
Time[Time == 'C2'] <- 4

# 100 points
library(maSigPro)
mat.count <- seurat.normal@assays$RNA@counts
Time <- as.character(seurat.normal$celltype)
Harmony2 <- seurat.normal@reductions$harmony@cell.embeddings[, 'harmony_2']
quantile.H2 <- quantile(Harmony2, probs = seq(0, 1, 0.01))
for (i in 1:100) {
    Time[Harmony2 >= quantile.H2[i] & Harmony2 <= quantile.H2[i+1]] <- i
}

# time series analysis
Time <- as.numeric(Time)
Replicates <- as.numeric(as.factor(seurat.normal$batch))
Children <- rep(0, length(seurat.normal$type))
Children[seurat.normal$Phase == 'Children'] <- 1
Adults <- rep(0, length(seurat.normal$type))
Adults[seurat.normal$Phase == 'Adults'] <- 1
edesign <- cbind(Time,Replicates,Children,Adults)
# rownames(edesign) = paste("Cell",c(1:nrow(edesign)), sep= "_")
edesign <- data.frame(Time,Replicates,Children,Adults,
                      row.names = colnames(mat.count))
d = make.design.matrix(edesign, degree = 1)

# filter genes
file.gsea.marker <- '/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_inChildrenAdults_GSEA.RDS'
list.marker.gsea <- readRDS(file.gsea.marker)
diff.genes <- c()
for (cell in names(list.marker.gsea)) {
    sub.marker <- list.marker.gsea[[cell]]
    sub.marker.1 <- sub.marker[sub.marker$pct.1 > 0.1 | sub.marker$pct.2 > 0.1,]
    sub.marker.1 <- sub.marker[sub.marker$p_val_adj < 0.01 & 
                                   abs(sub.marker$avg_logFC) > 0.3,]
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
mat.gene <- seurat.normal@assays$RNA@data
cond_times <- paste(seurat.normal$type, Time, sep = '_')
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

sigs.new <- sigs$sig.genes$AdultsvsChildren
new.genes <- intersect(rownames(sigs.new$sig.profiles), rownames(mat.gene))
sigs.new$sig.profiles <- sigs.new$sig.profiles[new.genes,]
sigs.new$coefficients <- sigs.new$coefficients[new.genes,]
sigs.new$group.coeffs <- sigs.new$group.coeffs[new.genes,]
sigs.new$sig.pvalues <- sigs.new$sig.pvalues[new.genes,]
length(new.genes)

# df.rsq <- tstep$sol
df.pval <- sigs.new$sig.pvalues
df.coff <- sigs.new$coefficients

pdf(file = paste(path.change,"cluster_gene_1.pdf",sep=''), width = 20, height = 20);
res.cluster <- see.genes(
    sigs$sig.genes$AdultsvsChildren,
    show.fit = T,
    dis =d$dis,
    cluster.method="hclust" ,
    cluster.data = 1,
    k = 9,
    # k.mclust = T,
    ylim = c(0, 20),
    newX11 = F)
dev.off()
file.cluster <- paste0(path.change, 'cluster.Rdata')
saveRDS(res.cluster, file.cluster)


View(merge(res.cluster$cut, df.pval, by = 'row.names'))
View(merge(res.cluster$cut, df.coff, by = 'row.names'))


# plot by clust
df.pc.gene.sel <- df.pc.gene[, c('idx', 'Harmony2', 'celltype', 'status', new.genes)]
clusters <- names(table(res.cluster$cut))
for (cluster in clusters) {
    sub.gene <- intersect(names(res.cluster$cut[res.cluster$cut == cluster]), new.genes)
    if (length(sub.gene) < 1) {next()}
    sub.df.gene <- df.pc.gene.sel[,sub.gene]
    sub.df <- data.frame(idx = df.pc.gene.sel$idx, status = df.pc.gene.sel$status, 
                         celltype = df.pc.gene.sel$celltype, 
                         mean.exp = rowMeans(sub.df.gene))
    p.mean <-
        ggplot(data = sub.df, aes(x = idx, y = mean.exp, linetype = status)) + 
        geom_point(aes(color = celltype), size = 0.3) + 
        scale_color_manual(labels = c('CSC', 'TC', 'C1', 'C2'),
                           values = colors) + 
        # xlim(-30, 10) + 
        geom_smooth(color = '#696969') + 
        labs(x = '', y = '') + 
        # theme(panel.background=element_rect(fill='transparent', color='black',size = 1),
        #       axis.text = element_blank(),
        #       axis.ticks = element_blank(),
        #       legend.position = 'none') +
        theme(panel.background=element_rect(fill='transparent', color='black',size = 1)) +
        annotate('text', label = paste0('Cluster_', cluster), 
                 x = 22000, y = max(sub.df$mean.exp), 
                 hjust = 1, vjust = 1, size = 7)
    ggsave(p.mean, path = path.change, 
           filename = paste0('Cluster_', cluster, '.png'),
           height = 6, width = 9, units = 'cm')
}


# plot
Time <- as.character(seurat.normal$celltype)
Harmony2 <- seurat.normal@reductions$harmony@cell.embeddings[, 'harmony_2']
quantile.H2 <- quantile(Harmony2, probs = seq(0, 1, 1/101))
for (i in 0:100) {
    Time[Harmony2 >= quantile.H2[i+1] & Harmony2 <= quantile.H2[i+2]] <- i
}
Time <- as.numeric(Time)

# mat.sel <- seurat.normal@assays$RNA@scale.data[new.genes,]
# length(new.genes)
# dim(mat.sel)
mat.gene <- seurat.normal@assays$RNA@data
df.sig.genes <- data.frame(t(as.matrix(mat.gene[new.genes,])), check.names = F)
df.sig.genes$Harmony2 <- Harmony2
df.sig.genes$celltype <- seurat.normal$celltype
df.sig.genes$status <- seurat.normal$Phase
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
                         cluster = rep(paste0('Cluster_', cluster), nrow(agge.group)), 
                         mean.exp = agge.group$x)
    df.time <- rbind(df.time, sub.df)
    for (st in unique(df.sig.genes$status)) {
        vec.single.st <- vec.single[df.sig.genes$status == st]
        vec.time.st <- df.sig.genes$Time.point[df.sig.genes$status == st]
        agge.group.st <- aggregate(vec.single.st, by = list(vec.time.st), FUN = mean)
        sub.df.st <- data.frame(Time.point = agge.group.st$Group.1, 
                                cluster = rep(paste0('Cluster_', cluster), nrow(agge.group.st)), 
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
                   clustering_method = 'average',
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
sort.clusters <- c('Cluster_6', 'Cluster_7', 'Cluster_8', 
                   'Cluster_2', 'Cluster_5', 'Cluster_9',
                   'Cluster_1', 'Cluster_3', 'Cluster_4')
sort.row <- c()
for (sub in rev(sort.clusters)) {
    sort.row <- c(sort.row, c(paste0('Adults_', sub), paste0('Children_', sub)))
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
df.heatmap$status <- factor(df.heatmap$status, levels = c("Adults", "Children"))
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
          axis.text.x = element_text(size = 12, color = "black", face = "bold"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.x = element_text(
              size = 12, color = "black", face = "bold"), 
          strip.background = element_rect(
              color = "transparent", fill = "transparent"),
          strip.text.y = element_text(angle = 0, size = 12, color = "black", face = "bold"),
          panel.spacing = unit(0.1, "cm"),
          strip.placement = "inside",
          axis.title.x = element_text(
              size = 14, color = "black", face = "bold"),
          legend.title = element_text(
              size = 13, color = "black", face = "bold"),
          legend.text = element_text(
              size = 12, color = "black", face = "bold"), 
          legend.position = 'right',
          plot.margin = unit(c(0, 0.6, 0.6, -0.6),"cm"))
ggsave(filename = 'heatmap_time_change_status.png', 
       path = path.change, plot = plot.heatmap,
       units = 'cm', height = 9, width = 20)

library(tidyverse)
library(aplot)
library(patchwork)
df.status <- data.frame(status=rep(c("Microtia","Normal"), times=9),
                        p = rep('white', 18), row_name = sort.row, 
                        cluster = rep(rev(sort.clusters), c(rep(2, 9))))
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
              size = 13, color = "black", face = "bold"),
          legend.text = element_text(
              size = 12, color = "black", face = "bold"),
          legend.position = 'left',
          plot.margin = unit(c(0, 0, 0.6, 0.6),"cm"))+
    labs(fill = "Status")
ggsave(filename = 'heatmap_bar.png', 
       path = path.change, plot = plot.Status,
       units = 'cm', height = 12, width = 4.2)

plot.final <- plot.Status + plot.heatmap + plot_layout(widths = c(1, 40),
                                                       guides = 'collect')

ggsave(filename = 'heatmap_time_change_final.png', 
       path = path.change, plot = plot.final,
       units = 'cm', height = 11, width = 30)

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
              size = 10, color = "black", face = "bold",margin=margin(0,-20,0,0)),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_text(
              size = 13, color = "black", face = "bold"),
          legend.text = element_text(
              size = 12, color = "black", face = "bold"),
          plot.margin = unit(c(0.6, 0, 0, 0.6),"cm"))
ggsave(filename = 'heatmap_density.png', 
       path = path.change, plot = plot.dens,
       units = 'cm', height = 4, width = 27)

# percentage of TF
df.TF <- read.table('/home/disk/drizzle/DataBase/Human_TF/TF_names_v_1.01.txt')
length(intersect(new.genes, df.TF$V1))
percent.TF <- data.frame()
for (cluster in clusters) {
    sub.gene <- intersect(names(res.cluster$cut[res.cluster$cut == cluster]), new.genes)
    percent.TF <- 
        rbind(percent.TF, 
              data.frame(Cluster = paste0('Cluster_', cluster),
                         Percent = length(intersect(sub.gene, df.TF$V1))/length(sub.gene)))
}
percent.TF$Cluster <- factor(percent.TF$Cluster, levels = rev(sort.clusters))
plot.TF <- 
    ggplot(percent.TF, aes(x = Cluster, y = Percent)) + 
    geom_bar(stat = 'identity', color = 'transparent', fill = '#3CB371') + 
    labs(y = 'Percentage of TF') + 
    scale_y_continuous(breaks = c(0, 0.25, 0.5), labels = c(0, 0.25, 0.5)) + 
    theme_classic() + 
    coord_flip() + 
    theme(panel.background = element_rect(color = 'transparent',
                                          fill = 'transparent'),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(
              size = 14, color = "black", face = "bold"),
          axis.text.x = element_text(
              size = 12, color = "black", face = "bold"),
          plot.margin = unit(c(0.6, 0.6, 0.6, -0.1),"cm"))
ggsave(filename = 'heatmap_percent_TF.png', 
       path = path.change, plot = plot.TF,
       units = 'cm', height = 11.2, width = 5)









