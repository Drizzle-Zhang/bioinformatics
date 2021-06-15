
################################
## 1. EXP/TF AUC 变化趋势图
################################
library(Seurat)
library(ggplot2)
library(SCENIC)
require("RColorBrewer")
library(maSigPro)

seurat.chon <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Chond.Rdata')
seurat.child <- subset(seurat.chon, subset = batch %in% c('C4', 'C6', 'M1', 'M2', 'M3'))
rm(seurat.chon)

# plot single gene
#sample.cells <- sample(colnames(seurat.child),10000)
#sample.seurat.child <- subset(seurat.child,cells=sample.cells)
Harmony2 <- seurat.child@reductions$harmony@cell.embeddings[, 'harmony_2']
mat.gene <- seurat.child@assays$RNA@data
# AUC
regulonAUC <- readRDS(file='/home/yzj/JingMA_NEW/res/SCENIC_main/int/3.4_regulonAUC.Rds')
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), colnames(mat.gene)]
mat.auc <- as.matrix(regulonAUC@assays@data@listData$AUC)
df.pc.gene <- data.frame(t(rbind(as.matrix(mat.gene), mat.auc)), check.names = F)
df.pc.gene$Harmony2 <- Harmony2
df.pc.gene$celltype <- seurat.child$celltype
df.pc.gene$status <- seurat.child$type
df.pc.gene <- df.pc.gene[order(Harmony2, decreasing = F),]
df.pc.gene$idx <- 1:nrow(df.pc.gene)
colors <- c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
names(colors) <- c('CSC', 'C0', 'C1', 'C2')

########################
# TF AUC
########################
# vec.TF <- c('SOX8 (158g)','SOX5 (218g)', 'DBP (45g)')
# type='down'

vec.TF <- c('ATF3 (83g)','EGR1 (264g)', 'EGR3 (53g)')
type='up'
for (i in 1:length(vec.TF)) {
    TF <- vec.TF[i]
    df.plot <- df.pc.gene[, c('idx', 'Harmony2', 'celltype', 'status', TF)]
    names(df.plot) <- c('idx', 'Harmony2', 'celltype', 'status', 'TF')
    p1 <- ggplot(data = df.plot, aes(x = Harmony2, 
                                     linetype = status, 
                                     y = TF)) + 
        geom_point(aes(color = celltype), size = 0.0000001) + 
        scale_color_manual(labels = c('CSC', 'C0', 'C1', 'C2'),values = colors) + 
        xlim(-30, 10) +
        geom_smooth(color = '#696969',size=0.5) + theme_classic()+
        labs(x = '', y = '') + 
        theme(panel.background=element_rect(fill='transparent', color='black',size = 0.5),
              axis.text = element_blank(),axis.ticks = element_blank(),plot.margin = unit(c(0.1,0.1,-0.5,-0.5), "cm"),
              legend.position = 'none') +
        annotate('text', label = TF, x = 9, y = max(df.plot$TF), 
                 hjust = 1, vjust = 1, size = 2)
    if (i == 1) {
        p <- p1
    } else {
        p <- p / p1
    }
}
# ggsave(paste('JingMA_NEW/res/compMicrotia/MicrotiavsNormal_inChildren/FIG/lineage_AUC/Fig5A_TF_AUC_',type,'.pdf',sep=''),p,
#        height = 6, width = 3.5, units = 'cm')

ggsave(paste('JingMA_NEW/res/compMicrotia/MicrotiavsNormal_inChildren/FIG/lineage_AUC/Fig5B_TF_AUC_',type,'.pdf',sep=''),p,
       height = 6, width = 3.5, units = 'cm')

########################
## gene
########################
# vec.TF.exp <- c('SOX8','SOX5', 'DBP')
# type='down'

vec.TF.exp <- c('ATF3','EGR1', 'EGR3')
type='up'
for (i in 1:length(vec.TF.exp)) {
    TF <- vec.TF.exp[i]
    df.plot <- df.pc.gene[, c('idx', 'Harmony2', 'celltype', 'status', TF)]
    names(df.plot) <- c('idx', 'Harmony2', 'celltype', 'status', 'TF')
    p1 <- ggplot(data = df.plot, aes(x = Harmony2, 
                                     linetype = status, 
                                     y = TF)) + 
        geom_point(aes(color = celltype), size = 0.0000001) +  theme_classic()+
        scale_color_manual(labels = c('CSC', 'TC', 'C1', 'C2'),values = colors) + 
        xlim(-30, 10) +
        geom_smooth(color = '#696969',size=0.5) +
        labs(x = '', y = '') + 
        theme(panel.background=element_rect(fill='transparent', color='black',size = 0.5),plot.margin = unit(c(0.1,0.1,-0.5,-0.5), "cm"),
              axis.text = element_blank(),axis.ticks = element_blank(),legend.position = 'none') +
        annotate('text', label = TF, x = 9, y = max(df.plot$TF), 
                 hjust = 1, vjust = 1, size = 2)
    if (i == 1) {
        p <- p1
    } else {
        p <- p / p1
    }
}
# ggsave(paste('JingMA_NEW/res/compMicrotia/MicrotiavsNormal_inChildren/FIG/lineage_EXP/Fig5A_TF_EXP_',type,'.pdf',sep=''),p,
#        height = 6, width = 3.5, units = 'cm')
ggsave(paste('JingMA_NEW/res/compMicrotia/MicrotiavsNormal_inChildren/FIG/lineage_EXP/Fig5B_TF_EXP_',type,'.pdf',sep=''),p,
       height = 6, width = 3.5, units = 'cm')


################################
## 2. network图
################################
.libPaths('/home/zy/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)
library(patchwork)
library(tidyverse)
library(SCENIC)
library(igraph)
library(ggraph)
library(tidygraph)

path.TF.net <- '/home/disk/drizzle/wgk/TF_net/'

fc.cutoff <- 0.4
path.M123 <- '/home/disk/drizzle/wgk/microtia_chon_child_M1M2M3/'
path.cutoff <- paste0(path.M123, 'cutoff_', fc.cutoff, '/')
file.go.BP <- paste0(path.cutoff, 'GO_BP_all.Rdata')
list.go.BP <- readRDS(file.go.BP)
file.go.MF <- paste0(path.cutoff, 'GO_MF_all.Rdata')
list.go.MF <- readRDS(file.go.MF)

# down
file_CSC_down <- paste0(path.TF.net, 'igraph_CSC_down.RDS')
igraph_down_2 <- readRDS(file_CSC_down)

sel.MF_SCS_down <- c('antioxidant activity',
                     'S100 protein binding',
                     'extracellular matrix structural constituent')
sel.BP_SCS_down <- c('extracellular matrix organization',
                     'cellular response to zinc ion',
                     'skeletal system development',
                     'cartilage development')

BP_SCS_down <- list.go.BP[['CSC_Microtia_decrease']]
rownames(BP_SCS_down) <- BP_SCS_down$Description
MF_SCS_down <- list.go.MF[['CSC_Microtia_decrease']]
rownames(MF_SCS_down) <- MF_SCS_down$Description
df_GO_pre <- rbind(MF_SCS_down[sel.MF_SCS_down, c('Description', 'geneID')], 
                   BP_SCS_down[sel.BP_SCS_down, c('Description', 'geneID')])

vec_desc <- c('Oxidoreductase activity',
              'S100 protein binding',
              'Extracellular matrix organization',
              'Cellular response to zinc/copper ion',
              'Cartilage development')
vec_genes <- c(paste0(paste(df_GO_pre[1, 'geneID'], collapse = '/'), '/TXN'),
               paste(df_GO_pre[2, 'geneID'], collapse = '/'),
               paste0(paste(df_GO_pre[3:4, 'geneID'], collapse = '/'), '/COL11A1/'),
               paste(df_GO_pre[5, 'geneID'], collapse = '/'),
               paste0(paste(df_GO_pre[6:7, 'geneID'], collapse = '/'), '/WWP2/SCRG1/FMOD/CTGF'))
df_GO <- data.frame(Description = vec_desc, geneID = vec_genes)

Carti <- unlist(strsplit(df_GO[df_GO$Description == 'Cartilage development', 'geneID'],  '/'))
ECM <- unlist(strsplit(df_GO[df_GO$Description == 'Extracellular matrix organization', 'geneID'],  '/'))
Ion <- unlist(strsplit(df_GO[df_GO$Description == 'Cellular response to zinc/copper ion', 'geneID'],  '/'))
ROS <- unlist(strsplit(df_GO[df_GO$Description == 'Oxidoreductase activity', 'geneID'],  '/'))
S100 <- unlist(strsplit(df_GO[df_GO$Description == 'S100 protein binding', 'geneID'],  '/'))

library(RColorBrewer)
mycolor=brewer.pal(10,"Set3")
col=mycolor[c(1,3:6)]
group <- c('Carti','Ion','ECM','ROS','S100')
names(col) <- group

node_type <- V(igraph_down_2)$name
node_group <- rep('7', length(node_type))
node_group[node_type %in% c(Carti)] <- '1'
node_group[node_type %in% c(Ion)] <- '2'
node_group[node_type %in% c(ECM)] <- '3'
node_group[node_type %in% c(ROS)] <- '4'
node_group[node_type %in% c(S100)] <- '5'
node_group[node_type %in% c('SOX5', 'SOX8', 'DBP', 'ARID5B')] <- '6'
V(igraph_down_2)$group <- node_group

# gene expression FC
FCs <- V(igraph_down_2)$FC
names(FCs) <- V(igraph_down_2)$name
FC_TF <- FCs[names(FCs) %in% c('SOX5', 'SOX8', 'DBP', 'ARID5B')]
FC_alpha <- (scale(FC_TF) + 1.5) / 3
FCs[rownames(FC_alpha)] <- FC_alpha
FCs[!(names(FCs) %in% c('SOX5', 'SOX8', 'DBP', 'ARID5B'))] <- 1
V(igraph_down_2)$FC <- as.numeric(FCs)


ggraph_CSC_down <- as_tbl_graph(igraph_down_2)

plot_CSC_down <-
    ggraph(ggraph_CSC_down, layout = 'stress') +
    geom_edge_link(aes(edge_width=weight, alpha = weight),color="gray",
                   arrow = arrow(length = unit(2, 'mm')),
                   end_cap = circle(1, 'mm'), 
                   start_cap = circle(0.3, 'mm')) + 
    scale_edge_width(range=c(0.5,1)) +
    scale_edge_alpha(range=c(0.2,1)) + 
    scale_size_continuous(range = c(2,10)) +
    geom_node_point(aes(size = page_rank, fill = group, alpha = FC),
                    shape=21, color = 'transparent') +
    scale_color_manual(values = c(col,'#4169E1', 'gray')) + 
    scale_fill_manual(values = c(col,'#4169E1', 'gray')) + 
    # scale_alpha_manual(values = c(1,1,1,1,1,1,1, 0.1)) + 
    geom_node_text(aes(filter = (group %in% c(1, 2, 3, 4, 5)),label=name),size=2, repel = T) + 
    geom_node_text(aes(filter = group == 6,label=name),size=3) +
    theme_void() + 
    theme(legend.position = 'none') + 
    guides(size = F, edge_width = F, alpha = F)
# ggsave(filename = '/home/yzj/JingMA_NEW/res/compMicrotia/MicrotiavsNormal_inChildren/FIG/Fig5A_TF_CSC_down.pdf',
#        plot_CSC_down,height = 10, width = 10, units = 'cm')
ggsave(plot = plot_CSC_down, path = path.TF.net,
       filename = 'TF_CSC_down.png',
       height = 10, width = 10, units = 'cm')

# color bar
df_plot <- data.frame(FC = scale(FC_TF), TF = rownames(scale(FC_TF)))
df_plot <- rbind(df_plot, data.frame(FC = -1.5, TF = 'DOWN'))
df_plot <- rbind(df_plot, data.frame(FC = 1.5, TF = 'UP'))
df_plot$NUL <- rep('1', nrow(df_plot))

plot_bar <- 
    ggplot(data = df_plot, aes(x = TF, y = NUL, fill = FC)) + 
    geom_tile() + 
    scale_fill_gradient(low = 'transparent', high = '#4169E1', breaks = c(-1.5, 0, 1.5)) + 
    labs(fill = expression(paste("Scaled FC"['TF']))) + 
    theme(legend.title = element_text(size = 6, color = "black"),
          legend.text = element_text(size = 6, color = "black")) 
ggsave(plot = plot_CSC_down, path = path.TF.net,
       filename = 'TF_CSC_down_bar.png',
       height = 5, width = 5, units = 'cm')


## 查看gene所属类别
show.genes <- intersect(union(union(union(union(Carti,ECM),Ion),ROS),S100),node_type)
for(i in 1:length(show.genes)){
    g=show.genes[i]
    print(paste('!!!Gene: ',g))
    if(g %in% Carti){print('Carti')}
    if(g %in%  ECM){print('ECM')}
    if(g %in%  Ion){print('Ion')}
    if(g %in%  ROS){print('ROS')}
    if(g %in%  S100){print('S100')}
}

### 画饼图
vec_desc <- c('Oxidoreductase activity',
              'S100 protein binding',
              'Extracellular matrix organization',
              'Cellular response to zinc/copper ion',
              'Cartilage development')
group <- c('Carti','Ion','ECM','ROS','S100')

set <- c('ECM','ROS')
pct <- data.frame(group=set,pct=rep(100/length(set),length(set)),ncol = 1)
p<- ggplot(pct,aes(x="",y=pct,fill=group)) +
    geom_bar(stat = "identity",color="white",size =0.1) + 
    scale_fill_manual(values = col[set]) +
    coord_polar(theta = "y") +
    theme(axis.text.x = element_blank(),axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),panel.background =  element_blank())+guides(fill=FALSE)
p
ggsave(paste(paste(set,collapse = '_'),'.pdf',sep=''),p,width = 3,height = 3,units = 'cm')


#####################################################################################
#####################################################################################

# up
file_CSC_up <- paste0(path.TF.net, 'igraph_CSC_up.RDS')
igraph_up_2 <- readRDS(file_CSC_up)

sel.BP_SCS_up <- c('cell cycle arrest', 
                   'negative regulation of cell growth',
                   'negative regulation of stem cell differentiation',
                   'response to oxidative stress',
                   'p38MAPK cascade',
                   'I-kappaB kinase/NF-kappaB signaling',
                   'intrinsic apoptotic signaling pathway',
                   'extrinsic apoptotic signaling pathway',
                   'response to unfolded protein',
                   'regulation of RNA stability',
                   'activation of innate immune response',
                   'cellular response to tumor necrosis factor',
                   'regulation of inflammatory response', 
                   'cellular response to interleukin-1')

BP_SCS_up <- list.go.BP[['CSC_Microtia_increase']]
rownames(BP_SCS_up) <- BP_SCS_up$Description

df_GO_pre <- BP_SCS_up[sel.BP_SCS_up, c('Description', 'geneID')]
vec_desc <- c('Reduction of stem cell ability',
              'Response to oxidative stress',
              'NF-kappaB signaling and p38MAPK cascade',
              'Apoptotic signaling pathway',
              'Stability of protein and RNA',
              'Immune and inflammatory response')
vec_genes <- c(paste(df_GO_pre[1:3, 'geneID'], collapse = '/'),
               paste(df_GO_pre[4, 'geneID'], collapse = '/'),
               paste(df_GO_pre[5:6, 'geneID'], collapse = '/'),
               paste(df_GO_pre[7:8, 'geneID'], collapse = '/'),
               paste(df_GO_pre[9:10, 'geneID'], collapse = '/'),
               paste(df_GO_pre[11:14, 'geneID'], collapse = '/'))
df_GO <- data.frame(Description = vec_desc, geneID = vec_genes)

Stem <- unlist(strsplit(df_GO[df_GO$Description == 'Reduction of stem cell ability', 'geneID'],  '/'))
ROS <- unlist(strsplit(df_GO[df_GO$Description == 'Response to oxidative stress', 'geneID'],  '/'))
NFK <- unlist(strsplit(df_GO[df_GO$Description == 'NF-kappaB signaling and p38MAPK cascade', 'geneID'],  '/'))
Apop <- unlist(strsplit(df_GO[df_GO$Description == 'Apoptotic signaling pathway', 'geneID'],  '/'))
Stab <- unlist(strsplit(df_GO[df_GO$Description == 'Stability of protein and RNA', 'geneID'],  '/'))
IL <- unlist(strsplit(df_GO[df_GO$Description == 'Immune and inflammatory response', 'geneID'],  '/'))


library(RColorBrewer)
mycolor=brewer.pal(10,"Set3")
col=mycolor[c(1,3:7)]
group <- c('Apop','IL','NFK','Stem','ROS','Stab')
names(col) <- group

node_type <- V(igraph_up_2)$name
node_group <- rep('8', length(node_type))
node_group[node_type %in% c(Apop)] <- '1'
node_group[node_type %in% c(IL)] <- '2'
node_group[node_type %in% c(NFK)] <- '3'
node_group[node_type %in% c(Stem)] <- '4'
node_group[node_type %in% c(ROS)] <- '5'
node_group[node_type %in% c(Stab)] <- '6'
node_group[node_type %in% c('EGR1', 'KLF10', 'JUNB', 'REL', 'EGR3','ATF3', 'HIVEP3', 'IRX2', 'EGR2', 'KLF2', 'BCL3',
                               'CEBPB', 'CEBPD', 'STAT5A')] <- '7'
V(igraph_up_2)$group <- node_group

# gene expression FC
FCs <- V(igraph_up_2)$FC
names(FCs) <- V(igraph_up_2)$name
FC_TF <- FCs[names(FCs) %in% c('EGR1', 'KLF10', 'JUNB', 'REL', 'EGR3','ATF3', 'HIVEP3', 'IRX2', 'EGR2', 'KLF2', 'BCL3',
                               'CEBPB', 'CEBPD', 'STAT5A')]
FC_alpha <- (scale(FC_TF) + 1.5) / 3
FCs[rownames(FC_alpha)] <- FC_alpha
FCs[!(names(FCs) %in% c('EGR1', 'KLF10', 'JUNB', 'REL', 'EGR3','ATF3', 'HIVEP3', 'IRX2', 'EGR2', 'KLF2', 'BCL3',
                        'CEBPB', 'CEBPD', 'STAT5A'))] <- 1
V(igraph_up_2)$FC <- as.numeric(FCs)

ggraph_CSC_up <- as_tbl_graph(igraph_up_2)

plot_CSC_up <-
    ggraph(ggraph_CSC_up, layout = 'stress') +
    geom_edge_link(aes(edge_width=weight, alpha = weight),color="gray",
                   arrow = arrow(length = unit(2, 'mm')),
                   end_cap = circle(1, 'mm'), 
                   start_cap = circle(0.3, 'mm')) +
    scale_edge_width(range=c(0.5,1)) +
    scale_edge_alpha(range=c(0.2,1)) + 
    scale_size_continuous(range = c(2,8)) +
    geom_node_point(aes(size = page_rank, fill = group, alpha = FC),shape=21, color = 'transparent') +
    scale_color_manual(values = c(col,'firebrick3', 'gray')) + 
    scale_fill_manual(values = c(col,'firebrick3', 'gray')) + 
    geom_node_text(aes(filter = (group %in% c(1, 2, 3, 4, 5, 6)),label=name),size=2, repel = T) + 
    geom_node_text(aes(filter = group == 7,label=name),size=3) +
    theme_void() + theme(legend.position = 'none')
ggsave(filename = '/home/yzj/JingMA_NEW/res/compMicrotia/MicrotiavsNormal_inChildren/FIG/Fig5B_TF_CSC_up.pdf',
       plot_CSC_up,height = 10, width = 10, units = 'cm')


# color bar
df_plot <- data.frame(FC = scale(FC_TF), TF = rownames(scale(FC_TF)))
df_plot <- rbind(df_plot, data.frame(FC = -1.5, TF = 'DOWN'))
df_plot <- rbind(df_plot, data.frame(FC = 1.5, TF = 'UP'))
df_plot$NUL <- rep('1', nrow(df_plot))

plot_bar <- 
    ggplot(data = df_plot, aes(x = TF, y = NUL, fill = FC)) + 
    geom_tile() + 
    scale_fill_gradient(low = 'transparent', high = 'firebrick3', breaks = c(-1.5, 0, 1.5)) + 
    labs(fill = expression(paste("Scaled FC"['TF']))) + 
    theme(legend.title = element_text(size = 6, color = "black"),
          legend.text = element_text(size = 6, color = "black")) 
ggsave(plot = plot_CSC_down, path = path.TF.net,
       filename = 'TF_CSC_down_bar.png',
       height = 5, width = 5, units = 'cm')


## 查看gene所属类别
show.genes <- intersect(union(union(union(union(union(Apop,IL),NFK),ROS),Stab),Stem),node_type)
for(i in 1:length(show.genes)){
    g=show.genes[i]
    print(paste('!!!Gene: ',g))
    if(g %in% Apop){print('Apop')}
    if(g %in%  IL){print('IL')}
    if(g %in%  NFK){print('NFK')}
    if(g %in%  ROS){print('ROS')}
    if(g %in%  Stab){print('Stab')}
    if(g %in%  Stem){print('Stem')}
}

### 画饼图
vec_desc <- c('Reduction of stem cell ability',
              'Response to oxidative stress',
              'NF-kappaB signaling and p38MAPK cascade',
              'Apoptotic signaling pathway',
              'Stability of protein and RNA',
              'Immune and inflammatory response')
group <- c('Apop','IL','NFK','Stem','ROS','Stab')

set <- c('ECM','ROS')
pct <- data.frame(group=set,pct=rep(100/length(set),length(set)),ncol = 1)
p<- ggplot(pct,aes(x="",y=pct,fill=group)) +
    geom_bar(stat = "identity",color="white",size =0.1) + 
    scale_fill_manual(values = col[set]) +
    coord_polar(theta = "y") +
    theme(axis.text.x = element_blank(),axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),panel.background =  element_blank())+guides(fill=FALSE)
p
ggsave(paste(paste(set,collapse = '_'),'.pdf',sep=''),p,width = 3,height = 3,units = 'cm')




