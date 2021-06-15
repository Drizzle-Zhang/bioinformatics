#########
## Fig2主图: 比较儿童和成年人年龄段(C4/C6 vs C1/C2/C3/C5)
#########
library(tibble)
library(dplyr)
library(Seurat)
library(pheatmap)
library(ggplot2)
library(ggrepel)
#.libPaths(c("/home/yzj/R/x86_64-pc-linux-gnu-library/4.0","/home/zy/tools/R-4.0.0/library"))

########################
#### step2.1 DEG画heatmap,参照卵巢衰老的文章
########################
save_pheatmap_pdf <- function(x, filename, width=4, height=5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

pbmc_chond <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Chond.Rdata')
Idents(pbmc_chond) <- pbmc_chond$type
pbmc_C <- subset(pbmc_chond,idents = 'Normal')
pbmc_C@meta.data$Phase <- 'Adults'
pbmc_C$Phase[pbmc_C$batch %in% c('C4','C6')] <- 'Children'
pbmc_C$Phase <- factor(pbmc_C$Phase,levels = c('Children','Adults'))

MK.lst <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_inChildrenAdults.RDS')
print(names(MK.lst))

data <- MK.lst[['CSC']]
up_CSC <- rownames(data)[data$avg_logFC > log(1.5) & data$p_val_adj < 0.05]
dn_CSC <- rownames(data)[data$avg_logFC < -log(1.5) & data$p_val_adj < 0.05]

data <- MK.lst[["TC"]]
up_TC <- rownames(data)[data$avg_logFC > log(1.5) & data$p_val_adj < 0.05]
dn_TC <- rownames(data)[data$avg_logFC < -log(1.5) & data$p_val_adj < 0.05]

data <- MK.lst[["C1"]]
up_C1 <- rownames(data)[data$avg_logFC > log(1.5) & data$p_val_adj < 0.05]
dn_C1 <- rownames(data)[data$avg_logFC < -log(1.5) & data$p_val_adj < 0.05]

data <- MK.lst[["C2"]]
up_C2 <- rownames(data)[data$avg_logFC > log(1.5) & data$p_val_adj < 0.05]
dn_C2 <- rownames(data)[data$avg_logFC < -log(1.5) & data$p_val_adj < 0.05]

get_values <- function(sigCSC,sigTC,sigC1,sigC2){
  sigGene <- unique(c(sigCSC,sigTC,sigC1,sigC2))
  values <- matrix(c(rep(0,4*length(sigGene))),ncol = 4,dimnames = list(sigGene,c('CSC','TC','C1','C2')))
  for(i in 1:length(sigGene)){
    g=sigGene[i]
    if(g %in% sigCSC){values[i,1] <-1};
    if(g %in% sigTC){values[i,2] <-1};
    if(g %in% sigC1){values[i,3] <-1};
    if(g %in% sigC2){values[i,4] <-1};
  }
  values_sum <- apply(values, 1, sum)
  values <- values[order(values_sum,decreasing = T),]
  return(values)
}

## 对成人来说，下调矩阵
upValues_mtx <- get_values(up_CSC,up_TC,up_C1,up_C2)
up_sum <- apply(upValues_mtx,1,sum)
up_df <- upValues_mtx[-(which(up_sum>1)),]

annotation_col = data.frame(CellType = factor(c("CSC", "TC","C1","C2")))
rownames(annotation_col) <- colnames(upValues_mtx)

annotation_row = data.frame(GeneClass = factor(rep(c("Common", "CSC", "TC","C1","C2"), 
                                                   c(length(which(up_sum>1)), 
                                                     length(which(up_df[,1]==1)), 
                                                     length(which(up_df[,2]==1)),
                                                     length(which(up_df[,3]==1)),
                                                     length(which(up_df[,4]==1))))))
rownames(annotation_row) = rownames(upValues_mtx)

ann_colors = list( CellType = c(CSC="#EE9572",TC="#B2DF8A",C1="#A6CEE3",C2="#9999FF"),
                   GeneClass = c(Common='grey',CSC="#EE9572",TC="#B2DF8A",C1="#A6CEE3",C2="#9999FF"))

p_UP <- pheatmap(upValues_mtx,cluster_rows = F,cluster_cols = F,color =  colorRampPalette(c("#EFEFEF", "white", "#7F99CE"))(100),
                 border_color ='transparent',show_rownames = F,angle_col='45',
                 annotation_row = annotation_row,annotation_colors = ann_colors,legend=F,annotation_legend = FALSE)
save_pheatmap_pdf(p_UP,'JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DEGHeatmap_UP.pdf',height = 4,width = 2)
saveRDS(upValues_mtx,'JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DEGHeatmap_UPmtx.RDS')

## 对成人来说上调矩阵
dnValues_mtx <- get_values(dn_CSC,dn_TC,dn_C1,dn_C2)
dn_sum <- apply(dnValues_mtx,1,sum)
dn_df <- dnValues_mtx[-(which(dn_sum>1)),]

annotation_col = data.frame(CellType = factor(c("CSC", "TC","C1","C2")))
rownames(annotation_col) <- colnames(dnValues_mtx)

annotation_row = data.frame(GeneClass = factor(rep(c("Common", "CSC", "TC","C1","C2"), 
                                                   c(length(which(dn_sum>1)), 
                                                     length(which(dn_df[,1]==1)), 
                                                     length(which(dn_df[,2]==1)),
                                                     length(which(dn_df[,3]==1)),
                                                     length(which(dn_df[,4]==1))))))
rownames(annotation_row) = rownames(dnValues_mtx)

ann_colors = list( CellType = c(CSC="#EE9572",TC="#B2DF8A",C1="#A6CEE3",C2="#9999FF"),
                   GeneClass = c(Common='grey',CSC="#EE9572",TC="#B2DF8A",C1="#A6CEE3",C2="#9999FF"))

p_DN <- pheatmap(dnValues_mtx,cluster_rows = F,cluster_cols = F,color =  colorRampPalette(c("#EFEFEF", "white","#B15E72"))(100),
                 border_color ='transparent',show_rownames = F,legend=F,angle_col='45',
                 annotation_row = annotation_row,annotation_colors = ann_colors,annotation_legend = FALSE)
save_pheatmap_pdf(p_DN,'JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DEGHeatmap_DN.pdf',height = 4,width = 2)
saveRDS(dnValues_mtx,'JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DEGHeatmap_DNmtx.RDS')


up_sum <- apply(upValues_mtx,1,sum)
length(which(up_sum==4))
length(which(up_sum>1))
up_df <- upValues_mtx[-(which(up_sum>1)),]
length(which(up_df[,1]==1))
length(which(up_df[,2]==1))
length(which(up_df[,3]==1))
length(which(up_df[,4]==1))

dn_sum <- apply(dnValues_mtx,1,sum)
length(which(dn_sum==4))
length(which(dn_sum>1))
dn_df <- dnValues_mtx[-(which(dn_sum>1)),]
length(which(dn_df[,1]==1))
length(which(dn_df[,2]==1))
length(which(dn_df[,3]==1))
length(which(dn_df[,4]==1))



#########
### 2.2 挑选term组合画图
#########
term_down_BP <- c("cartilage development","chondrocyte differentiation","cartilage condensation","chondrocyte morphogenesis",
                  "extracellular matrix organization","collagen fibril organization",
                  "NAD metabolic process")
term_down_MF <- c("cartilage development","chondrocyte differentiation","cartilage condensation","chondrocyte morphogenesis",
                  "extracellular matrix organization","collagen fibril organization",
                  "NAD metabolic process")


library(xlsx)
CSC_DN <- read.xlsx('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro/FC1.5_adjP0.05/CSC_BP.xlsx',sheetName = 'DN')
CSC_DN <- CSC_DN[CSC_DN$p.adjust < 0.1,]
print(CSC_DN$Description)
index <- c(8,17,22,24,29,49)
pickCSC_DN <- CSC_DN[index,]
geneCSC_DN <- c()
for(i in 1:nrow(pickCSC_DN)){
  geneCSC_DN <- c(geneCSC_DN,unlist(strsplit(pickCSC_DN[i,8],'/')))
}
geneCSC_DN <-unique(geneCSC_DN)
print(length(geneCSC_DN))


CSC_UP <- read.xlsx('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro/FC1.5_adjP0.05/CSC_BP.xlsx',sheetName = 'UP')
CSC_UP <- CSC_UP[CSC_UP$p.adjust < 0.1,]
print(CSC_UP$Description)
index <- c(12,14,27,30,35,89,123,145,153,154,193,217,302)
pickCSC_UP <- CSC_UP[index,]
geneCSC_UP <- c()
for(i in 1:nrow(pickCSC_UP)){
  geneCSC_UP <- c(geneCSC_UP,unlist(strsplit(pickCSC_UP[i,8],'/')))
}
geneCSC_UP <-unique(geneCSC_UP)
print(length(geneCSC_UP))


####
Trans_DN <- read.xlsx('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro/FC1.5_adjP0.05/TC_BP.xlsx',sheetName = 'DN')
Trans_DN <- Trans_DN[Trans_DN$p.adjust < 0.1,]
print(Trans_DN$Description)
index <- c(11,13,18,27,29)
pickTrans_DN <- Trans_DN[index,]
geneTrans_DN <- c()
for(i in 1:nrow(pickTrans_DN)){
  geneTrans_DN <- c(geneTrans_DN,unlist(strsplit(pickTrans_DN[i,8],'/')))
}
geneTrans_DN <-unique(geneTrans_DN)
print(length(geneTrans_DN))


Trans_UP <- read.xlsx('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro/FC1.5_adjP0.05/TC_BP.xlsx',sheetName = 'UP')
Trans_UP <- Trans_UP[Trans_UP$p.adjust < 0.1,]
print(Trans_UP$Description)
index <- c(12,28,73,118,123,130)
pickTrans_UP <- Trans_UP[index,]
geneTrans_UP <- c()
for(i in 1:nrow(pickTrans_UP)){
  geneTrans_UP <- c(geneTrans_UP,unlist(strsplit(pickTrans_UP[i,8],'/')))
}
geneTrans_UP <-unique(geneTrans_UP)
print(length(geneTrans_UP))


####
Chond1_DN <- read.xlsx('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro/FC1.5_adjP0.05/C1_BP.xlsx',sheetName = 'DN')
Chond1_DN <- Chond1_DN[Chond1_DN$p.adjust < 0.1,]
print(Chond1_DN$Description)
index <- c(40)
pickChond1_DN <- Chond1_DN[index,]
geneChond1_DN<- c()
for(i in 1:nrow(pickChond1_DN)){
  geneChond1_DN <- c(geneChond1_DN,unlist(strsplit(pickChond1_DN[i,8],'/')))
}
geneChond1_DN <-unique(geneChond1_DN)
print(length(geneChond1_DN))


Chond1_UP <- read.xlsx('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro/FC1.5_adjP0.05/C1_BP.xlsx',sheetName = 'UP')
Chond1_UP <- Chond1_UP[Chond1_UP$p.adjust < 0.1,]
print(Chond1_UP$Description)
index <- c(12,48,90,138,143,165,200)
pickChond1_UP <- Chond1_UP[index,]
geneChond1_UP <- c()
for(i in 1:nrow(pickChond1_UP)){
  geneChond1_UP <- c(geneChond1_UP,unlist(strsplit(pickChond1_UP[i,8],'/')))
}
geneChond1_UP <-unique(geneChond1_UP)
print(length(geneChond1_UP))


####
Chond2_DN <- read.xlsx('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro/FC1.5_adjP0.05/C2_BP.xlsx',sheetName = 'DN')
Chond2_DN <- Chond2_DN[Chond2_DN$p.adjust < 0.1,]
print(Chond2_DN$Description)
index <- c(63)
pickChond2_DN <- Chond2_DN[index,]
geneChond2_DN<- c()
for(i in 1:nrow(pickChond2_DN)){
  geneChond2_DN <- c(geneChond2_DN,unlist(strsplit(pickChond2_DN[i,8],'/')))
}
geneChond2_DN <-unique(geneChond2_DN)
print(length(geneChond2_DN))


Chond2_UP <- read.xlsx('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro/FC1.5_adjP0.05/C2_BP.xlsx',sheetName = 'UP')
Chond2_UP <- Chond2_UP[Chond2_UP$p.adjust < 0.1,]
print(Chond2_UP$Description)
index <- c(13,14,43,51,74,144,166,199,204)
pickChond2_UP <- Chond2_UP[index,]
geneChond2_UP <- c()
for(i in 1:nrow(pickChond2_UP)){
  geneChond2_UP <- c(geneChond2_UP,unlist(strsplit(pickChond2_UP[i,8],'/')))
}
geneChond2_UP <-unique(geneChond2_UP)
print(length(geneChond2_UP))

bar_DN <- rbind(pickCSC_DN,pickTrans_DN,pickChond1_DN,pickChond2_DN)
bar_DN$CellType <- c(rep('CSC',nrow(pickCSC_DN)),rep('TC',nrow(pickTrans_DN)),
                     rep('C1',nrow(pickChond1_DN)),rep('C2',nrow(pickChond2_DN)))
bar_DN$CellType <- factor(bar_DN$CellType,levels = c('CSC','TC','C1','C2'))
bar_DN$Group <- 'Children'
bar_DN$log10Pval <- -log(bar_DN$p.adjust,10)

bar_UP <- rbind(pickCSC_UP,pickTrans_UP,pickChond1_UP,pickChond2_UP)
bar_UP$CellType <- c(rep('CSC',nrow(pickCSC_UP)),rep('TC',nrow(pickTrans_UP)),
                     rep('C1',nrow(pickChond1_UP)),rep('C2',nrow(pickChond2_UP)))
bar_UP$CellType <- factor(bar_UP$CellType,levels = c('CSC','TC','C1','C2'))
bar_UP$Group <- 'Adults'
bar_UP$log10Pval <- log(bar_UP$p.adjust,10)
bar_df <- rbind(bar_DN,bar_UP)

levels_DN=rev(c("cartilage development","chondrocyte differentiation","cartilage condensation","chondrocyte morphogenesis",
                "extracellular matrix organization","collagen fibril organization",
                "NAD metabolic process"))
setdiff(unique(bar_DN$Description),levels_DN)

levels_UP=rev(c("negative regulation of stem cell differentiation","cell cycle arrest",
                "autophagy","aging","cellular senescence",
                "response to oxidative stress","reactive oxygen species metabolic process","reactive oxygen species biosynthetic process",
                "cell death in response to oxidative stress",
                "DNA damage response, signal transduction by p53 class mediator",
                "ERK1 and ERK2 cascade", 'p38MAPK cascade',
                # "positive regulation of p38MAPK cascade",
                # "response to interleukin-6","extrinsic apoptotic signaling pathway",
                "intrinsic apoptotic signaling pathway"))
setdiff(unique(bar_UP$Description),levels_UP)

bar_df$Description <- factor(bar_df$Description,levels = rev(c(levels_UP,levels_DN)))
bar_df$Count <- as.numeric(bar_df$Count)


library(reshape2)
mat.plot <- bar_df[,c('Description','CellType','Group','log10Pval')]
mat.plot <- dcast(mat.plot,Description~CellType+Group)
mat.plot[is.na(mat.plot)] <- 0
rownames(mat.plot) <- mat.plot$Description
mat.plot <- mat.plot[,-1]
colNames <- c('CSC_Children','TC_Children','C1_Children','C2_Children','CSC_Adults','TC_Adults','C1_Adults','C2_Adults')
mat.plot <- dplyr::select(mat.plot,colNames)

# col annotation
annotation_col = data.frame(
  CellType = factor(c(rep('CSC', 2),rep('TC', 2),rep('C1', 2),rep('C2', 2)), 
                    levels = c('CSC', 'TC','C1', 'C2')), 
  Phase = factor(rep(c('Children', 'Adults'), 4), levels = c('Children', 'Adults')),
  row.names = colNames
)

annotation_col = data.frame(
  CellType = factor(rep(c('CSC', 'TC','C1', 'C2'), 2), 
                    levels = c('CSC', 'TC','C1', 'C2')), 
  Phase = factor(rep(c('Children', 'Adults'), each=4), levels = c('Children', 'Adults')),
  row.names = colNames
)

ann_colors = list(
  CellType = c(CSC="#EE9572",TC="#B2DF8A",C1="#A6CEE3",C2="#9999FF"),
  Phase = c(Children = "#6C6C6C", Adults = "#637FBF")
) 

bk <- c(seq(-8,-0.1,by=0.01),seq(0,8,by=0.01))
plot.heatmap <- pheatmap::pheatmap(mat.plot,
                   cluster_rows = F, cluster_cols = F, scale = "none",
                   display_numbers = F,
                   annotation_col = annotation_col ,annotation_colors = ann_colors,
                   show_rownames = T, show_colnames = F, legend = T, 
                   gaps_col = c(4),
                   color = c(colorRampPalette(colors = c("red","white"))(length(bk)/2),colorRampPalette(colors = c("white","blue"))(length(bk)/2)),
                   legend_breaks=seq(-8,8,2),
                   breaks=bk
)

ggsave('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/Fig3C_pickHeatmap.pdf',plot.heatmap,width = 8,height = 5)


########################
#### 2.2 挑选基因画vlnplot
########################
library(ggpubr)
get_vlnplot <- function(gene){
  pickEXP <- data.frame(cells=colnames(EXP),exp=as.numeric(EXP[rownames(EXP) ==gene,]),celltype=pbmc_C$celltype,phase=pbmc_C$Phase)
  p <- ggplot(pickEXP, aes(x=celltype, y=exp,fill=phase)) + 
    geom_violin(trim=FALSE,color="white") + 
    geom_boxplot(width=0.2,position=position_dodge(0.9))+
    scale_fill_manual(values = c("#6C6C6C", "#637FBF"))+ 
    theme_bw()+ 
    labs(title=gene)+
    theme(axis.text.x=element_blank(),axis.ticks.x =element_blank(),
          axis.text.y=element_text(size=12,colour="black"),axis.title.y=element_text(size = 12,colour="black"), axis.ticks.y =element_line(colour="black"),
          legend.text=element_text(colour="black", size=12),legend.title=element_blank(),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5,size = 12,face = 'bold.italic'))+ 
    ylab("")+xlab("")+ 
    facet_wrap(~celltype,ncol = 4,scales= "free_x")+
    theme(strip.background = element_rect(color = "black", fill = "#8C90C6",size = 1.2),
          strip.text.x = element_text(size = 10, color = "black",face = 'bold'),
          panel.grid = element_blank(),panel.border = element_rect(color = 'black',size = 2))+
    stat_compare_means(label = "p.signif",label.x=1.5)
  return(p)
}

EXP <- as.data.frame(pbmc_C@assays$RNA@data)

pick_genes <- c('COL2A1','COL11A1','COL11A2','COL9A1','COL11A2','COL9A2','COL9A2','ELN','TIMP4',
                'MATN3','VIT','CYTL1','PTX3','PTGS2','GPX3','SOD2','MGP','MMP3','CDKN1A','IL6')
pdf('JingMA_NEW/res/compControl/ChildrenvsAdults/pickGene_vlnplot.pdf',width = 6,height = 3)
for(i in 1:length(pick_genes)){
  gene=pick_genes[i]
  p <- get_vlnplot(gene)
  print(p)
}
dev.off()
