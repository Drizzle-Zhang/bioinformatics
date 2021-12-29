library(Seurat)
library(ggplot2)
library(dplyr)

pbmc <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')

subpbmc <- subset(pbmc,cells = colnames(pbmc)[pbmc$celltype.abbr %in% c('CSC','TC','Chondrocyte1','Chondrocyte2')
                                              & pbmc$batch %in% c('C4','C6','M1','M2','M3')])
subpbmc$celltype.abbr <- factor(subpbmc$celltype.abbr,levels = rev(c('CSC','TC','Chondrocyte1','Chondrocyte2')))

### A. Proportion.
phylog_df <- subpbmc@meta.data[,c('type',"celltype.abbr")]
phylog_df <- table(phylog_df$type,phylog_df[,"celltype.abbr"])
phylog_df <- data.frame(phylog_df)
colnames(phylog_df) <- c('Status','CellType','Freq')
phylog_df$CellType <- factor(phylog_df$CellType)

Color <- c("#33A02C","#B2DF8A" ,"#1F78B4","#A6CEE3")
p1 <- ggplot(phylog_df,aes(x=Status,y=Freq,fill=CellType))+
  geom_col(position = "fill", width = 0.9)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(size = 10,colour = "black"),
        axis.line = element_line(size=0.7, colour = "black"),
        axis.title.y = element_text(size=10))+
  labs(x='',y='')+theme(legend.position="right")+
  scale_fill_manual(values = rev(Color))
p1

ggsave('/home/yzj/JingMA_NEW/res/Harmony/ALL/FIG/Barplot_PropChond_Microtia.pdf',p1,width = 4,height = 4)


### B. Heatmap of GO
save_pheatmap_pdf <- function(x, filename, width=4, height=5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

fc.cutoff <- 0.5
path.M123 <- '/home/disk/drizzle/wgk/microtia_child_M1M2M3/'
path.cutoff <- paste0(path.M123, 'cutoff_', fc.cutoff, '/')
# path.M12 <- '/home/disk/drizzle/wgk/microtia_child_M1M2/'
# path.cutoff <- paste0(path.M12, 'cutoff_', fc.cutoff, '/')
file.marker.go <- paste0(path.cutoff, 'marker_go.Rdata')
list.marker.go <- readRDS(file.marker.go)

file.go.BP <- paste0(path.cutoff, 'GO_BP.Rdata')
list.go.BP <- readRDS(file.go.BP)
file.go.MF <- paste0(path.cutoff, 'GO_MF.Rdata')
list.go.MF <- readRDS(file.go.MF)


# select GO
df.GO <- data.frame(stringsAsFactors = F)
# Chondral stem cell
GO.BP.CSC.M <- c('response to oxidative stress', 'response to unfolded protein', 
                 'response to tumor necrosis factor', 'RNA splicing',
                 'RNA localization', 'positive regulation of defense response')
GO.BP.CSC.N <- c('ribosome biogenesis', 'response to copper ion', 'oxidative phosphorylation',
                 'extracellular matrix organization', 'skeletal system development',
                 'cell aggregation', 'cellular zinc ion homeostasis')
sel.GO.BP <- c('response to oxidative stress', 
               'response to unfolded protein', 
               'positive regulation of defense response', 
               # 'regulation of inflammatory response',
               'p38MAPK cascade', 'ERK1 and ERK2 cascade',
               # 'regulation of ERK1 and ERK2 cascade', 
               'intrinsic apoptotic signaling pathway', 
               'cell cycle arrest',
               'negative regulation of stem cell differentiation',
               'negative regulation of cell growth',
               'RNA splicing', 'RNA localization', 
               'vascular endothelial growth factor production', 'angiogenesis', 
               'positive regulation of vasculature development',
               'positive regulation of cell migration',
               'negative regulation of cell adhesion',
               'translational initiation', 'ribosome biogenesis', 
               'cartilage condensation', 
               'extracellular matrix organization', 
               # 'skeletal system development',
               'connective tissue development',
               'zinc ion homeostasis')
sel.GO.MF <- c('extracellular matrix structural constituent', 
               'extracellular matrix binding', 'S100 protein binding')

terms <- c("Chondral stem cell_Microtia_increase",
           "Chondral stem cell_Microtia_decrease",
           "Transitional chondrocyte_Microtia_increase",
           "Transitional chondrocyte_Microtia_decrease",
           "Chondrocyte1_Microtia_increase",
           "Chondrocyte1_Microtia_decrease",
           "Chondrocyte2_Microtia_increase",
           "Chondrocyte2_Microtia_decrease")

df.plot <- data.frame(stringsAsFactors = F)
for (term in terms) {
  cell <- strsplit(term, split = '_')[[1]][1]
  status <- strsplit(term, split = '_')[[1]][3]
  sub.BP <- list.go.BP[[term]]
  rownames(sub.BP) <- sub.BP$Description
  sub.BP <- sub.BP[sub.BP$p.adjust < 0.1,]
  sel.BP <- sub.BP[sel.GO.BP, c('Description', 'pvalue', 'geneID')]
  sel.BP$Description <- sel.GO.BP
  sub.MF <- list.go.MF[[term]]
  rownames(sub.MF) <- sub.MF$Description
  sub.MF <- sub.MF[sub.MF$p.adjust < 0.1,]
  sel.MF <- sub.MF[sel.GO.MF, c('Description', 'pvalue', 'geneID')]
  sel.MF$Description <- sel.GO.MF
  sub.plot <- rbind(sel.BP, sel.MF)
  sub.plot$pvalue[is.na(sub.plot$pvalue)] <- 1
  sub.plot$CellType <- rep(cell, nrow(sub.plot))
  if (status == 'increase') {
    sub.plot$Status <- rep('Microtia', nrow(sub.plot))
    sub.plot$Coeff <- rep(1, nrow(sub.plot))
  } else {
    sub.plot$Status <- rep('Normal', nrow(sub.plot))
    sub.plot$Coeff <- rep(-1, nrow(sub.plot))
  }
  df.plot <- rbind(df.plot, sub.plot)
}
df.plot$log10Pval <- -log10(df.plot$pvalue)
df.plot$log10Pval[abs(df.plot$log10Pval) > 10] = 10
df.plot$log10Pval <- df.plot$log10Pval * df.plot$Coeff
df.plot$abs_log10Pval <- abs(df.plot$log10Pval)
df.plot$Description <- factor(df.plot$Description, levels = c(sel.GO.BP, sel.GO.MF))
df.plot$CellType_raw <- df.plot$CellType
df.plot$CellType[df.plot$CellType_raw=='Chondral stem cell']  <- 'CSC'
df.plot$CellType[df.plot$CellType_raw=='Transitional chondrocyte']  <- 'TC'
df.plot$CellType <- factor(df.plot$CellType,levels = c('CSC','TC','Chondrocyte1','Chondrocyte2'))


df.plot$col_name <- paste(df.plot$CellType, df.plot$Status, sep = '_')
mat.plot <- reshape2::dcast(df.plot, Description ~ col_name, value.var = 'log10Pval')
row.names(mat.plot) <- mat.plot$Description
mat.plot$Description <- NULL

ann_colors = list(
  CellType = c(CSC="#33A02C",TC="#B2DF8A",Chondrocyte1="#1F78B4",Chondrocyte2="#A6CEE3"),
  Status = c(Normal = "#637FBF", Microtia = "#6C6C6C")
) 

# col annotation
annotation_col = data.frame(
  CellType = factor(c(rep('CSC', 2),
                      rep('Chondrocyte1', 2),
                      rep('Chondrocyte2', 2),
                      rep('TC', 2)), 
                    levels = c('CSC', 'TC',
                               'Chondrocyte1', 'Chondrocyte2')), 
  Status = factor(rep(c('Microtia', 'Normal'), 4), levels = c('Normal', 'Microtia')),
  row.names = colnames(mat.plot)
)

cols <- c("CSC_Microtia", "TC_Microtia", 
          "Chondrocyte1_Microtia", "Chondrocyte2_Microtia",
          "CSC_Normal", "TC_Normal",
          "Chondrocyte1_Normal", "Chondrocyte2_Normal")
mat.plot <- mat.plot[rev(rownames(mat.plot)), cols]
annotation_col <- annotation_col[cols,]

ann_colors = list(
  CellType = c(CSC="#33A02C",TC="#B2DF8A",Chondrocyte1="#1F78B4",Chondrocyte2="#A6CEE3"),
  Status = c(Normal = "#637FBF", Microtia = "#6C6C6C")
) 

p2 <- pheatmap::pheatmap(mat.plot,
                   color = colorRampPalette(c('blue', 'white', 'red'))(100),
                   cluster_rows = F, cluster_cols = F, scale = "none",
                   display_numbers = F,
                   annotation_col = annotation_col ,annotation_colors = ann_colors,
                   show_rownames = T, show_colnames = F, legend = T, 
                   # fontsize_row = 18, fontsize_col = 15,
                   gaps_col = c(4)
)
save_pheatmap_pdf('JingMA_NEW/res/compMicrotia/Fig3B.pdf',plot.heatmap,width = 2,height = 4)

