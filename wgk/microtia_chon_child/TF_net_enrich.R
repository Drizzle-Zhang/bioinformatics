.libPaths(c('/home/zy/R/x86_64-pc-linux-gnu-library/4.0', 
            .libPaths(),
            '/home/yzj/R/x86_64-pc-linux-gnu-library/4.0'))
library(Seurat)
library(patchwork)
library(tidyverse)
library(SCENIC)
library(igraph)
library(ggraph)
library(tidygraph)
library(xlsx)
data(c2BroadSets)
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(AnnotationHub)
library(org.Hs.eg.db)   #人类注释数据库
library(clusterProfiler)
library(DOSE)
library(dplyr)
library(tidyverse)
library(reshape2)

regulon2Targets <- readRDS('/home/yzj/JingMA_NEW/res/SCENIC_main/int/2.5_regulonTargetsInfo.Rds')
TF <- unique(regulon2Targets$TF)
path.TF.net <- '/home/disk/drizzle/wgk/TF_net/'

# down
file_CSC_down <- paste0(path.TF.net, 'igraph_CSC_down.RDS')
igraph_down_2 <- readRDS(file_CSC_down)
dnGene<- V(igraph_down_2)$name

path.down <- paste0(path.TF.net, 'TF_down/')
if (!file.exists(path.down)) {
    dir.create(path.down)
}


##
TFlst <- c('SOX8','SOX5','ARID5B','DBP')
tf='SOX8'
list.go.down <- list()
for(tf in TFlst) {
    sub.regulon2Targets <- regulon2Targets[regulon2Targets$TF %in% tf & regulon2Targets$gene %in% dnGene & regulon2Targets$highConfAnnot =='TRUE',]
    testgenes <- sub.regulon2Targets$gene
    
    
    symbol2id=mapIds(org.Hs.eg.db,testgenes,"ENTREZID",'SYMBOL')
    id=symbol2id[which(symbol2id!='')] #提取出非NA的ENTREZID
    
    if(length(id) > 0){
        #GO BP 富集分析#
        ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "BP", pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
        # ego <- clusterProfiler::simplify(ego,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")
        ego_res <- as.data.frame(ego)
        if(nrow(ego_res)>0){
            write.xlsx(ego_res,paste(path.down,tf,"_BP.xlsx",sep=''),row.names = FALSE,sheetName = type,append = TRUE)
        } 
    }
    print(ego_res$Description)
    list.go.down[[tf]] <- ego_res
}

list.sel.GO <- list()
list.sel.GO$SOX8 <- c('extracellular matrix organization', 
                      'cartilage development',
                      'chondrocyte differentiation')
list.sel.GO$SOX5 <- c('extracellular matrix organization', 
                      'cartilage development',
                      'chondrocyte differentiation')
list.sel.GO$ARID5B <- c('extracellular matrix organization', 
                        "response to copper ion", 
                        "keratan sulfate biosynthetic process", 
                        "cartilage development",
                        "response to zinc ion")
list.sel.GO$DBP <- c("response to transforming growth factor beta", 
                     "stem cell differentiation", 
                     "positive regulation of chondrocyte differentiation", 
                     "positive regulation of cartilage development",
                     "developmental induction")

## barplot
sort.TFs <- TFlst
colors <- c("#BC80BD", "#80B1D3", "#F4A460", "#FB8072")
df.plot <- data.frame()
GO_terms <- c()
i = 0
for (TF in sort.TFs) {
    i = i + 1
    sub.go <- list.go.down[[TF]]
    sel.go.term <- list.sel.GO[[TF]]
    sel.go <- sub.go[sub.go$Description %in% sel.go.term, 
                     c('Description', 'pvalue')]
    sel.go$log10Pval <- -log10(sel.go$pvalue)
    sel.go$TFgene <- rep(TF, nrow(sel.go))
    # sel.go$Description <- factor(sel.go$Description, levels = rev(sel.go.term))
    df.plot <- rbind(df.plot, sel.go)
    GO_terms <- c(GO_terms, setdiff(sel.go$Description, GO_terms))
}
col_name <- paste(df.plot$TFgene, df.plot$Description, sep = '_')
df.plot$col_name <- factor(col_name, levels = rev(col_name))
df.plot$TFgene <- factor(df.plot$TFgene, levels = sort.TFs)
df.plot$Description <- factor(df.plot$Description, levels = (GO_terms))

p <- ggplot(df.plot, aes(x = TFgene, y = Description, 
                         color = TFgene, size = log10Pval)) + 
    geom_point(fill = 'cornsilk') +
    scale_color_manual(breaks = sort.TFs,
                       values = colors) + 
    scale_size_continuous(range = c(2,5)) +
    labs(x = '', y = 'GO term', color = 'TF',
         size = expression(paste("-log"[10], "(adj", italic("P"), "-value)"))) +
    theme(panel.background = element_rect(color = 'black',
                                          fill = 'transparent'), 
          panel.grid.major = element_line(colour = 'gray', size = 0.2, linetype = 5),
          axis.title = element_text(size = 7, color = 'black'), 
          axis.text.y = element_text(size = 6, color = 'black'), 
          axis.text.x = element_text(size = 6, color = 'black'),
          legend.text = element_text(size = 6, color = 'black'),
          legend.title = element_text(size = 7, color = 'black'),
          legend.key = element_blank()) + 
    guides(colour = guide_legend(override.aes = list(size=2)))
ggsave(plot = p, path = path.down, 
       filename = paste0('GO_down.png'),
       height = 7.5, width = 12, units = 'cm')


# up
file_CSC_up <- paste0(path.TF.net, 'igraph_CSC_up.RDS')
igraph_up_2 <- readRDS(file_CSC_up)
dnGene<- V(igraph_up_2)$name

path.up <- paste0(path.TF.net, 'TF_up/')
if (!file.exists(path.up)) {
    dir.create(path.up)
}


##
TFlst <- c('EGR1','ATF3','EGR3','KLF10', 'REL', 'BCL3', 'JUNB', 'CEBPB', 'CEBPD')
# tf='SOX8'
list.go.up <- list()
for(tf in TFlst) {
    sub.regulon2Targets <- regulon2Targets[regulon2Targets$TF %in% tf & regulon2Targets$gene %in% dnGene & regulon2Targets$highConfAnnot =='TRUE',]
    testgenes <- sub.regulon2Targets$gene
    
    
    symbol2id=mapIds(org.Hs.eg.db,testgenes,"ENTREZID",'SYMBOL')
    id=symbol2id[which(symbol2id!='')] #提取出非NA的ENTREZID
    
    if(length(id) > 0){
        #GO BP 富集分析#
        ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "BP", pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
        # ego <- clusterProfiler::simplify(ego,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")
        ego_res <- as.data.frame(ego)
        if(nrow(ego_res)>0){
            write.xlsx(ego_res,paste(path.up,tf,"_BP.xlsx",sep=''),row.names = FALSE,sheetName = type,append = TRUE)
        } 
    }
    # print(ego_res$Description)
    list.go.up[[tf]] <- ego_res
}

list.sel.GO <- list()
list.sel.GO$EGR1 <- c('cellular response to oxidative stress', 
                      'I-kappaB kinase/NF-kappaB signaling',
                      'response to unfolded protein')
list.sel.GO$ATF3 <- c('response to unfolded protein', 
                      'negative regulation of cell growth',
                      'extrinsic apoptotic signaling pathway',
                      'intrinsic apoptotic signaling pathway')
list.sel.GO$EGR3 <- c('response to unfolded protein', 
                        "extrinsic apoptotic signaling pathway", 
                        "I-kappaB kinase/NF-kappaB signaling")
list.sel.GO$KLF10 <- c("extrinsic apoptotic signaling pathway", 
                     "intrinsic apoptotic signaling pathway")
list.sel.GO$JUNB <- c("intrinsic apoptotic signaling pathway",
                      'cellular response to oxidative stress',
                      'negative regulation of cell growth',
                      'p38MAPK cascade')
list.sel.GO$REL <- c("regulation of inflammatory response", 
                     "p38MAPK cascade", 
                     "intrinsic apoptotic signaling pathway",
                     "extrinsic apoptotic signaling pathway",
                     'cellular response to oxidative stress', 
                     'cellular response to tumor necrosis factor',
                     'I-kappaB kinase/NF-kappaB signaling', 
                     'cell cycle arrest',
                     'cellular response to interleukin-1')
list.sel.GO$BCL3 <- c("I-kappaB kinase/NF-kappaB signaling", 
                     "regulation of inflammatory response", 
                     "cellular response to tumor necrosis factor", 
                     "cellular response to interleukin-1",
                     "developmental induction")
list.sel.GO$CEBPB <- c("cellular senescence")
list.sel.GO$CEBPD <- c("regulation of inflammatory response")

## barplot
sort.TFs <- TFlst
mycolor=brewer.pal(10,"Set3")
colors <- mycolor[c(1,3:10)]
df.plot <- data.frame()
GO_terms <- c()
i = 0
for (TF in sort.TFs) {
    i = i + 1
    sub.go <- list.go.up[[TF]]
    sel.go.term <- list.sel.GO[[TF]]
    sel.go <- sub.go[sub.go$Description %in% sel.go.term, 
                     c('Description', 'pvalue')]
    sel.go$log10Pval <- -log10(sel.go$pvalue)
    sel.go$TFgene <- rep(TF, nrow(sel.go))
    # sel.go$Description <- factor(sel.go$Description, levels = rev(sel.go.term))
    df.plot <- rbind(df.plot, sel.go)
    GO_terms <- c(GO_terms, setdiff(sel.go$Description, GO_terms))
}
col_name <- paste(df.plot$TFgene, df.plot$Description, sep = '_')
df.plot$col_name <- factor(col_name, levels = rev(col_name))
df.plot$TFgene <- factor(df.plot$TFgene, levels = sort.TFs)
df.plot$Description <- factor(df.plot$Description, levels = (GO_terms))

p <- ggplot(df.plot, aes(x = TFgene, y = Description, 
                         color = TFgene, size = log10Pval)) + 
    geom_point(fill = 'cornsilk') +
    scale_color_manual(breaks = sort.TFs,
                       values = colors) + 
    scale_size_continuous(range = c(2,5)) +
    labs(x = '', y = 'GO term', color = 'TF',
         size = expression(paste("-log"[10], "(adj", italic("P"), "-value)"))) +
    theme(panel.background = element_rect(color = 'black',
                                          fill = 'transparent'), 
          panel.grid.major = element_line(colour = 'gray', size = 0.2, linetype = 5),
          axis.title = element_text(size = 7, color = 'black'), 
          axis.text.y = element_text(size = 6, color = 'black'), 
          axis.text.x = element_text(size = 6, color = 'black'),
          legend.text = element_text(size = 6, color = 'black'),
          legend.title = element_text(size = 7, color = 'black'),
          legend.key = element_blank()) + 
    guides(colour = guide_legend(override.aes = list(size=2), ncol = 2),
           size = guide_legend(ncol = 2))
ggsave(plot = p, path = path.up, 
       filename = paste0('GO_up.pdf'),
       height = 7.5, width = 15.5, units = 'cm')
