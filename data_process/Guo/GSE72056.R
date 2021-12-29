library(dplyr)
library(Seurat)
library(ggplot2)

# setwd('/home/drizzle_zhang/data_process/GSE72056')
setwd('/home/zy/data_processing/guo/GSE72056')

df.input <- read.table('GSE72056_melanoma_single_cell_revised_v2.txt', sep = '\t',
                       stringsAsFactors = F)
malignant.pre <- as.numeric(df.input[3,2:dim(df.input)[2]])
malignant <- malignant.pre[malignant.pre != 0]
malignant <- factor(malignant, levels = c(1, 2, 0), 
                    labels = c("N", "T", "unresolved"))
cell.type <- as.numeric(df.input[4,2:dim(df.input)[2]])
cell.type <- cell.type[malignant.pre != 0]
cell.type <- factor(cell.type, levels = c(1, 2, 3, 4, 5, 6), 
                    labels = c("T", "B", "Macro", "Endo.", "CAF", "NK"))
gene <- df.input[5:dim(df.input)[1], 1]
gene <- sub("_", "-", gene, fixed = T)
cell <- as.character(as.factor(df.input[1, 2:dim(df.input)[2]]))
df.mtx <- df.mtx[5:dim(df.input)[1], 2:dim(df.input)[2]]
df.mtx <- df.mtx[!duplicated(gene),]
df.mtx <- df.mtx[,malignant.pre != 0]
rownames(df.mtx) <- gene[!duplicated(gene)]
colnames(df.mtx) <- cell[malignant.pre != 0]


pbmc <- CreateSeuratObject(counts = df.mtx, project = "pbmc3k", 
                           min.cells = 0, min.features = 0)
pbmc@meta.data$batch <- malignant
pbmc@meta.data$group <- cell.type

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes, vars.to.regress = c("batch",'nCount_RNA'))


pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)


saveRDS(pbmc,file = 'pbmc.RDS')
# pbmc <- readRDS('pbmc.RDS')

umap.plot <- DimPlot(pbmc, reduction = "umap", group.by = 'batch', label = T)
ggsave(
    plot = umap.plot, path = './', filename = "UMAP.png",
    units = 'cm', width = 30, height = 30)


pbmc@meta.data$type = 
    paste(as.character(pbmc@meta.data$batch),
          as.character(pbmc@meta.data$group), sep = '_')
feature.plot <- 
    FeaturePlot(pbmc, features = c("CD3G",'CD3D','CD3E','CD4','FOXP3'))
ggsave(
    plot = feature.plot, path = './', filename = "feature.png",
    units = 'cm', width = 20, height = 30)


#show_gene=c('CD3g', 'CD4','FOXP3','CDC42','WAS','GATA3',"CA1","IL4","IL17")

#VlnPlot(pbmc, features = show_gene,slot='data',group.by='type',idents=c('T_Tcell','N_Tcell'))



used <- which(pbmc@meta.data$type %in% c('T_T','N_T'))
EXP <- as.matrix(pbmc@assays$RNA@data[,used])
colnames(EXP) <- paste0(pbmc@meta.data$batch[used],'_',colnames(EXP))


CD3G <- which(rownames(EXP) %in% c('CD3G'))
CD3D <- which(rownames(EXP) %in% c('CD3D'))
CD3E <- which(rownames(EXP) %in% c('CD3E'))


CD4 <- which(rownames(EXP) %in% c('CD4'))
FOXP3 <- which(rownames(EXP) %in% c('FOXP3'))

TMP <- EXP[cbind(CD3G,CD3D,CD3E,CD4,FOXP3),]
TMP[which(TMP > 0)] <- 1
#which()

#PAT <- EXP[which(rownames(EXP) %in% show_gene),]
pdf('heatmap.pdf', width = 10, height = 7)
heatmap(TMP,scale = 'none',margins = c(15,10),Colv = F)
dev.off()


table(pbmc@meta.data$type)
#T_Tcell <- 4
#N_Tcell <- 2040
MIN <- apply(TMP[c(4,5),],2,min)

POS <- which(MIN > 0)
NEG <- which(MIN == 0)
EXP <- cbind(EXP[,NEG], EXP[,POS])
colnames(EXP)[ncol(EXP)] <- names(POS)


#EXP[,which(colnames(EXP)==POS)]
PT <- rep('NEG',ncol(EXP))
PT[which(colnames(EXP) == names(POS))] <- 'POS'
PT <- t(as.matrix(PT))


OUT <- cbind(rownames(EXP),rep('NO',nrow(EXP)),EXP)
colnames(OUT)[c(1,2)] <- c('GENE','DESCRIPTION')
write.table(OUT,'EXP.txt',sep = '\t',quote = F,row.names = F,col.names = T)
write.table(PT,'PT.cls',sep = ' ',quote = F,row.names = F,col.names = F)
