## 在生科院服务器上产生.loom文件(JingMA/rawData/bam)

.libPaths(c("/home/yzj/R/x86_64-pc-linux-gnu-library/4.0","/home/zy/tools/R-4.0.0/library"))
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

file.seurat <- '/home/disk/drizzle/wgk/data/marker_1.5_merge/seurat_first.Rdata'
pbmc <- readRDS(file.seurat)
Idents(pbmc) <- pbmc$sample

batch='C1'
if(batch=='C'){
  pbmc_c1 <- subset(pbmc,idents='C1')
  cells_c1 <- gsub('_',':',colnames(pbmc_c1))
  cells_c1 <- gsub('-1','x',cells_c1)
  
  pbmc_c2 <- subset(pbmc,idents='C2')
  cells_c2 <- gsub('_',':',colnames(pbmc_c2))
  cells_c2 <- gsub('-1','x',cells_c2)
  
  pbmc_c3 <- subset(pbmc,idents='C3')
  cells_c3 <- gsub('_',':',colnames(pbmc_c3))
  cells_c3 <- gsub('-1','x',cells_c3)
  
  pbmc_c4 <- subset(pbmc,idents='C4')
  cells_c4 <- gsub('_',':',colnames(pbmc_c4))
  cells_c4 <- gsub('-1','x',cells_c4)
  
  pbmc <- subset(pbmc,idents=c('C1','C2','C3','C4'))
  sub_cells <- c(cells_c1,cells_c2,cells_c3,cells_c4)
  rm(pbmc_c1);rm(pbmc_c2);rm(pbmc_c3);rm(pbmc_c4)
}else{
  pbmc <- subset(pbmc,idents=batch)
  sub_cells <- gsub('_',':',colnames(pbmc))
  sub_cells <- gsub('-1','x',sub_cells)
}

pbmc <- subset(pbmc,idents=c('C1','C2','C3'))
sub_cells <- gsub('_',':',colnames(pbmc))
sub_cells <- gsub('-1','x',sub_cells)

input <- ReadVelocity(file = paste('/home/yzj/JingMA_ORI/res/RNAvelocity/',batch,'/',batch,'.loom',sep = ''))
RV <- as.Seurat(x = input)
RV <- subset(RV,cells=sub_cells)
RV <- SCTransform(object = RV, assay = "spliced")
saveRDS(RV, '/home/disk/drizzle/wgk/trajectory/RNA_volocity/RV_SCT_first_normal.RDS')


RV <- readRDS('/home/disk/drizzle/wgk/trajectory/RNA_volocity/RV_SCT_first_normal.RDS')
RV <- RunPCA(object = RV, seed.use=123, npcs=100, features = VariableFeatures(object = RV), ndims.print=1,nfeatures.print=1)
RV <- RunUMAP(RV, dims = 1:100,n.components=2,n.neighbors = 30)

# RV@reductions$pca <- pbmc@reductions$pca
# RV@reductions$umap <- pbmc@reductions$umap

RV@reductions$umap@cell.embeddings[,1] <- pbmc@reductions$umap@cell.embeddings[,1]
RV@reductions$umap@cell.embeddings[,2] <- pbmc@reductions$umap@cell.embeddings[,2]
#RV@meta.data$Cluster <- as.character(pbmc@meta.data$CellType_ID)
#RV@meta.data$Cluster <- factor(RV@meta.data$Cluster,levels=c('C0','C1','C2','C3','C4'))

DimPlot(pbmc, group.by = 'vec.cluster', label = T)

RV@meta.data$Cluster <- pbmc@meta.data$vec.cluster
Idents(RV) <- RV@meta.data$Cluster

DimPlot(RV, group.by = 'Cluster', label = T)
RV <- RunVelocity(object = RV, deltaT = 1, kCells = 25, fit.quantile = 0.1,n.cores=25)
file.volocity <- '/home/disk/drizzle/wgk/trajectory/RNA_volocity/RV_volocity2.RDS'
saveRDS(RV,file.volocity)


RV <- readRDS(file.volocity)

RV$cluster_1.5 <- pbmc@meta.data$RNA_snn_res.1.5
DimPlot(RV, group.by = 'cluster_1.5', label = T)
RV <- subset(RV, subset = cluster_1.5 != 23)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = RV)))
names(x = ident.colors) <- levels(x = RV)
cell.colors <- ident.colors[Idents(object = RV)]
names(x = cell.colors) <- colnames(x = RV)

for(n in c(100)){
  N=n
  for(s in c('sqrt')){
    scale=s
    pdf(paste('/home/disk/drizzle/wgk/trajectory/RNA_volocity/RV_',n,'_',scale,'2.pdf',sep=''))
    p <- show.velocity.on.embedding.cor(emb = Embeddings(object = RV, reduction = "umap"),
                                        vel = Tool(object = RV, slot = "RunVelocity"),
                                        n = N, scale = scale, cell.colors = ac(x = cell.colors, alpha = 0.5),
                                        cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5,
                                        grid.n = 40, arrow.lwd = 1, n.cores=30,
                                        do.par = FALSE, cell.border.alpha = 0.1)
    dev.off()
  } 
}



