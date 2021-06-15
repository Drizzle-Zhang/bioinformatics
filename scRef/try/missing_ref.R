library(Seurat)
#source('scRef.R')
source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')

# To get Seurat object "pbmc", please follow the instruction of Seurat:
# https://satijalab.org/seurat/pbmc3k_tutorial.html 
setwd('/home/drizzle_zhang/my_git/scRef/zy_scripts')

pbmc.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")
# data preparing
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = F)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10, n.neighbors = 25)
# DimPlot(pbmc, reduction = "umap", label = T)

# assign cell type
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)

# find all markers of cluster 1
# cluster1.markers <- FindMarkers(pbmc, ident.1 = "Memory CD4 T", min.pct = 0.25)
# head(cluster1.markers, n = 5)

# saving
# saveRDS(pbmc, file = "./pbmc.RData", version = '3.6.0')

# loading
# load("./pbmc.RData")

COL = c()
i = 1
while(i <= length(pbmc@active.ident)) {
    this_col=which(colnames(pbmc@assays$RNA@counts)==names(pbmc@active.ident)[i])
    COL=c(COL,this_col)
    i=i+1
}      
exp_sc_mat=as.matrix(pbmc@assays$RNA@counts)[,COL]

############# MCA ref ## pbmc
exp_ref_mat=read.table(
    './PeripheralBlood_ref_human.txt',
    header=T,row.names=1,sep='\t',check.name=F)    
tag=SCREF(exp_sc_mat, exp_ref_mat)$tag2
type.tags <- as.character(tag[,2])
pbmc@meta.data$scref=tag[,2]
UMAPPlot(object = pbmc, label=T, label.size=3, group.by ='scref')

##############DEV#################
exp_ref_mat=read.table('exp_ref_mat_human_brain_dev',header=T,row.names=1,sep='\t',check.name=F)
tag=SCREF(exp_sc_mat, exp_ref_mat)$tag2
pbmc@meta.data$scref=tag[,2]
COLOR=heat.colors(n=length(table(pbmc@meta.data$scref))+2)
TSNEPlot(object = pbmc, colors.use=COLOR, group.by ='scref')

##################################      