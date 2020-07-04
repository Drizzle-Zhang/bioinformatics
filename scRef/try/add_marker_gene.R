library(Seurat)
library("DESeq2")
library(foreach)
library(doParallel)
#source('scRef.R')

# path of python
Sys.setenv(RETICULATE_PYTHON = "/home/zy/tools/anaconda3/bin/python3")
source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')

# setwd('/home/drizzle_zhang/my_git/bioinformatics/scRef/zy_scripts')
setwd('/home/zy/scRef/try_data')

######################### unlabeled data
file.data <- "./summary/Zeisel_exp_sc_mat.txt"
data.unlabeled <- read.delim(file.data, row.names=1)
# data preparing
seurat.unlabeled <- CreateSeuratObject(counts = data.unlabeled, project = "pbmc3k", min.cells = 3, min.features = 200)
seurat.unlabeled <- NormalizeData(seurat.unlabeled, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.unlabeled <- FindVariableFeatures(seurat.unlabeled, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat.unlabeled)
seurat.unlabeled <- ScaleData(seurat.unlabeled, features = all.genes)

# add label
file.label.unlabeled <- './summary/Zeisel_exp_sc_mat_cluster_merged.txt'
label.unlabeled <- read.delim(file.label.unlabeled, row.names=1)
row.names(label.unlabeled) <- paste0('X', row.names(label.unlabeled))
seurat.unlabeled@meta.data$original.label <- 
    label.unlabeled[dimnames(seurat.unlabeled@assays$RNA@scale.data)[[2]], 1]

# PCA
seurat.unlabeled <- RunPCA(seurat.unlabeled, features = VariableFeatures(object = seurat.unlabeled), verbose = F)
seurat.unlabeled <- FindNeighbors(seurat.unlabeled, dims = 1:10)
seurat.unlabeled <- FindClusters(seurat.unlabeled, resolution = 0.5)

# UMAP
seurat.unlabeled <- RunUMAP(seurat.unlabeled, dims = 1:10, n.neighbors = 30)
DimPlot(seurat.unlabeled, reduction = "umap", label = T, group.by = 'original.label')

######################### reference data
file.ref <- './scRef/Reference/MouseBrain_Bulk_Zhang2014/Reference_expression.txt'
exp_ref_mat <- read.table(file.ref, header=T, row.names=1, sep='\t', check.name=F)
exp_ref_mat <- exp_ref_mat[, setdiff(names(exp_ref_mat), c("Neuron"))]

######################### run scRef
# select non-filtered cells
COL = c()
i = 1
while (i <= length(seurat.unlabeled@active.ident)) {
    this_col <- 
        which(colnames(seurat.unlabeled@assays$RNA@counts) == 
                  names(seurat.unlabeled@active.ident)[i])
    COL <- c(COL,this_col)
    i <- i + 1
}      
exp_sc_mat <- as.matrix(seurat.unlabeled@assays$RNA@counts)[,COL]

result.scref <- SCREF(exp_sc_mat, exp_ref_mat, CPU = 8)
tag <- result.scref$tag2

## try
# cell.precurcor <- tag[tag[,2] == 'Oligodendrocyte Precursor Cell','cell_id']
# corr.precurcor <- result.scref$out1[,cell.precurcor]
# out2.precurcor <- result.scref$out2[,cell.precurcor]

type.tags <- as.character(tag[,2])
seurat.unlabeled@meta.data$scref <- tag[,2]
UMAPPlot(object = seurat.unlabeled, label = T, label.size = 3, group.by = 'scref')

# different expression analysis 
# get overlap genes
exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
gene_sc=rownames(exp_sc_mat)
gene_ref=rownames(exp_ref_mat)
gene_over= gene_sc[which(gene_sc %in% gene_ref)]
exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
print('Number of overlapped genes:')
print(nrow(exp_sc_mat))

# limma function
getDEgeneF <- function(esetm=NULL,group=NULL,pair=FALSE,block=NULL,p_adj="fdr",fpkm=T){
    if (is.null(esetm)) {
        cat("esetm: gene expression matrix",
            "group: factor: \"c\"/\"d\"",
            "pair: TRUE/FALSE*",
            "block: e.g.1 2 2 1 if paired; blank if not",
            "p_adj: p.adjust, fdr* ",
            "fpkm: TRUE/FALSE*",        
            sep="\n")
    }else{
        library(limma)
        if(pair){
            design<-model.matrix(~block+group)
        }else{
            design<-model.matrix(~group)
        }
        fit<-lmFit(esetm,design)
        if(fpkm){
            fit<-eBayes(fit,trend=T,robust=T)
        }else{
            fit<-eBayes(fit)
        }
        x<-topTable(fit,number=nrow(esetm),adjust.method=p_adj,coef="group2")
        x<-x[!is.na(row.names(x)),]
        x<-x[!duplicated(row.names(x)),]
        return(x)    
    }
}

######################## DEseq
# countdata <- as.matrix(exp_ref_mat)
# cells <- names(exp_ref_mat)
# # class(countdata)
# # head(countdata)
# coldata <- DataFrame(row.names = cells)
# for (cell in cells) {
#     bool.cell <- cells
#     bool.cell[bool.cell != cell] <- 'other'
#     bool.cell <- as.factor(bool.cell)
#     coldata$bool.cell <- bool.cell
#     ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
#                                      colData = coldata,
#                                      design = ~ 1 + bool.cell)
#     
# }

######################### limma
# cells <- names(exp_ref_mat)
# cutoff.fc <- 2
# cutoff.pval <- 0.01
# list.cell.genes <- list()
# for (cell in cells) {
#     bool.cell <- cells
#     bool.cell[bool.cell == cell] <- '2'
#     bool.cell[bool.cell != '2'] <- '1'
#     bool.cell <- as.factor(bool.cell)
#     res.limma <- getDEgeneF(exp_ref_mat, bool.cell)
#     df.diff <- res.limma[
#         ((res.limma$logFC > cutoff.fc) & (res.limma$adj.P.Val < cutoff.pval)),]
#     genes.diff <- row.names(df.diff)
#     list.cell.genes[[cell]] <- genes.diff
# }

###### regard MCA as reference of DEG
file.MCA <- '/home/disk/scRef/MouseAtlas_SingleCell_Han2018/combinedMCA/MCA_combined_mouse.txt'
df.MCA <- read.table(file.MCA, header=T, row.names=1, sep='\t', check.name=F)
coldata.MCA <- DataFrame(row.names = names(df.MCA))
obj.DESeq.MCA <- DESeqDataSetFromMatrix(countData = df.MCA, colData = coldata.MCA, 
                                        design = ~ 1)
fpm.MCA <- fpm(obj.DESeq.MCA, robust = T)
# overlap genes
fpm.MCA=fpm.MCA[order(rownames(fpm.MCA)),]
exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
gene_MCA=rownames(fpm.MCA)
gene_ref=rownames(exp_ref_mat)
gene_over= gene_MCA[which(gene_MCA %in% gene_ref)]
fpm.MCA=fpm.MCA[which(gene_MCA %in% gene_over),]
exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
print('Number of overlapped genes:')
print(nrow(exp_ref_mat))

cell.MCA <- dimnames(fpm.MCA)[[2]]
cell.ref <- names(exp_ref_mat)
cell.overlap <- intersect(cell.MCA, cell.ref)
# combat
library(sva)
mtx.in <- cbind(fpm.MCA, exp_ref_mat)
names(mtx.in) <- c(paste0('MCA.', cell.MCA), paste0('Ref.', cell.ref))
batch <- c(rep(1, dim(fpm.MCA)[2]), rep(2, dim(exp_ref_mat)[2]))
cov.cell <- c(cell.MCA, "Astrocyte", "Neuron", "Oligodendrocyte precursor cell",
              "Newly Formed Oligodendrocyte", "Myelinating oligodendrocyte", 
              "Microglia")
mod <- model.matrix(~ as.factor(cov.cell))
# bpparam <- MulticoreParam(6, RNGseed = 7739465)
mtx.combat <- ComBat(mtx.in, batch, mod, par.prior = T)

# b <- mtx.combat[,c('MCA.Microglia', 'Ref.Microglia')]
# cor.fk(mtx.combat[,'MCA.Microglia'], mtx.combat[,'Ref.Microglia'])

mtx.MCA <- mtx.combat[,paste0('MCA.', cell.MCA)]

cells <- names(exp_ref_mat)
cutoff.fc <- 2
cutoff.pval <- 0.01
list.cell.genes <- list()
for (cell in cells) {
    mtx.in <- cbind(mtx.MCA, exp_ref_mat[,cell])
    bool.cell <- as.factor(c(rep('1', dim(mtx.MCA)[2]), '2'))
    res.limma <- getDEgeneF(mtx.in, bool.cell)
    df.diff <- res.limma[
        ((res.limma$logFC > cutoff.fc) & (res.limma$adj.P.Val < cutoff.pval)),]
    genes.diff <- row.names(df.diff)
    list.cell.genes[[cell]] <- genes.diff
}

# similarity of DEGs
# mtx.similar <- matrix(1:36, nrow = 6)
# dimnames(mtx.similar)[[1]] <- cells
# dimnames(mtx.similar)[[2]] <- cells
# for (cell1 in cells) {
#     for (cell2 in cells) {
#         mtx.similar[cell1, cell2] <- 
#             length(intersect(list.cell.genes[[cell1]], list.cell.genes[[cell2]])) / 
#             min(length(list.cell.genes[[cell1]]), length(list.cell.genes[[cell2]]))
#     }
# }

# confirm label
exp_sc_mat <- exp_sc_mat[gene_over,]
ori.tag = seurat.unlabeled@meta.data$original.label
scRef.tag = seurat.unlabeled@meta.data$scref
meta.tag <- data.frame(ori.tag, scRef.tag, row.names = dimnames(exp_sc_mat)[[2]])
registerDoParallel(10)
confirm.classify <- function(exp_sc_mat, list.cell.genes, meta.tag, barcode) {
    expression.barcode <- exp_sc_mat[, barcode]
    bool.mark.gene <- rep(1, dim(exp_sc_mat)[1])
    cell <- meta.tag[barcode, 'scRef.tag']
    genes.marker <- list.cell.genes[[cell]]
    bool.mark.gene[names(expression.barcode) %in% genes.marker] <- 2
    wlicox.in <- cbind(expression.barcode, as.factor(bool.mark.gene))
    wlicox.in <- as.data.frame(wlicox.in)
    names(wlicox.in) <- c('expression.level', 'factor.mark.gene')
    out.wilcox <- wilcox.test(expression.level ~ factor.mark.gene, data = wlicox.in, 
                              alternative = 'less')
    pvalue <- out.wilcox$p.value
    expression.not0 <- expression.barcode[expression.barcode != 0]
    percent.marker <- 
        length(intersect(names(expression.not0), genes.marker)) / length(genes.marker)
    return(data.frame(pvalue = pvalue, percent.marker = percent.marker))
}

out.par <- foreach(barcode = dimnames(exp_sc_mat)[[2]], .combine = rbind) %dopar% confirm.classify(exp_sc_mat, list.cell.genes, meta.tag, barcode)
meta.tag <- cbind(meta.tag, out.par)
meta.tag$qvalue <- p.adjust(meta.tag$pvalue, method = 'BH')
new.tag <- meta.tag$scRef.tag
new.tag[meta.tag$qvalue > 0.05] <- 'unknown'
meta.tag$new.tag <- new.tag

# new UMAP
seurat.unlabeled@meta.data$scref.v2 <- new.tag
UMAPPlot(object = seurat.unlabeled, label = T, label.size = 3, group.by = 'scref.v2')
