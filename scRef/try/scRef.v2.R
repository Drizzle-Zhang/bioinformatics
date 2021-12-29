library(Seurat)
library(DESeq2)
library(foreach)
library(doParallel)
library(stringr)
library(mclust)
library(coin)
# library(pcaPP)

# path of python
Sys.setenv(RETICULATE_PYTHON = "/home/zy/tools/anaconda3/bin/python3")
# source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
source('/home/zy/my_git/bioinformatics/scRef/try/scRef.R')

# setwd('/home/drizzle_zhang/my_git/bioinformatics/scRef/zy_scripts')
setwd('/home/zy/scRef/try_data')

# input file
file.ref <- './scRef/Reference/MouseBrain_Bulk_Zhang2014/Reference_expression.txt'
file.data <- "./summary/Zeisel_exp_sc_mat.txt"
file.label.unlabeled <- './summary/Zeisel_exp_sc_mat_cluster_merged.txt'
# file.data <- "./summary/Tasic_exp_sc_mat.txt"
# file.label.unlabeled <- './summary/Tasic_exp_sc_mat_cluster_merged.txt'
# parameters
num.cpu <- 10

######################### unlabeled data
file.data <- file.data
data.unlabeled <- read.delim(file.data, row.names=1)
names(data.unlabeled) <- str_replace_all(names(data.unlabeled), '_', '.')
names(data.unlabeled) <- str_replace_all(names(data.unlabeled), '-', '.')
# read label file
file.label.unlabeled <- file.label.unlabeled
label.unlabeled <- read.delim(file.label.unlabeled, row.names=1)
row.names(label.unlabeled) <- str_replace_all(row.names(label.unlabeled), '_', '.')
row.names(label.unlabeled) <- str_replace_all(row.names(label.unlabeled), '-', '.')
col.name1 <- names(data.unlabeled)[1]
if (substring(col.name1, 1, 1) == 'X') {
    row.names(label.unlabeled) <- paste0('X', row.names(label.unlabeled))
}
# filter data
use.cols <- row.names(label.unlabeled)[label.unlabeled[,1] != 'Unclassified']
data.filter <- data.unlabeled[,use.cols]
label.filter <- data.frame(label.unlabeled[use.cols,], row.names = use.cols)
# data preparing
seurat.unlabeled <- CreateSeuratObject(counts = data.filter, project = "pbmc3k", min.cells = 3, min.features = 200)
seurat.unlabeled <- NormalizeData(seurat.unlabeled, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.unlabeled <- FindVariableFeatures(seurat.unlabeled, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat.unlabeled)
seurat.unlabeled <- ScaleData(seurat.unlabeled, features = all.genes)

# add label
seurat.unlabeled@meta.data$original.label <- 
    label.filter[dimnames(seurat.unlabeled@assays$RNA@scale.data)[[2]], 1]

# PCA
seurat.unlabeled <- RunPCA(seurat.unlabeled, features = VariableFeatures(object = seurat.unlabeled), verbose = F)

# UMAP
seurat.unlabeled <- RunUMAP(seurat.unlabeled, dims = 1:20, n.neighbors = 30)
# figure1: ture label
DimPlot(seurat.unlabeled, reduction = "umap", label = T, group.by = 'original.label')

######################### reference data
file.ref <- file.ref
exp_ref_mat <- read.table(file.ref, header=T, row.names=1, sep='\t', check.name=F)
# MCA cell name
names(exp_ref_mat) <- c("Astrocyte", "Neuron", "Oligodendrocyte precursor cell",
                        "Newly Formed Oligodendrocyte", "Myelinating oligodendrocyte",
                        "Microglia", "Endothelial cell")

######################### delete specific cell
cell.delete <- "Astrocyte"
exp_ref_mat <- exp_ref_mat[, setdiff(names(exp_ref_mat), c(cell.delete))]

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

result.scref <- SCREF(exp_sc_mat, exp_ref_mat, CPU = num.cpu)
tag <- result.scref$tag2

## try
# cell.precurcor <- tag[tag[,2] == 'Oligodendrocyte Precursor Cell','cell_id']
# corr.precurcor <- result.scref$out1[,cell.precurcor]
# out2.precurcor <- result.scref$out2[,cell.precurcor]

type.tags <- as.character(tag[,2])
seurat.unlabeled@meta.data$scref <- tag[,2]
# figure2: scRef.v1
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
cov.cell <- c(cell.MCA, names(exp_ref_mat))
mod <- model.matrix(~ as.factor(cov.cell))
mtx.combat <- ComBat(mtx.in, batch, mod, par.prior = T)

mtx.MCA <- mtx.combat[,paste0('MCA.', cell.MCA)]

cells <- names(exp_ref_mat)
cutoff.fc <- 1.5
cutoff.pval <- 0.05
list.cell.genes <- list()
for (cell in cells) {
    vec.cell <- exp_ref_mat[, cell]
    exp.top10 <- quantile(vec.cell, 0.9)
    genes.high <- gene_over[vec.cell > exp.top10]
    mtx.in <- cbind(mtx.MCA, vec.cell)
    bool.cell <- as.factor(c(rep('1', dim(mtx.MCA)[2]), '2'))
    res.limma <- getDEgeneF(mtx.in, bool.cell)
    df.diff <- res.limma[
        ((res.limma$logFC > cutoff.fc) & (res.limma$adj.P.Val < cutoff.pval)),]
    genes.diff <- row.names(df.diff)
    list.cell.genes[[cell]] <- genes.diff
}

# confirm label
exp_sc_mat <- exp_sc_mat[gene_over,]
ori.tag = seurat.unlabeled@meta.data$original.label
scRef.tag = seurat.unlabeled@meta.data$scref
meta.tag <- data.frame(ori.tag, scRef.tag, row.names = dimnames(exp_sc_mat)[[2]])
registerDoParallel(10)
confirm.classify <- function(exp_sc_mat, list.cell.genes, meta.tag, method.test, barcode) {
    expression.barcode <- exp_sc_mat[, barcode]
    bool.mark.gene <- rep(1, dim(exp_sc_mat)[1])
    cell <- meta.tag[barcode, 'scRef.tag']
    genes.marker <- list.cell.genes[[cell]]
    bool.mark.gene[names(expression.barcode) %in% genes.marker] <- 2
    test.in <- cbind(expression.barcode, bool.mark.gene)
    test.in <- as.data.frame(test.in)
    names(test.in) <- c('expression.level', 'factor.mark.gene')
    test.in$factor.mark.gene <- as.factor(test.in$factor.mark.gene)
    if (method.test == 'wilcox') {
        out.test <- wilcox.test(expression.level ~ factor.mark.gene, data = test.in, 
                                alternative = 'less')
        pvalue <- out.test$p.value
    }
    if (method.test == 't.test') {
        out.test <- t.test(expression.level ~ factor.mark.gene, data = test.in, 
                           alternative = 'less')
        pvalue <- out.test$p.value
    }
    if (method.test == 'oneway_test') {
        out.test <- oneway_test(formula = expression.level ~ factor.mark.gene, data = test.in, 
                                alternative = 'less')
        pvalue <- pvalue(out.test)
    }
    # if (method.test == 'kendall_permutation') {
    #     
    # }
    # expression.not0 <- expression.barcode[expression.barcode != 0]
    # percent.marker <- 
    #     length(intersect(names(expression.not0), genes.marker)) / length(genes.marker)
    return(data.frame(pvalue = pvalue))
    
}

method.test <- 'oneway_test'
out.par <- foreach(barcode = dimnames(exp_sc_mat)[[2]], .combine = rbind) %dopar% confirm.classify(exp_sc_mat, list.cell.genes, meta.tag, method.test, barcode)
meta.tag <- cbind(meta.tag, out.par)
meta.tag$qvalue <- p.adjust(meta.tag$pvalue, method = 'BH')

# evaluation
# uniform tags
scRef.tag[scRef.tag == "Oligodendrocyte precursor cell"] <- "Oligodendrocyte Precursor Cell"
scRef.tag[scRef.tag == "Newly Formed Oligodendrocyte"] <- "Oligodendrocyte"
scRef.tag[scRef.tag == "Myelinating oligodendrocyte"] <- "Oligodendrocyte"
scRef.tag[scRef.tag == "Endothelial cell"] <- "Endothelial Cell" 
meta.tag$scRef.tag <- scRef.tag
vec.cutoff <- c(seq(0.005, 0.1, 0.005), seq(0.15, 0.95, 0.05))
df.ARI <- data.frame(cutoff = vec.cutoff, ARI = rep(0, length(vec.cutoff)))
for (i in 1:length(vec.cutoff)) {
    cutoff <- vec.cutoff[i]
    true.tag <- meta.tag$ori.tag
    true.tag[true.tag == cell.delete] <- 'unknown'
    our.tag <- meta.tag$scRef.tag
    our.tag[meta.tag$qvalue > cutoff] <- 'unknown'
    sub.ARI <- adjustedRandIndex(true.tag, our.tag)
    df.ARI[i, 'ARI'] <- sub.ARI
}

# best cutoff
best.cutoff <- df.ARI$cutoff[df.ARI$ARI == max(df.ARI$ARI)][1]
new.tag <- meta.tag$scRef.tag
new.tag[meta.tag$qvalue > best.cutoff] <- 'unknown'
meta.tag$new.tag <- new.tag
print(best.cutoff)

# no scRef plus
true.tag <- meta.tag$ori.tag
true.tag[true.tag == cell.delete] <- 'unknown'
adjustedRandIndex(true.tag, scRef.tag)
adjustedRandIndex(true.tag, new.tag)

# new UMAP
seurat.unlabeled@meta.data$scref.v2 <- new.tag
UMAPPlot(object = seurat.unlabeled, label = T, label.size = 3, group.by = 'scref.v2')

# remove neuron
meta.tag.remove <- meta.tag[meta.tag$ori.tag != "Neuron", c("ori.tag", "scRef.tag", "qvalue")]
vec.cutoff <- c(seq(0.005, 0.1, 0.005), seq(0.15, 0.95, 0.05))
df.ARI.remove <- data.frame(cutoff = vec.cutoff, ARI = rep(0, length(vec.cutoff)))
for (i in 1:length(vec.cutoff)) {
    cutoff <- vec.cutoff[i]
    true.tag <- meta.tag.remove$ori.tag
    true.tag[true.tag == cell.delete] <- 'unknown'
    our.tag <- meta.tag.remove$scRef.tag
    our.tag[meta.tag.remove$qvalue > cutoff] <- 'unknown'
    sub.ARI <- adjustedRandIndex(true.tag, our.tag)
    df.ARI.remove[i, 'ARI'] <- sub.ARI
}
# best cutoff
best.cutoff.remove <- df.ARI.remove$cutoff[df.ARI.remove$ARI == max(df.ARI.remove$ARI)][1]
new.tag.remove <- meta.tag.remove$scRef.tag
new.tag.remove[meta.tag.remove$qvalue > best.cutoff] <- 'unknown'
meta.tag.remove$new.tag <- new.tag.remove
print(best.cutoff.remove)

# no scRef plus
true.tag.remove <- meta.tag.remove$ori.tag
true.tag.remove[true.tag.remove == cell.delete] <- 'unknown'
adjustedRandIndex(true.tag.remove, meta.tag.remove$scRef.tag)
adjustedRandIndex(true.tag.remove, new.tag.remove)



