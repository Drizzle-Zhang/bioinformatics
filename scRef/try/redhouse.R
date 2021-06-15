library(stringr)
library(reticulate)

Sys.setenv(RETICULATE_PYTHON = "/home/zy/tools/anaconda3/bin/python3")
source('/home/zy/my_git/bioinformatics/scRef/try/scRef.R')

num.cpu <- 10

# reference (science advance)
exp_ref_mat=read.table('/home/yjingjing/project/hfz2020adjust/pla_adjust/4th_combNA_SA_scRef/SA.villi.csv',
                       header=T,row.names=1,sep='\t')
exp_ref_mat.origin <- exp_ref_mat
# HCA cell name
names(exp_ref_mat.origin) <- c("Villous trophoblast cell", "Syncytiotrophoblast cell", 
                               "Extravillous trophoblast", "Hofbauer cell", 
                               "Erythroid cell", "Fibroblast cell", "Fibroblast cell", 
                               "Fibroblast cell", "Vascular endothelial cell")

find.markers <- function(exp_ref_mat) {
    ###### regard MCA as reference of DEG
    file.MCA <- '/home/disk/scRef/HumanAtlas_SingleCell_Han2020/combinedHCA/HCA_combined_by_cell.txt'
    df.MCA <- read.table(file.MCA, header=T, row.names=1, sep='\t', check.name=F)
    df.MCA[is.na(df.MCA)] <- 0
    df.MCA <- df.MCA[rowSums(df.MCA) != 0,]
    library(DESeq2)
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
    # mod <- model.matrix(~ 1)
    mtx.combat <- ComBat(mtx.in, batch, mod, par.prior = T)
    
    mtx.MCA <- mtx.combat[,paste0('MCA.', cell.MCA)]
    mtx.ref <- mtx.combat[,paste0('Ref.', cell.ref)]
    dimnames(mtx.ref)[[2]] <- cell.ref
    
    cells <- cell.ref
    cutoff.fc <- 1
    cutoff.pval <- 0.05
    list.cell.genes <- list()
    for (cell in cells) {
        vec.cell <- mtx.ref[, cell]
        exp.top10 <- quantile(vec.cell, 0.9)
        genes.high <- gene_over[vec.cell > exp.top10]
        mtx.in <- cbind(mtx.MCA, vec.cell)
        bool.cell <- as.factor(c(rep('1', dim(mtx.MCA)[2]), '2'))
        res.limma <- .getDEgeneF(mtx.in, bool.cell)
        # df.diff <- res.limma[
        #     ((res.limma$logFC > cutoff.fc) & (res.limma$adj.P.Val < cutoff.pval)),]
        genes.diff <- row.names(res.limma)[1:100]
        list.cell.genes[[cell]] <- genes.diff
    }
    
    out <- list()
    out[['list.cell.genes']] <- list.cell.genes
    out[['exp_ref_mat']] <- mtx.ref
    return(out)
    
}

out.markers <- find.markers(exp_ref_mat.origin)
list.cell.genes <- out.markers[['list.cell.genes']]
genes.ref <- dimnames(out.markers[['exp_ref_mat']])[[1]]

######################### unlabeled data
file.data <- '/home/yjingjing/project/hfz2020adjust/pla_adjust/4th_combNA_SA_scRef/6701.deci_pla'
file.label.unlabeled <- '/home/yjingjing/project/hfz2020adjust/pla_adjust/4th_combNA_SA_scRef/E-MTAB-6701.processed.2'
file.data <- file.data
data.unlabeled <- read.delim(file.data, row.names=1)
colnames <- names(data.unlabeled)
list.colnames <- strsplit(colnames, '.', fixed = T)
colnames <- c()
for (i in 1:length(list.colnames)) {
    sub_colname <- list.colnames[[i]]
    colnames <- c(colnames, paste0(strsplit(sub_colname[1], '_')[[1]][1:2], collapse = '.'))
}
names(data.unlabeled) <- colnames
# read label file
file.label.unlabeled <- file.label.unlabeled
label.unlabeled <- read.delim(file.label.unlabeled, row.names=1)
row.names(label.unlabeled) <- str_replace_all(row.names(label.unlabeled), '_', '.')
# filter data
use.cols <- row.names(label.unlabeled)[label.unlabeled$location == 'Placenta']
data.filter <- data.unlabeled[,use.cols]
label.filter <- data.frame(label.unlabeled[use.cols, 'annotation'], row.names = use.cols)
# label.filter <- data.frame(label.unlabeled[use.cols,], row.names = use.cols)

# get overlap genes
exp_sc_mat = data.filter
exp_ref_mat <- exp_ref_mat.origin[genes.ref,]
exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
gene_sc=rownames(exp_sc_mat)
gene_ref=rownames(exp_ref_mat)
gene_over= gene_sc[which(gene_sc %in% gene_ref)]
exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
print('Number of overlapped genes:')
print(nrow(exp_sc_mat))

######################### run scRef
result.scref <- SCREF(exp_sc_mat, exp_ref_mat, CPU = num.cpu)
tag <- result.scref$tag2

# confirm label
exp_sc_mat <- exp_sc_mat[gene_over,]
ori.tag = label.filter[names(exp_sc_mat), 1]
scRef.tag = tag[,2]
method.test <- 'wilcox'
# method.test <- 't.test'
# method.test <- 'oneway_test'
meta.tag <- comfirm.label(exp_sc_mat, ori.tag, scRef.tag, method.test)

# evaluation
# import python package: sklearn.metrics
use_condaenv("/home/zy/tools/anaconda3")
# py_config()
metrics <- import('sklearn.metrics')
# uniform tags
scRef.tag[scRef.tag == "Astrocyte"] <- "astrocytes_ependymal"
scRef.tag[scRef.tag == "Newly Formed Oligodendrocyte"] <- "oligodendrocytes"
scRef.tag[scRef.tag == "Myelinating oligodendrocyte"] <- "oligodendrocytes"
scRef.tag[scRef.tag == "Endothelial cell"] <- "endothelial-mural" 
scRef.tag[scRef.tag == "Neuron"] <- "neurons" 
scRef.tag[scRef.tag == "Microglia"] <- "microglia" 
cell.delete <- "astrocytes_ependymal"
meta.tag$scRef.tag <- scRef.tag
vec.cutoff <- c(seq(0.005, 0.1, 0.005), seq(0.15, 0.95, 0.05))
df.metrics <- data.frame(cutoff = vec.cutoff, weighted.f1 = rep(0, length(vec.cutoff)))
for (i in 1:length(vec.cutoff)) {
    cutoff <- vec.cutoff[i]
    true.tag <- meta.tag$ori.tag
    true.tag[true.tag == cell.delete] <- 'unknown'
    our.tag <- meta.tag$scRef.tag
    our.tag[meta.tag$qvalue > cutoff] <- 'unknown'
    sub.weighted.f1 <- metrics$f1_score(true.tag, our.tag, average = 'weighted')
    df.metrics[i, 'weighted.f1'] <- sub.weighted.f1
}

# best cutoff
best.cutoff <- df.metrics$cutoff[df.metrics$weighted.f1 == max(df.metrics$weighted.f1)][1]
best.cutoff <- 0.01
new.tag <- meta.tag$scRef.tag
new.tag[meta.tag$qvalue > best.cutoff] <- 'unknown'
meta.tag$new.tag <- new.tag
print(best.cutoff)

# no scRef plus
true.tag <- meta.tag$ori.tag
true.tag[true.tag == cell.delete] <- 'unknown'
meta.tag$ori.tag <- true.tag
metrics$f1_score(true.tag, scRef.tag, average = 'weighted')
metrics$f1_score(true.tag, new.tag, average = 'weighted')

###
fb=meta.tag[meta.tag$ori.tag %in% c('fFB1', 'fFB2'),]

### plot
library(Seurat)
# data preparing
seurat.unlabeled <- CreateSeuratObject(counts = data.filter, project = "pbmc3k", min.cells = 3, min.features = 200)
seurat.unlabeled <- NormalizeData(seurat.unlabeled, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.unlabeled <- FindVariableFeatures(seurat.unlabeled, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat.unlabeled)
seurat.unlabeled <- ScaleData(seurat.unlabeled, features = all.genes)

# add label
seurat.unlabeled@meta.data$original.label <- ori.tag
seurat.unlabeled@meta.data$scRef.tag <- scRef.tag
seurat.unlabeled@meta.data$new.tag <- new.tag

# PCA
seurat.unlabeled <- RunPCA(seurat.unlabeled, features = VariableFeatures(object = seurat.unlabeled), verbose = F)

# UMAP
seurat.unlabeled <- RunUMAP(seurat.unlabeled, dims = 1:15, n.neighbors = 30)
# figure1: ture label
DimPlot(seurat.unlabeled, reduction = "umap", label = T, group.by = 'original.label')
# figure2: scRef label
DimPlot(seurat.unlabeled, reduction = "umap", label = T, group.by = 'scRef.tag')
# figure3: scRef plus label
DimPlot(seurat.unlabeled, reduction = "umap", label = T, group.by = 'new.tag')

