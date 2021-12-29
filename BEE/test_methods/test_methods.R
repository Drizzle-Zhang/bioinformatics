setwd("~/my_git/bioinformatics/BEE/test_methods")
source('functions_evaluator.R')
library(ggplot2)

# read data
root.path <- '/home/drizzle_zhang/Desktop/single_cell/BEE/real_data'
folder.data <- '10X_balance_subset1'
pc.num <- 50
path.in <- paste(root.path, folder.data, sep = '/')
path <- paste(path.in, 'output5', sep = '/')

dir.create(path)
file.mtx <- paste(path.in, 'matrix.txt', sep = '/')
file.meta <- paste(path.in, 'meta.txt', sep = '/')
file.info <- 
    '/home/drizzle_zhang/Desktop/single_cell/BEE/real_data/lab_info.txt'
df.mtx <- read.table(file.mtx, sep = '\t', header = T, stringsAsFactors = F)
row.names(df.mtx) <- df.mtx[, 'Gene']
df.mtx['Gene'] <- NULL
df.mtx <- na.omit(df.mtx)

df.meta <- read.table(file.meta, sep = '\t', header = T, stringsAsFactors = F)
df.meta$column <- names(df.mtx)
# select data
# select_folders <- c('MouseArcME_SingleCell_JN2017',
#                     'MouseCerebralcortex_SingleCell_Loo2018')
select_folders <- c('MouseEmbryo_SingleCell_PijuanSala2019',
                    'MouseGastrulation_SingleCell_Ximena2017')
df.meta <- df.meta[df.meta$folder %in% select_folders,]
cells <- df.meta$column
df.mtx <- df.mtx[, cells]

groups <- df.meta$cell
df.info <- read.delim(file.info, stringsAsFactors = F)
batches <- c()
for (batch in df.meta$folder) {
    batches <- c(batches, df.info[df.info$folder == batch, 'author'])
}


label <- c()
timing <- c()
df.result <- data.frame()
# origin
object <- preprocess.data(df.mtx, batches, groups, num.pc = pc.num)
plot.origin <- DimPlot(
    object, reductions = 'umap', 
    group.by = "group", shape.by = 'batch', pt.size = 1.5
)
ggsave(plot = plot.origin, filename = 'UMAP_origin.png', path = path,
       units = 'cm', height = 15, width = 25)
# result.origin <- evaluate.two.dims(object)
result.origin <- evaluate.sil(object)
df.result <- rbind(df.result, result.origin)
timing <- c(timing, 0)
label <- c(label, 'Origin')

# normalized data
mtx.norm <- as.matrix(object@assays$RNA@data)[object@assays$RNA@var.features,]

# combat
time1 <- Sys.time()
library(sva)
mtx.combat <- ComBat(dat = mtx.norm, batch = object@meta.data$batch)
time2 <- Sys.time()
object.combat <- CreateSeuratObject(counts = mtx.combat)
object.combat@assays$RNA <- CreateAssayObject(data = mtx.combat)
object.combat@commands$NormalizeData.RNA <- object@commands$NormalizeData.RNA
object.combat@commands$FindVariableFeatures.RNA <- 
    object@commands$FindVariableFeatures.RNA
object.combat@assays$RNA@var.features <- object@assays$RNA@var.features
object.combat <- 
    preprocess.data(object.combat, batches, groups, num.pc = pc.num)
plot.combat <- DimPlot(
    object.combat, reductions = 'umap', 
    group.by = "group", shape.by = 'batch', pt.size = 1.5
)
ggsave(plot = plot.combat, filename = 'UMAP_combat.png', path = path,
       units = 'cm', height = 15, width = 25)
# result.combat <- evaluate.two.dims(object.combat)
result.combat <- evaluate.sil(object.combat)
df.result <- rbind(df.result, result.combat)
timing <- c(timing, difftime(time2, time1, units = 'secs'))
label <- c(label, 'ComBat')

# BEER
time3 <- Sys.time()
source('/home/drizzle_zhang/tools/BEER/BEER.R')
mybeer <- BEER(mtx.norm, object@meta.data$batch, GNUM = 30, PCNUM = pc.num, 
               ROUND = 1, GN = 2000, SEED = 1, COMBAT = TRUE)
time4 <- Sys.time()
object.beer <- CreateSeuratObject(counts = mtx.norm)
object.beer@reductions$pca <- CreateDimReducObject(
    embeddings = mybeer$seurat@reductions$pca@cell.embeddings,
    loadings = mybeer$seurat@reductions$pca@feature.loadings,
    key = 'PC_', assay = 'RNA'
)
object.beer <- RunUMAP(
    object.beer, reduction = 'pca', dims = 1:pc.num, n.components = 2, 
    verbose = F)
object.beer@meta.data <- object@meta.data
plot.beer <- DimPlot(
    object.beer, reductions = 'umap', 
    group.by = "group", shape.by = 'batch', pt.size = 1.5
)
ggsave(plot = plot.beer, filename = 'UMAP_beer.png', path = path,
       units = 'cm', height = 15, width = 25)
# result.beer <- evaluate.two.dims(object.beer)
result.beer <- evaluate.sil(object.beer)
df.result <- rbind(df.result, result.beer)
timing <- c(timing, difftime(time4, time3, units = 'secs'))
label <- c(label, 'BEER')

# MNN
time5 <- Sys.time()
library(scran)
# library(psych)
# library(batchelor)
# library(SingleCellExperiment)
# sce <- SingleCellExperiment(assays = list(logcounts = mtx.norm))
# out.mnn <- fastMNN(sce, batch = as.factor(object@meta.data$batch),
#                     d = pc.num, k = 20)
input.mnn <- list()
cell.names <- dimnames(object@assays$RNA@counts)[[2]]
input.cells <- c()
for (lab in unique(object@meta.data$batch)) {
    lab.cells <- cell.names[object@meta.data$batch == lab]
    input.mnn[[lab]] <- mtx.norm[,lab.cells]
    input.cells <- c(input.cells, lab.cells)
}
if (length(input.mnn) == 2) {
    out.mnn <- fastMNN(input.mnn[[1]], input.mnn[[2]], d = pc.num, k = 20)
}
if (length(input.mnn) == 3) {
    out.mnn <- fastMNN(input.mnn[[1]], input.mnn[[2]], input.mnn[[3]],
                       d = pc.num, k = 20)
}
if (length(input.mnn) == 4) {
    out.mnn <- fastMNN(input.mnn[[1]], input.mnn[[2]], input.mnn[[3]],
                       input.mnn[[4]], d = pc.num, k = 20)
}
if (length(input.mnn) == 5) {
    out.mnn <- fastMNN(input.mnn[[1]], input.mnn[[2]], input.mnn[[3]],
                       input.mnn[[4]], input.mnn[[5]], d = pc.num, k = 20)
}
if (length(input.mnn) == 6) {
    out.mnn <- fastMNN(input.mnn[[1]], input.mnn[[2]], input.mnn[[3]],
                       input.mnn[[4]], input.mnn[[5]], input.mnn[[6]], 
                       d = pc.num, k = 20)
}

time6 <- Sys.time()

mtx.mnn <- out.mnn[["corrected"]]
dimnames(mtx.mnn)[[1]] <- input.cells
dimnames(mtx.mnn)[[2]] <- paste('MNN', 1:(pc.num), sep = '-')

object.mnn <- CreateSeuratObject(counts = mtx.norm)
object.mnn@reductions$pca <- CreateDimReducObject(
    embeddings = mtx.mnn, key = 'PC_', assay = 'RNA'
)
object.mnn <- RunUMAP(
    object.mnn, reduction = 'pca', dims = 1:pc.num, n.components = 2, 
    verbose = F)

if (sum(input.cells == cell.names) == length(cell.names)) {
    object.mnn@meta.data <- object@meta.data
} else {
    print("Error: MNN Cell order not match")
}
plot.mnn <- DimPlot(
    object.mnn, reductions = 'umap', 
    group.by = "group", shape.by = 'batch', pt.size = 1.5
)
ggsave(plot = plot.mnn, filename = 'UMAP_mnn.png', path = path,
       units = 'cm', height = 15, width = 25)
# result.mnn <- evaluate.two.dims(object.mnn)
result.mnn <- evaluate.sil(object.mnn)
df.result <- rbind(df.result, result.mnn)
timing <- c(timing, difftime(time6, time5, units = 'secs'))
label <- c(label, 'MNN')


### BBKNN
time7 <- Sys.time()
library(reticulate)
use_python("~/tools/anaconda3/bin/python3")
source_python('bbknnpy.py')
pca.input <- object@reductions$pca@cell.embeddings
umap.bbknn <- bbknn_py(pca.input, batches, pc.num)
time8 <- Sys.time()
# use_python("~/tools/anaconda3/bin/python3")
# 
# anndata = import("anndata", convert = FALSE)
# bbknn = import("bbknn", convert = FALSE)
# sc = import("scanpy.api", convert = FALSE)
# 
# pca.input <- object@reductions$pca@cell.embeddings
# adata <- anndata$AnnData(X = pca.input, obs = batches)
# sc$tl$pca(adata)
# adata$obsm$X_pca <- r_to_py(pca.input)
# sc$pp$bbknn(adata, batch_key = 0)
# sc$tl$umap(adata)
# umap.bbknn <- py_to_r(adata$obsm['X_umap'])
dimnames(umap.bbknn)[[1]] <- dimnames(mtx.norm)[[2]]

object.bbknn <- CreateSeuratObject(counts = mtx.norm)
object.bbknn@reductions$umap <- CreateDimReducObject(
    embeddings = umap.bbknn, key = 'UMAP_', assay = 'RNA'
)
object.bbknn@meta.data <- object@meta.data
plot.bbknn <- DimPlot(
    object.bbknn, reductions = 'umap', 
    group.by = "group", shape.by = 'batch', pt.size = 1.5
)
ggsave(plot = plot.bbknn, filename = 'UMAP_bbknn.png', path = path,
       units = 'cm', height = 15, width = 25)
# result.bbknn <- evaluate.two.dims(object.bbknn)
result.bbknn <- evaluate.sil(object.bbknn)
df.result <- rbind(df.result, result.bbknn)
timing <- c(timing, difftime(time8, time7, units = 'secs'))
label <- c(label, 'BBKNN')


### Seurat
time9 <- Sys.time()
batch.list <- list()
cell.names <- dimnames(object@assays$RNA@counts)[[2]]
for (batch in unique(batches)) {
    sub.cells <- cell.names[object@meta.data$batch == batch]
    object.sub <- CreateSeuratObject(counts = mtx.norm[,sub.cells])
    object.sub@assays$RNA <- CreateAssayObject(data = mtx.norm[,sub.cells])
    object.sub@commands$NormalizeData.RNA <- object@commands$NormalizeData.RNA
    object.sub@commands$FindVariableFeatures.RNA <- 
        object@commands$FindVariableFeatures.RNA
    object.sub@assays$RNA@var.features <- object@assays$RNA@var.features
    batch.list[[batch]] <- object.sub
}
anchors <- FindIntegrationAnchors(
    object.list = batch.list, dims = 1:pc.num, k.filter = 50)
object.integrated <- IntegrateData(anchorset = anchors, dims = 1:pc.num)
time10 <- Sys.time()
object.integrated <- ScaleData(
    object.integrated, assay = 'integrated', 
    feature = object@assays$RNA@var.features, 
    verbose = F)
object.integrated <- RunPCA(object.integrated, npcs = pc.num, verbose = F)
object.integrated <- RunUMAP(
    object.integrated, reduction = 'pca', dims = 1:pc.num)
if (sum(dimnames(object.integrated@assays$integrated@data)[[2]] == 
        cell.names) == length(cell.names)) {
    object.integrated@meta.data <- object@meta.data
} else {
    print("Error: Seurat Cell order not match")
}
plot.integrated <- DimPlot(
    object.integrated, reductions = 'umap', 
    group.by = "group", shape.by = 'batch', pt.size = 1.5
)
ggsave(plot = plot.integrated, filename = 'UMAP_integrated.png', path = path,
       units = 'cm', height = 15, width = 25)
# result.integrated <- evaluate.two.dims(object.integrated)
result.integrated <- evaluate.sil(object.integrated)
df.result <- rbind(df.result, result.integrated)
timing <- c(timing, difftime(time10, time9, units = 'secs'))
label <- c(label, 'Seurat')


# result
df.result$label <- label
df.result$timing <- timing
df.result$dataset <- folder.data
path.write <- paste(path, paste0(folder.data, '.txt'), sep = '/')
write.table(df.result, file = path.write, quote = F, sep = '\t')

# plot
df.plot <- data.frame()
for (i in 1:dim(df.result)[1]) {
    df.plot <- rbind(
        df.plot,
        data.frame(score = df.result[i, 'cell.distance'],
                   method = df.result[i, 'label'],
                   evaluator = 'Cell Distance'))
        df.plot <- rbind(
            df.plot,
            data.frame(score = df.result[i, 'batch.effect.factor'],
                       method = df.result[i, 'label'],
                       evaluator = 'Batch Mixing Degree'))
        df.plot <- rbind(
            df.plot,
            data.frame(score = df.result[i, 'harmonic.mean'],
                       method = df.result[i, 'label'],
                       evaluator = 'Harmonic Mean'))
        if (i == 1) {
            df.plot <- rbind(
                df.plot,
                data.frame(score = 0,
                           method = df.result[i, 'label'],
                           evaluator = 'Timing'))
        } else {
            df.plot <- rbind(
                df.plot,
                data.frame(score = df.result[i, 'timing'],
                           method = df.result[i, 'label'],
                           evaluator = 'Timing'))
        }
}

bar.plot <-
    ggplot(df.plot[df.plot$evaluator != 'Timing', ],
           aes(x = evaluator, y = score, fill = method)) +
    geom_bar(position = 'dodge', stat = 'identity') +
    labs(title = "Evaluate multiple methods of removing batch effect",
         y = 'Score', x = 'Evaluator') +
    scale_fill_discrete(name = 'Method')

ggsave(
    plot = bar.plot, path = path, filename = "test_methods.png",
    units = 'cm', width = 15, height = 10)


time.plot <-
    ggplot(df.plot[df.plot$evaluator == 'Timing', ],
           aes(x = evaluator, y = score, fill = method)) +
    geom_bar(position = 'dodge', stat = 'identity') +
    labs(title = "Run time of methods",
         y = 'Score', x = 'Evaluator') +
    scale_fill_discrete(name = 'Method')

ggsave(
    plot = time.plot, path = path, filename = "calc_timing.png",
    units = 'cm', width = 8, height = 10)




