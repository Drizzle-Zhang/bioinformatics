setwd("~/my_git/bioinformatics/BEE/test_methods")
source('functions_evaluator.R')
library(ggplot2)

# read data
path <- '/home/drizzle_zhang/Desktop/single_cell/BEE/real_data/10X_balance'
file.mtx <- paste(path, 'matrix.txt', sep = '/')
file.meta <-  paste(path, 'meta.txt', sep = '/')
file.info <- 
    '/home/drizzle_zhang/Desktop/single_cell/BEE/real_data/lab_info.txt'
df.mtx <- read.table(file.mtx, sep = '\t', header = T, stringsAsFactors = F)
row.names(df.mtx) <- df.mtx[, 'Gene']
df.mtx['Gene'] <- NULL
df.mtx <- na.omit(df.mtx)
df.meta <- read.table(file.meta, sep = '\t', header = T, stringsAsFactors = F)
groups <- df.meta$cell
df.info <- read.delim(file.info, stringsAsFactors = F)
batches <- c()
for (batch in df.meta$folder) {
    batches <- c(batches, df.info[df.info$folder == batch, 'author'])
}
pc.num <- 50


label <- c()
timeing <- c()
df.result <- data.frame()
# origin
time1 <- Sys.time()
object <- preprocess.data(df.mtx, batches, groups, num.pc = pc.num)
plot.origin <- DimPlot(
    object, reductions = 'umap', 
    group.by = "group", shape.by = 'batch', pt.size = 1.5
)
ggsave(plot = plot.origin, filename = 'UMAP_origin.png', path = path,
       units = 'cm', height = 15, width = 25)
result.origin <- evaluate.two.dims(object)
df.result <- rbind(df.result, result.origin)
time2 <- Sys.time()
timeing <- c(timeing, time2 - time1)
label <- c(label, 'Origin')

# normalized data
mtx.norm <- as.matrix(object@assays$RNA@data)[object@assays$RNA@var.features,]

# combat
time3 <- Sys.time()
library(sva)
mtx.combat <- ComBat(dat = mtx.norm, batch = object@meta.data$batch)
time4 <- Sys.time()
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
result.combat <- evaluate.two.dims(object.combat)
df.result <- rbind(df.result, result.combat)
timeing <- c(timeing, time4 - time3)
label <- c(label, 'ComBat')

# BEER
time5 <- Sys.time()
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
mybeer <- BEER(mtx.norm, object@meta.data$batch, GNUM = 30, PCNUM = pc.num, 
               ROUND = 1, GN = 2000, SEED = 1, COMBAT = TRUE)
time6 <- Sys.time()
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
result.beer <- evaluate.two.dims(object.beer)
df.result <- rbind(df.result, result.beer)
timeing <- c(timeing, time6 - time5)
label <- c(label, 'BEER')

# MNN
time7 <- Sys.time()
library(scran)
# library(psych)
# library(batchelor)
# library(SingleCellExperiment)
# sce <- SingleCellExperiment(assays = list(logcounts = mtx.norm))
# out.mnn <- fastMNN(sce, batch = as.factor(object@meta.data$batch), 
#                     d = pc.num, k = 20)
input.mnn <- list()
cell.names <- dimnames(object@assays$RNA@counts)[[2]]
for (lab in unique(object@meta.data$batch)) {
    lab.cells <- cell.names[object@meta.data$batch == lab]
    input.mnn[[lab]] <- mtx.norm[,lab.cells]
}
out.mnn <- fastMNN(input.mnn[[1]], input.mnn[[2]], input.mnn[[3]], 
                   input.mnn[[4]], input.mnn[[5]], d = pc.num, k = 20)
time8 <- Sys.time()

mtx.mnn <- out.mnn[["corrected"]]
dimnames(mtx.mnn)[[1]] <- dimnames(mtx.norm)[[2]]
dimnames(mtx.mnn)[[2]] <- paste('MNN', 1:(pc.num), sep = '-')

# object.mnn <- CreateSeuratObject(counts = t(mtx.mnn))
# object.mnn@assays$RNA <- CreateAssayObject(data = t(mtx.mnn))
# object.mnn@commands$NormalizeData.RNA <- object@commands$NormalizeData.RNA
# object.mnn@commands$FindVariableFeatures.RNA <-
#     object@commands$FindVariableFeatures.RNA
# object.mnn@assays$RNA@var.features <- dimnames(mtx.mnn)[[2]]
# object.mnn <-
#     preprocess.data(object.mnn, batches, groups, num.pc = pc.num)
object.mnn <- CreateSeuratObject(counts = mtx.norm)
object.mnn@reductions$pca <- CreateDimReducObject(
    embeddings = mtx.mnn, key = 'PC_', assay = 'RNA'
)
object.mnn <- RunUMAP(
    object.mnn, reduction = 'pca', dims = 1:pc.num, n.components = 2, 
    verbose = F)
object.mnn@meta.data <- object@meta.data
plot.mnn <- DimPlot(
    object.mnn, reductions = 'umap', 
    group.by = "group", shape.by = 'batch', pt.size = 1.5
)
ggsave(plot = plot.mnn, filename = 'UMAP_mnn.png', path = path,
       units = 'cm', height = 15, width = 25)
result.mnn <- evaluate.two.dims(object.mnn)
df.result <- rbind(df.result, result.mnn)
timeing <- c(timeing, time8 - time7)
label <- c(label, 'MNN')


# result
df.result$label <- label
df.result$timeing <- timeing

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
                data.frame(score = df.result[i, 'timeing'],
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




