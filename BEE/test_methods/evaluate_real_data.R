source('functions_evaluator.R')

# read data
path <- '/home/drizzle_zhang/Desktop/single_cell/BEE/real_data/drop_balance'
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

object <- preprocess.data(df.mtx, batches, groups)
result <- evaluate.two.dims(object)
DimPlot(
    object, reductions = 'umap', 
    group.by = "group", shape.by = 'batch', pt.size = 1.5
)

select_folders <- c('Blanca Pijuan-Sala', 'Ximena Ibarra-Soria')
cells <- dimnames(object@assays$RNA@counts)[[2]]
cells2 <- cells[object@meta.data$batch %in% select_folders]
object2 <- subset(object, cells = cells2)
DimPlot(
    object2, reductions = 'umap', 
    group.by = "group", shape.by = 'batch', pt.size = 2
)
result2 <- evaluate.two.dims(object2)

