tran_mat <- function(path){
    library(ggplot2)
    library(cowplot)
    library(Matrix)
    library(Seurat)
    pbmc.data <- Read10X(data.dir = path)
    mnt <- as.matrix(pbmc.data)
    df <- as.data.frame(mnt)
    outfile <- paste(path, 'res_mnt_pre.txt', sep = '/')
    write.table(df, outfile, sep = '\t', quote = F)
}
args <- commandArgs(T)
tran_mat(args[1])
