library(Seurat)
library(ggplot2)

meta.in <- '/local/zy/PEI/data/DHS/GRCh38tohg19_cluster/metadata.simple.tsv'
file.in <- '/local/zy/PEI/data/DHS/GRCh38tohg19_cluster/accession_matrix.txt'
path.out <- '/local/zy/PEI/data/DHS/GRCh38tohg19_cluster/'

df.lable.peak <- read.delim(file.in, sep = '\t', header = T,
                            row.names = 'peak_id', stringsAsFactors = F)
df.meta <- read.delim(meta.in, sep = '\t', header = T, stringsAsFactors = F)
df.meta <- df.meta[df.meta$File.accession %in% names(df.lable.peak),]

object <- CreateSeuratObject(mtx.count)
object@meta.data$organ <- df.meta$Biosample.organ
object@meta.data$lifestage <- df.meta$Biosample.life.stage
object <- NormalizeData(object, verbose = F)


