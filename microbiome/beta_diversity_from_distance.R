# PcoA from distance
library(amplicon)

# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.out.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)

# distance matrix
setwd('/home/drizzle_zhang/microbiome/result/5.Beta_Diversity/Distance')
type.distance <- 'bray_curtis'
file.distance <- paste0('./', type.distance, '_otu_table_even.txt')
mat.distance <- read.table(file.distance, sep = '\t', header = T, row.names = 1)

# time series
path.plot <- paste0(
    '/home/drizzle_zhang/microbiome/result/5.Beta_Diversity/PCoA/', type.distance)
series.time <- unique(df.meta$Time)
for (sub.time in series.time) {
    # select meta
    sel.meta <- df.meta
    sel.meta <- df.meta[df.meta$Time == sub.time,]
    sel.meta <- sel.meta[sel.meta$Dose %in% c(0, 3),]
    row.names(sel.meta) <- sel.meta$Sample
    
    # select sample
    use.sample <- sel.meta$Sample
    mat.select <- mat.distance[use.sample, use.sample]
    
    # PcoA
    plot.beta <- beta_pcoa(mat.select, sel.meta)
    ggsave(plot = plot.beta, path = path.plot, 
           filename = paste0(sub.time, '_0vs3.png'))
    
}
