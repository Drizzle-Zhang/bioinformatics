# if (!requireNamespace("devtools", quietly=TRUE))    
#     install.packages("devtools")
# library(devtools)
# if (!requireNamespace("amplicon", quietly=TRUE))    
#     install_github("microbiota/amplicon")
# suppressWarnings(suppressMessages(library(amplicon)))

library(ggplot2)
library(gg3D)

# input file 
file_unweighted_unifrac <- 
    '/home/drizzle_zhang/microbiome/result/5.Beta_Diversity/PCoA/pcoa_unweighted_unifrac_otu_table_even.txt'
file_bray_curtis <- 
    '/home/drizzle_zhang/microbiome/result/5.Beta_Diversity/PCoA/pcoa_bray_curtis_otu_table_even.txt'

mtx_unweighted_unifrac <- read.table(file_unweighted_unifrac, header = T, 
                                     row.names = 1, sep = '\t')
mtx_bray_curtis <- read.table(file_bray_curtis, header = T, 
                              row.names = 1, sep = '\t')

# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.out.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)

sel.meta <- df.meta
# sel.meta <- df.meta[df.meta$Time == 'A',]
sel.meta <- df.meta[df.meta$Dose %in% c(0, 3),]

mtx.in <- mtx_bray_curtis
use.pc <- 3
col.name <- paste0('X', 1:use.pc)
mat.plot <- mtx.in[sel.meta$Sample, col.name]
names(mat.plot) <- paste0('PC', 1:use.pc)
mat.plot <- cbind(mat.plot, sel.meta)
### 2D plot
# plot.unweighted_unifrac <- 
ggplot(mat.plot, aes(x = PC1, y = PC3, color = Group)) + 
    geom_point(size = 2)
ggplot(mat.plot, aes(x = PC1, y = PC3, color = Time)) + 
    geom_point(size = 2)

# theta=60
# phi=30
# # plot.unweighted_unifrac <- 
# ggplot(mat.plot, aes(x = PC1, y = PC2, z = PC3, color = Group)) + 
#     axes_3D(theta=theta, phi=phi) + stat_3D(theta=theta, phi=phi) + 
#     labs_3D(theta=theta, phi=phi)