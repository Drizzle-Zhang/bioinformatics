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
df.meta$Dose <- as.factor(df.meta$Dose)

sel.meta <- df.meta
# sel.meta <- df.meta[df.meta$Time == 'A',]
# sel.meta <- df.meta[df.meta$Dose %in% c(1, 2, 3),]

mtx.in <- mtx_bray_curtis
use.pc <- 10
col.name <- paste0('X', 1:use.pc)
mat.plot <- mtx.in[sel.meta$Sample, col.name]
names(mat.plot) <- paste0('PCo', 1:use.pc)
mat.plot <- cbind(mat.plot, sel.meta)
mat.plot.dose <- mat.plot[mat.plot$Dose %in% c(1, 2, 3),]
### 2D plot
# plot.unweighted_unifrac <- 
ggplot(mat.plot, aes(x = PCo1, y = PCo5, color = Group)) + 
    geom_point(size = 2)
ggplot(mat.plot, aes(x = PCo3, y = PCo5, color = Group)) + 
    geom_point(size = 2)
ggplot(mat.plot.dose, aes(x = PCo1, y = PCo5, color = Dose)) +
    geom_point(size = 2)
ggplot(mat.plot.dose, aes(x = PCo5, y = PCo10, color = Dose)) +
    geom_point(size = 2)
ggplot(mat.plot, aes(x = PCo1, y = PCo5, color = Time)) + 
    geom_point(size = 2)
ggplot(mat.plot, aes(x = PCo5, y = PCo6, color = Gender)) +
    geom_point(size = 2)
ggplot(mat.plot, aes(x = PCo1, y = PCo5, color = Gender)) + 
    geom_point(size = 2)

df.pc <- data.frame()
for (sub.pc in paste0('PCo', 1:use.pc)) {
    wilcox.Group <- wilcox.test(
        formula = as.formula(paste0(sub.pc, ' ~ Group')), data = mat.plot)
    wilcox.Gender <- wilcox.test(
        formula = as.formula(paste0(sub.pc, ' ~ Gender')), data = mat.plot)
    kruskal.Dose <- kruskal.test(
        formula = as.formula(paste0(sub.pc, ' ~ Dose')), data = mat.plot.dose)
    df.pc <- rbind(df.pc, 
                   data.frame(Pco = sub.pc, 
                              Pval.Group = wilcox.Group$p.value,
                              Pval.Gender = wilcox.Gender$p.value,
                              Pval.Dose = kruskal.Dose$p.value))
}

# theta=60
# phi=30
# # plot.unweighted_unifrac <- 
# ggplot(mat.plot, aes(x = PC1, y = PC2, z = PC3, color = Group)) + 
#     axes_3D(theta=theta, phi=phi) + stat_3D(theta=theta, phi=phi) + 
#     labs_3D(theta=theta, phi=phi)