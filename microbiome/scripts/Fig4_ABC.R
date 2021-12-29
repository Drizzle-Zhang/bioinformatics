library(ggplot2)
library(gg3D)
library(amplicon)

# input file 
file_bray_curtis <- 
    '/home/drizzle_zhang/microbiome/result/5.Beta_Diversity/PCoA/pcoa_bray_curtis_otu_table_even.txt'

# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)
df.meta$Dose <- as.factor(df.meta$Dose)

# decomposition from distance matrix
type.distance <- 'bray_curtis'
file.distance <- 
    paste0('/home/drizzle_zhang/microbiome/result/5.Beta_Diversity/Distance/', 
           type.distance, '_otu_table_even.txt')
mat.distance <- read.table(file.distance, sep = '\t', header = T, row.names = 1)
pcoa = cmdscale(mat.distance, k = 10, eig = T)
points <- as.data.frame(pcoa$points)
eig <- pcoa$eig
mtx_bray_curtis <- points

mtx.in <- mtx_bray_curtis
use.pc <- 10
col.name <- paste0('V', 1:use.pc)
mat.plot <- mtx.in[df.meta$SampleName, col.name]
names(mat.plot) <- paste0('PCo', 1:use.pc)
mat.plot <- cbind(mat.plot, df.meta)
mat.plot.dose <- mat.plot[mat.plot$Dose %in% c(1, 2, 3),]
plot.unweighted_unifrac <-
ggplot(mat.plot, aes(x = PCo1, y = PCo4, color = Group)) + 
    geom_point(size = 2)
ggplot(mat.plot, aes(x = PCo4, y = PCo6, color = Group)) + 
    geom_point(size = 2)
ggplot(mat.plot.dose, aes(x = PCo1, y = PCo5, color = Dose)) +
    geom_point(size = 2)
ggplot(mat.plot.dose, aes(x = PCo5, y = PCo10, color = Dose)) +
    geom_point(size = 2)
ggplot(mat.plot, aes(x = PCo1, y = PCo2, color = as.factor(Time))) + 
    geom_point(size = 2)
ggplot(mat.plot, aes(x = PCo1, y = PCo4, color = as.factor(Time))) + 
    geom_point(size = 2)
ggplot(mat.plot, aes(x = PCo5, y = PCo6, color = Gender)) +
    geom_point(size = 2)
ggplot(mat.plot, aes(x = PCo1, y = PCo5, color = Gender)) + 
    geom_point(size = 2)
ggplot(mat.plot, aes(x = PCo1, y = PCo3, color = as.factor(Group))) + 
    geom_point(size = 2)

df.pc <- data.frame()
for (sub.pc in paste0('PCo', 1:use.pc)) {
    wilcox.Group <- wilcox.test(
        formula = as.formula(paste0(sub.pc, ' ~ Group')), data = mat.plot)
    wilcox.Gender <- wilcox.test(
        formula = as.formula(paste0(sub.pc, ' ~ Gender')), data = mat.plot)
    kruskal.Dose <- kruskal.test(
        formula = as.formula(paste0(sub.pc, ' ~ Dose')), data = mat.plot.dose)
    kruskal.Time <- kruskal.test(
        formula = as.formula(paste0(sub.pc, ' ~ Time')), data = mat.plot)
    df.pc <- rbind(df.pc, 
                   data.frame(Pco = sub.pc, 
                              Pval.Group = wilcox.Group$p.value,
                              Pval.Gender = wilcox.Gender$p.value,
                              Pval.Dose = kruskal.Dose$p.value,
                              Pval.Time = kruskal.Time$p.value))
}
# output fig
path.plot <- '/home/drizzle_zhang/microbiome/result/Figs/'

# time
getPalette = colorRampPalette(rev(brewer.pal(10, "Spectral")))
vec.color <- c(getPalette(length(series.time)))
plot.time <- 
    ggplot(mat.plot, aes(x = PCo1, y = PCo4, color = as.factor(Time))) + 
    geom_point(size = 1.5) + 
    scale_color_manual(values = rev(vec.color)) + 
    labs(x = paste0('PCo 1 (', format(100 * eig[1]/sum(eig), digits = 4), '%)'),
         y = paste0('PCo 4 (', format(100 * eig[4]/sum(eig), digits = 4), '%)'),
         color = 'Time (weeks)') + 
    theme_bw() +
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'bottom')
ggsave(path = path.plot,
       filename = paste0(type.distance, '_Time.png'),
       plot = plot.time, height = 10, width = 10, units = 'cm')

# group
plot.group <- 
    ggplot(mat.plot, aes(x = PCo1, y = PCo4, color = Group)) + 
    geom_point(size = 1.5) + 
    scale_color_manual(breaks = c('Control', 'Treat'),
                       values = c("#5AB4AC", "#BC80BD"),
                       labels = c('Control', 'Treatment')) +
    labs(x = paste0('PCo 1 (', format(100 * eig[1]/sum(eig), digits = 4), '%)'),
         y = paste0('PCo 4 (', format(100 * eig[4]/sum(eig), digits = 4), '%)'),
         color = 'Group') + 
    stat_ellipse(level = 0.95) +
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.key = element_blank(),
          legend.position = 'bottom')
ggsave(path = path.plot,
       filename = paste0(type.distance, '_Group.png'),
       plot = plot.group, height = 8.8, width = 10, units = 'cm')

# gender
plot.gender <- 
    ggplot(mat.plot, aes(x = PCo5, y = PCo8, color = Gender)) + 
    geom_point(size = 1.5) + 
    scale_color_manual(breaks = c('male', 'female'),
                       values = c("#80B1D3", "#FFA07A"),
                       labels = c('Male', 'Female')) +
    labs(x = paste0('PCo 5 (', format(100 * eig[1]/sum(eig), digits = 4), '%)'),
         y = paste0('PCo 8 (', format(100 * eig[8]/sum(eig), digits = 4), '%)')) + 
    stat_ellipse(level = 0.95) +
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.key = element_blank(),
          legend.position = 'bottom')
ggsave(path = path.plot,
       filename = paste0(type.distance, '_Gender.png'),
       plot = plot.gender, height = 8.8, width = 10, units = 'cm')

# dose
plot.time <- ggplot(mat.plot, aes(x = PCo3, y = PCo5, color = Dose)) + 
    geom_point(size = 2) + 
    labs(x = paste0('PCo 3 (', format(100 * eig[3]/sum(eig), digits = 4), '%)'),
         y = paste0('PCo 5 (', format(100 * eig[5]/sum(eig), digits = 4), '%)')) + 
    theme(panel.background = element_rect(color = 'gray', size = 2,
                                          fill = 'transparent'),) + 
    stat_ellipse(level = 0.95)
ggsave(path = '/home/drizzle_zhang/microbiome/result/5.Beta_Diversity/PCoA',
       filename = paste0(type.distance, '_Dose.png'),
       plot = plot.time)

# group
sub.time = '1'
vec.dose <- c(0, 3)
gender <- 'male'
df.meta <- df.meta
df.meta <- df.meta[df.meta$Time == sub.time,]
df.meta <- df.meta[df.meta$Dose %in% vec.dose,]
df.meta <- df.meta[df.meta$Gender == gender,]
row.names(df.meta) <- df.meta$SampleName
# select sample
use.sample <- df.meta$SampleName
mat.select <- mat.distance[use.sample, use.sample]
plot.beta <- 
    beta_pcoa(mat.select, df.meta, groupID = 'Group') + 
    theme_bw() + 
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'bottom')
    
ggsave(path = path.plot,
       filename = paste0('example_Group.png'),
       plot = plot.beta, height = 9, width = 8, units = 'cm')

