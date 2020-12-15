# PcoA from distance
library(amplicon)
library(cluster)
library(ggplot2)

# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.out.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)
df.meta$Dose <- as.factor(df.meta$Dose)

# distance matrix
setwd('/home/drizzle_zhang/microbiome/result/5.Beta_Diversity/Distance')
type.distance <- 'bray_curtis'
# type.distance <- 'euclidean'
# type.distance <- 'unweighted_unifrac'
gender <- 'female'
use.dim <- 5
vec.dose <- c(0, 1, 2, 3)
# vec.dose <- c(0, 3)
file.distance <- paste0('./', type.distance, '_otu_table_even.txt')
mat.distance <- read.table(file.distance, sep = '\t', header = T, row.names = 1)

# time series
path.plot <- paste0(
    '/home/drizzle_zhang/microbiome/result/5.Beta_Diversity/PCoA/', 
    type.distance, '_', gender)
if (!file.exists(path.plot)) {
    dir.create(path.plot)
}
series.time <- unique(df.meta$Time)
vector.sil <- c()
for (sub.time in series.time) {
    # select meta
    sel.meta <- df.meta
    sel.meta <- df.meta[df.meta$Time == sub.time,]
    sel.meta <- sel.meta[sel.meta$Dose %in% vec.dose,]
    sel.meta <- sel.meta[sel.meta$Gender == gender,]
    row.names(sel.meta) <- sel.meta$SampleName
    
    # select sample
    use.sample <- sel.meta$SampleName
    mat.select <- mat.distance[use.sample, use.sample]
    
    # PcoA
    plot.beta <- 
        beta_pcoa(mat.select, sel.meta, groupID = 'Dose') + 
        geom_text(aes(label = use.sample))
    ggsave(plot = plot.beta, path = path.plot, 
           filename = paste0(
           paste0(as.character(vec.dose), collapse = ''), '_', sub.time, '.png'))
    
    # calcaulate distance in plot
    pcoa = cmdscale(mat.select, k = use.dim, eig = T)
    dissimilar.dist <- dist(pcoa$points)
    sil.out <- silhouette(as.numeric(as.factor(sel.meta$Group)), dissimilar.dist)
    sil.group <- abs(summary(sil.out)$avg.width)
    vector.sil <- c(vector.sil, sil.group)
    
}

plot.fit <- data.frame(Distance = vector.sil[2:length(vector.sil)],
                       Time = as.numeric(series.time)[2:length(series.time)])
ggplot(data = plot.fit, aes(x = Time, y = Distance)) +
    geom_smooth(method = lm, formula = y ~ poly(x, 3)) +
    geom_point() +
    geom_hline(yintercept = vector.sil[1], color = 'gray', linetype = 5) +
    theme(panel.background = element_rect(color = 'gray',
                                          fill = 'transparent')) +
    ylab('Silhouette Distance')

df.save <- data.frame(Distance = vector.sil,
                      Time = series.time,
                      Gender = rep(gender, length(series.time)))
file.save <- paste0(
    path.plot, paste0('/diastance_', as.character(use.dim), '_',
                      paste0(as.character(vec.dose), collapse = ''), '.txt'))
write.table(df.save, file.save, sep = '\t', quote = F, row.names = F)
