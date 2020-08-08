# PcoA from distance
library(amplicon)
library(cluster)

# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.out.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)
df.meta$Dose <- as.factor(df.meta$Dose)

# distance matrix
setwd('/home/drizzle_zhang/microbiome/result/5.Beta_Diversity/Distance')
type.distance <- 'bray_curtis'
# type.distance <- 'weighted_unifrac'
gender <- 'male'
vec.dose <- c(0, 1, 2, 3)
file.distance <- paste0('./', type.distance, '_otu_table_even.txt')
mat.distance <- read.table(file.distance, sep = '\t', header = T, row.names = 1)

# time series
path.plot <- paste0(
    '/home/drizzle_zhang/microbiome/result/5.Beta_Diversity/PCoA/', 
    type.distance, '_gender_diff')
series.time <- unique(df.meta$Time)
vector.sil <- c()
for (sub.time in series.time) {
    # select meta
    sel.meta <- df.meta
    sel.meta <- df.meta[df.meta$Time == sub.time,]
    sel.meta <- sel.meta[sel.meta$Dose %in% vec.dose,]
    # sel.meta <- sel.meta[sel.meta$Gender == gender,]
    row.names(sel.meta) <- sel.meta$Sample
    
    # select sample
    use.sample <- sel.meta$Sample
    mat.select <- mat.distance[use.sample, use.sample]
    
    # PcoA
    pcoa <- cmdscale(mat.select, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
    points <- as.data.frame(pcoa$points) # get coordinate string, format to dataframme
    eig <- pcoa$eig
    colnames(points) <- c("x", "y", "z")
    points <- cbind(points, sel.meta[rownames(points),])
    # plot.beta <- beta_pcoa(mat.select, sel.meta, groupID = 'Gender')
    # plot.beta <- 
    #     ggplot(points, aes(x = x, y = y, shape = Group, color = Gender)) +
    #     labs(x=paste("PCo 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
    #          y=paste("PCo 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
    #          color='Gender', shape='Group') +
    #     geom_point(alpha=.7, size=2) + theme_classic() + theme(text=element_text(family="sans", size=7)) + 
    #     stat_ellipse(level = 0.68)
    plot.beta <- 
        ggplot(points, aes(x = x, y = y, color = Gender)) +
        labs(x=paste("PCo 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
             y=paste("PCo 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
             color='Gender') +
        geom_point(alpha=.7, size=2) + theme_classic() + theme(text=element_text(family="sans", size=7)) + 
        stat_ellipse(level = 0.68)
    
    ggsave(plot = plot.beta, path = path.plot, 
           filename = paste0(
               sub.time, '_', paste0(as.character(vec.dose), collapse = ''), 
               '.png'))
    
    # calcaulate distance in plot
    # pcoa = cmdscale(mat.select, k = 3, eig = T)
    # dissimilar.dist <- dist(pcoa$points)
    # sil.out <- silhouette(as.numeric(as.factor(sel.meta$Dose)), dissimilar.dist)
    # sil.group <- abs(summary(sil.out)$avg.width)
    # vector.sil <- c(vector.sil, sil.group)
    
}

plot.fit <- data.frame(Distance = vector.sil[2:length(vector.sil)],
                       Time = as.numeric(as.factor(series.time))[2:length(series.time)])
ggplot(data = plot.fit, aes(x = Time, y = Distance)) + 
    geom_smooth(method = lm, formula = y ~ poly(x, 3)) + 
    geom_point() + 
    geom_hline(yintercept = vector.sil[1], color = 'gray', linetype = 5) + 
    theme(panel.background = element_rect(color = 'gray', 
                                          fill = 'transparent')) + 
    ylab('Silhouette Distance')


