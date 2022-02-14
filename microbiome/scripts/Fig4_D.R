# PcoA from distance
library(amplicon)
library(cluster)
library(ggplot2)

# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)
df.meta$Dose <- as.factor(df.meta$Dose)

# distance matrix
setwd('/home/drizzle_zhang/microbiome/result/5.Beta_Diversity/Distance')
type.distance <- 'bray_curtis'
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

# data combination
type.distance <- 'bray_curtis'
use.dim <- '5'

vec.gender <- c('male', 'female')
vec.linetype <- c(1, 2)
vec.point <- c(16, 21)
vec.str.dose <- c('01', '02', '03')
vec.color <- c('firebrick', '#33A02C', '#1F78B4')

df.combine <- data.frame(stringsAsFactors = F)
for (i in 1:length(vec.gender)) {
    gender <- vec.gender[i]
    for (j in 1:length(vec.str.dose)) {
        str.dose <- vec.str.dose[j]
        path.plot <- paste0(
            '/home/drizzle_zhang/microbiome/result/5.Beta_Diversity/PCoA/', 
            type.distance, '_', gender)
        file.male <- paste0(path.plot, paste0('/diastance_', use.dim, '_', str.dose, '.txt'))
        df.sub <- read.delim(file.male, stringsAsFactors = F)
        df.sub$Gender[df.sub$Gender == 'female'] <- 'Female'
        df.sub$Gender[df.sub$Gender == 'male'] <- 'Male'
        df.sub$Dose <- rep(str.dose, nrow(df.sub))
        df.sub$LineType <- rep(vec.linetype[i], nrow(df.sub))
        df.sub$PointShape <- rep(vec.point[i], nrow(df.sub))
        df.sub$Color <- rep(vec.color[j], nrow(df.sub))
        df.combine <- rbind(df.combine, df.sub)
    }
}

# baseline
linewidth.base <- 0.2
df.baseline <- df.combine[df.combine$Time == -1,]
baseline.male.01 <- df.baseline[(df.baseline$Gender == 'Male') & (df.baseline$Dose == '01'), 'Distance']
baseline.male.02 <- df.baseline[(df.baseline$Gender == 'Male') & (df.baseline$Dose == '02'), 'Distance']
baseline.male.03 <- df.baseline[(df.baseline$Gender == 'Male') & (df.baseline$Dose == '03'), 'Distance']
baseline.female.01 <- df.baseline[(df.baseline$Gender == 'Female') & (df.baseline$Dose == '01'), 'Distance']
baseline.female.02 <- df.baseline[(df.baseline$Gender == 'Female') & (df.baseline$Dose == '02'), 'Distance']
baseline.female.03 <- df.baseline[(df.baseline$Gender == 'Female') & (df.baseline$Dose == '03'), 'Distance']
df.plot <- df.combine[df.combine$Time != -1, ]

# plot
plot.fit <- 
    ggplot(data = df.plot, aes(x = Time, y = Distance, color = Dose, 
                               shape = Gender, linetype = Gender)) +
    geom_smooth(method = lm, formula = y ~ poly(x, 3), se = F) +
    geom_point() +
    scale_linetype_manual(breaks = c('Male', 'Female'), values = vec.linetype) + 
    scale_color_manual(breaks = vec.str.dose, values = vec.color,
                       labels = c('0.0 vs. 0.375', '0.0 vs. 1.5',
                                  '0.0 vs. 6.0')) + 
    labs(color = 'Dosage (mg/kg IAA)') + 
    scale_shape_manual(breaks = c('Male', 'Female'), values = vec.point) + 
    geom_hline(yintercept = baseline.male.01, color = 'firebrick', 
               linetype = 1, size = linewidth.base) +
    geom_hline(yintercept = baseline.male.02, color = '#33A02C', 
               linetype = 1, size = linewidth.base) +
    geom_hline(yintercept = baseline.male.03, color = '#1F78B4', 
               linetype = 1, size = linewidth.base) +
    geom_hline(yintercept = baseline.female.01, color = 'firebrick', 
               linetype = 5, size = linewidth.base) +
    geom_hline(yintercept = baseline.female.02, color = '#33A02C', 
               linetype = 5, size = linewidth.base) +
    geom_hline(yintercept = baseline.female.03, color = '#1F78B4', 
               linetype = 5, size = linewidth.base) +
    theme_bw() + 
    theme(panel.background = element_rect(color = 'black',
                                          fill = 'transparent'),
          panel.grid = element_blank(),
          legend.key.height = unit(20, "pt"),
          legend.key.width = unit(30, "pt")) +
    ylab('Silhouette Distance') + xlab('Time (week)')
file.fit <- '/home/drizzle_zhang/microbiome/result/Figs/Beta.png'
ggsave(filename = file.fit, plot = plot.fit, units = 'cm', height = 10, width = 22)
