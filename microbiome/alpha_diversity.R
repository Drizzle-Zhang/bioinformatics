# alpha diversity analysis
library(amplicon)

file.in <- '/home/drizzle_zhang/microbiome/result/4.Alpha_Diversity/alpha_index_table/alpha_estimator_summary.csv'
mat.alpha <- read.table(file.in, sep = ',', header = T, row.names = 1)

# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.out.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)
df.meta$Dose <- as.factor(df.meta$Dose)
row.names(df.meta) <- df.meta$Sample

# normalization
mat.combine <- cbind(mat.alpha, df.meta)
type.dose <- unique(mat.combine$Dose)
mat.alpha.norm <- data.frame()
for (dose in type.dose) {
    baseline.alpha <- mat.combine[
        ((mat.combine$Dose == dose) & (mat.combine$Time == 'A')), 
        c("observed_species", "shannon", "simpson")]
    baseline.alpha <- apply(baseline.alpha, 2, median)
    sub.mat.alpha <- mat.combine[
        mat.combine$Dose == dose, c("observed_species", "shannon", "simpson")]
    sub.alpha.norm <- t(apply(sub.mat.alpha, 1, function(x) {x / baseline.alpha}))
    mat.alpha.norm <- rbind(mat.alpha.norm, sub.alpha.norm)
}


# time series
path.plot <- '/home/drizzle_zhang/microbiome/result/4.Alpha_Diversity/alpha_boxplot'
series.time <- unique(df.meta$Time)
df.plot.fit <- data.frame()
i = 1
for (sub.time in series.time) {
    # select meta
    sel.meta <- df.meta
    sel.meta <- df.meta[df.meta$Time == sub.time,]
    # sel.meta <- sel.meta[sel.meta$Dose %in% c(0, 3),]
    row.names(sel.meta) <- sel.meta$Sample
    
    # select sample
    use.sample <- sel.meta$Sample
    mat.plot.in <- data.frame()
    for (alpha_index in c("observed_species", "shannon", "simpson")) {
        sub.mat <- data.frame(
            value = mat.alpha.norm[use.sample, alpha_index],
            type.alpha = rep(alpha_index, length(use.sample)),
            row.names = use.sample)
        sub.mat <- cbind(sub.mat, sel.meta)
        mat.plot.in <- rbind(mat.plot.in, sub.mat)
    }

    # boxplot
    plot.alpha <- 
        ggplot(aes(x = Dose, y = value, color = Dose, shape = Dose), 
               data = mat.plot.in) + 
        geom_boxplot() + 
        facet_wrap(. ~ type.alpha, scales = 'free') + 
        labs(x = '', y = 'Alpha Diversity Measure') + 
        theme(panel.background = element_rect(color = 'gray', 
                                              fill = 'transparent'))
    
    ggsave(plot = plot.alpha, path = path.plot, 
           filename = paste0(sub.time, '_0123_norm.png'))

    # diff of alpha
    mat.shannon <- mat.plot.in[
        mat.plot.in$type.alpha == 'shannon', c('value', 'Dose')]
    for (dose in type.dose) {
        sub.shannon <- median(
            mat.shannon[mat.shannon$Dose == dose, 'value'])
        df.plot.fit <-
            rbind(df.plot.fit,
                  data.frame(Shannon = sub.shannon,
                             Dose = dose, Time = i))
    }

    i = i + 1
    
}

plot.fit <- 
    ggplot(data = df.plot.fit, aes(x = Time, y = Shannon, color = Dose)) + 
    geom_line() + 
    geom_point()

