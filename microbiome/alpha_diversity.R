# alpha diversity analysis
library(amplicon)

file.in <- '/home/drizzle_zhang/microbiome/result/4.Alpha_Diversity/alpha_index_table/alpha_estimator_summary.csv'
mat.alpha <- read.table(file.in, sep = ',', header = T, row.names = 1)

# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.out.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)

# time series
path.plot <- '/home/drizzle_zhang/microbiome/result/4.Alpha_Diversity/alpha_boxplot'
series.time <- unique(df.meta$Time)
for (sub.time in series.time) {
    # select meta
    sel.meta <- df.meta
    sel.meta <- df.meta[df.meta$Time == sub.time,]
    sel.meta <- sel.meta[sel.meta$Dose %in% c(0, 3),]
    row.names(sel.meta) <- sel.meta$Sample
    
    # select sample
    use.sample <- sel.meta$Sample
    mat.plot.in <- data.frame()
    for (alpha_index in c("observed_species", "shannon", "simpson")) {
        sub.mat <- data.frame(
            value = mat.alpha[use.sample, alpha_index],
            type.alpha = rep(alpha_index, length(use.sample)),
            row.names = use.sample)
        sub.mat <- cbind(sub.mat, sel.meta)
        mat.plot.in <- rbind(mat.plot.in, sub.mat)
    }

    # boxplot
    plot.alpha <- 
        ggplot(aes(x = Group, y = value, color = Group, shape = Group), 
               data = mat.plot.in) + 
        geom_boxplot() + 
        facet_wrap(. ~ type.alpha, scales = 'free') + 
        labs(x = '', y = 'Alpha Diversity Measure') + 
        theme(panel.grid.minor.y = element_line(colour = "gray"),
              panel.grid.minor.x = element_line(colour = "gray"),
              panel.background = element_rect(color = 'gray', 
                                              fill = 'transparent'))
    
    ggsave(plot = plot.alpha, path = path.plot, 
           filename = paste0(sub.time, '_0vs3.png'))
    
}



