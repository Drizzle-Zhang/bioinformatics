scatter.plot <- function(file.in, fig.out) {
    df.lable.peak <- read.delim(file.in, sep = '\t')
    row.names(df.lable.peak) <- df.lable.peak$peak_id
    df.lable.peak$peak_id <- NULL
    
    d <- dist(t(df.lable.peak), method = 'binary')
    fit.average <- hclust(d, method = 'average')
    pdf(fig.out)
    plot(fit.average, hang = 0.2, cex = 1)
    dev.off()
    
}

args <- commandArgs(T)
scatter.plot(args[1], args[2])

