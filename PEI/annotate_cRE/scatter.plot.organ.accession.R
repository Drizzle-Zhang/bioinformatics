scatter.plot <- function(file.in, meta.in, path.out) {
    df.lable.peak <- read.delim(file.in, sep = '\t', header = F,
                                stringsAsFactors = F)
    df.meta <- read.delim(meta.in, sep = '\t', stringsAsFactors = F)
    label <- as.character(df.lable.peak[1, 2:dim(df.lable.peak)[2]])
    names(df.lable.peak) <- as.character(df.lable.peak[1,])
    df.lable.peak <- df.lable.peak[-1,]
    row.names(df.lable.peak) <- df.lable.peak$peak_id
    df.lable.peak$peak_id <- NULL
    df.lable.peak <- data.frame(lapply(df.lable.peak, as.numeric))
    
    df.scale <- t(scale(df.lable.peak))
    res.pca <- prcomp(df.scale)
    df.plot <- data.frame(res.pca$x[,1:2])
    pdf(paste0(path.out, '/scatter.organ.accession.pdf'))
    par(pin = c(6, 4))
    plot(res.pca$x[,1:2], main = 'Scatter Plot')
    text(res.pca$x[,1:2], label, cex = .4, pos = 1)
    dev.off()
    df.out <- cbind(res.pca$x[,1:2], label)
    write.table(df.out, paste0(path.out, '/scatter.organ.accession.txt'), 
                sep = '\t', quote = F, row.names = F, col.names = F)
    
}

args <- commandArgs(T)
scatter.plot(args[1], args[2], args[3])

