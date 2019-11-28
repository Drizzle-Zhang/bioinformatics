scatter.plot <- function(file.in, fig.out) {
    df.lable.peak <- read.delim(file.in, sep = '\t', header = F,
                                stringsAsFactors = F)
    label <- as.character(df.lable.peak[1, 2:dim(df.lable.peak)[2]])
    list.label <- strsplit(label, '|', fixed = T)
    label <- unlist(
        lapply(list.label, function(x){paste(x[2:3], collapse = '_')})
        )
    names(df.lable.peak) <- as.character(df.lable.peak[1,])
    df.lable.peak <- df.lable.peak[-1,]
    row.names(df.lable.peak) <- df.lable.peak$peak_id
    df.lable.peak$peak_id <- NULL
    df.lable.peak <- data.frame(lapply(df.lable.peak, as.numeric))
    
    df.binary <- df.lable.peak
    df.binary[df.binary != 0] = 1
    df.scale <- t(scale(df.binary))
    res.pca <- prcomp(df.scale)
    pdf(fig.out)
    par(pin = c(6, 4))
    plot(res.pca$x[,1:2], main = 'Scatter Plot')
    text(res.pca$x[,1:2], label, cex = .6, pos = 1)
    dev.off()
    df.out <- cbind(res.pca$x[,1:2], label)
    write.table(df.out, paste0(fig.out, '.txt'), sep = '\t', quote = F, 
                row.names = F, col.names = F)
    
}

args <- commandArgs(T)
scatter.plot(args[1], args[2])

