library(ggplot2)
scatter.plot <- function(file.in, path.out) {
    df.lable.peak <- read.delim(file.in, sep = '\t', header = F,
                                stringsAsFactors = F)
    list.label <- strsplit(
        as.character(df.lable.peak[1, 2:dim(df.lable.peak)[2]]), '|', fixed = T)
    label <- as.character(unlist(
        lapply(list.label, function(x){paste(x[2:3], collapse = ' ')})
    ))
    
    df.lable.peak <- read.delim(file.in, sep = '\t', header = T,
                                row.names = 'peak_id', stringsAsFactors = F)
    
    df.scale <- t(scale(df.lable.peak))
    res.pca <- prcomp(df.scale)
    PC1 <- res.pca$x[,1]
    PC2 <- res.pca$x[,2]
    len.PC1 <- max(PC1) - min(PC1)
    down.x <- min(PC1) - 0.1*len.PC1
    up.x <- max(PC1) + 0.1*len.PC1
    len.PC2 <- max(PC2) - min(PC2)
    down.y <- min(PC2) - 0.1*len.PC2
    up.y <- max(PC2) + 0.1*len.PC2
    df.plot <- data.frame(PC1, PC2, label)
    obj.ggplot <- ggplot(
        data = df.plot, 
        aes(x = PC1, y = PC2)
    ) + 
        geom_point(size = 2, alpha = 0.5, shape = 1) + 
        geom_text(aes(label = label), size = 2.5) + 
        xlim(down.x, up.x) + ylim(down.y, up.y)
    ggsave(plot = obj.ggplot, filename = 'scatter.organ.pdf', 
           path = path.out, width = 10, height = 7)
    write.table(df.plot, paste0(path.out, '/scatter.organ.txt'), 
                sep = '\t', quote = F)

}

args <- commandArgs(T)
scatter.plot(args[1], args[2])

