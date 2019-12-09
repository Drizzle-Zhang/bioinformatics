library(ggplot2)
scatter.plot <- function(file.in, str.head, meta.in, path.out) {
    labels <- as.character(unlist(strsplit(str.head, '_', fixed = T)))
    list.label <- strsplit(labels, '|', fixed = T)
    df.label <- data.frame(
        matrix(unlist(list.label), ncol = 3, byrow = T), 
        stringsAsFactors = FALSE
        )
    names(df.label) <- c("Biosample.organ.1", "Biosample.life.stage",
                         "Biosample.term.name")
    df.label$labels <- labels
    
    df.lable.peak <- read.delim(file.in, sep = '\t', header = T,
                                row.names = 'peak_id', stringsAsFactors = F)
    df.meta <- read.delim(meta.in, sep = '\t', header = T, stringsAsFactors = F)
    
    meta.labels <- merge(
        df.meta, df.label, 
        by = c("Biosample.life.stage", "Biosample.term.name")
        )
    
    df.scale <- t(scale(df.lable.peak))
    res.pca <- prcomp(df.scale)
    PC1 <- res.pca$x[,1]
    PC2 <- res.pca$x[,2]
    df.plot <- data.frame(PC1, PC2, labels)
    df.plot <- merge(df.plot, meta.labels, by = 'labels')
    obj.ggplot <- ggplot(
        data = df.plot, 
        aes(x = PC1, y = PC2, color = Biosample.organ, 
            shape = Biosample.life.stage)
    ) + 
        geom_point(size = 3, alpha = 0.3) + 
        geom_text(aes(label = labels), size = 2)
    ggsave(plot = obj.ggplot, filename = 'PCA.all.pdf', 
           path = path.out, width = 20, height = 15)
    write.table(df.plot, paste0(path.out, '/PCA.all.txt'), 
                sep = '\t', quote = F)
    
}

args <- commandArgs(T)
scatter.plot(args[1], args[2], args[3], args[4])

