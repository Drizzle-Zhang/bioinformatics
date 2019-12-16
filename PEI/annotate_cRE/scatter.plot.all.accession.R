library(ggplot2)
scatter.plot <- function(file.in, meta.in, path.out) {
    df.lable.peak <- read.delim(file.in, sep = '\t', header = T,
                                row.names = 'peak_id', stringsAsFactors = F)
    df.meta <- read.delim(meta.in, sep = '\t', header = T, stringsAsFactors = F)
    
    df.scale <- t(scale(df.lable.peak))
    res.pca <- prcomp(df.scale)
    PC1 <- res.pca$x[,1]
    PC2 <- res.pca$x[,2]
    File.accession <- names(df.lable.peak)
    df.plot <- data.frame(PC1, PC2, File.accession)
    df.plot <- merge(df.plot, df.meta, by = 'File.accession')
    obj.ggplot <- ggplot(
        data = df.plot, 
        aes(x = PC1, y = PC2, color = Biosample.organ, 
            shape = Biosample.life.stage)
        ) + 
        geom_point(size = 3) + 
        geom_text(aes(label = File.accession), size = 2)
    ggsave(plot = obj.ggplot, filename = 'PCA.all.accession.pdf', 
           path = path.out, width = 20, height = 15)
    write.table(df.plot, paste0(path.out, '/PCA.all.accession.txt'), 
                sep = '\t', quote = F)
    
    
    # filter
    df.meta.filter <- df.meta[
        df.meta$Biosample.life.stage %in% c('child', 'adult', 'newborn'),]
    df.scale <- t(scale(df.lable.peak[, df.meta.filter$File.accession]))
    res.pca <- prcomp(df.scale)
    PC1 <- res.pca$x[,1]
    PC2 <- res.pca$x[,2]
    File.accession <- row.names(df.scale)
    df.plot <- data.frame(PC1, PC2, File.accession)
    df.plot <- merge(df.plot, df.meta.filter, by = 'File.accession')
    obj.ggplot <- ggplot(
        data = df.plot, 
        aes(x = PC1, y = PC2, color = Biosample.organ, 
            shape = Biosample.life.stage)
    ) + 
        geom_point(size = 3) + 
        geom_text(aes(label = File.accession), size = 2)
    ggsave(plot = obj.ggplot, filename = 'PCA.all.adult.accession.pdf', 
           path = path.out, width = 20, height = 15)
    write.table(df.plot, paste0(path.out, '/PCA.all.adult.accession.txt'), 
                sep = '\t', quote = F)
    
    
}

args <- commandArgs(T)
scatter.plot(args[1], args[2], args[3])

