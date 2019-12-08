library(ggplot2)
scatter.plot <- function(file.in, meta.in, path.out) {
    df.lable.peak <- read.delim(file.in, sep = '\t', header = T,
                                row.names = 'peak_id', stringsAsFactors = F)
    df.meta <- read.delim(meta.in, sep = '\t', header = T,
                          row.names = 'File.accession', stringsAsFactors = F)

    df.scale <- t(scale(df.lable.peak))
    res.pca <- prcomp(df.scale)
    PC1 <- res.pca$x[,1]
    PC2 <- res.pca$x[,2]
    accessions <- names(df.lable.peak)
    df.plot <- data.frame(PC1, PC2, accessions)
    df.plot <- merge(df.plot, df.meta, by = 'File.accession')
    obj.ggplot <- ggplot(
        data = df.plot, 
        aes(x = PC1, y = PC2, color = Biosample.term.id, 
            shape = Biosample.life.stage, text = File.accession)
        ) + 
        geom_point() + geom_text(aes(label = text))
    ggsave(filename = 'scatter.organ.accession.png', path = path.out, 
           width = 6, height = 4)
    write.table(df.plot, paste0(path.out, '/scatter.organ.accession.txt'), 
                sep = '\t', quote = F)
    
}

args <- commandArgs(T)
scatter.plot(args[1], args[2], args[3])

