library(plyr)

# path.in <-
#     '/home/drizzle_zhang/driver_mutation/cRE_plot/model_test/DHS_promoter_H3K4me3.origin'

fisher.combine <- function(vec.lgp, cutoff.lgp) {
    vec.lgp[vec.lgp == '.'] <- '0'
    vec.lgp.num <- as.numeric(vec.lgp)
    vec.lgp.num[vec.lgp.num <= cutoff.lgp] <- 0
    vec.lnp <- vec.lgp.num * log(10)
    chisq <- sum(2 * vec.lnp)
    df <- 2 * length(vec.lgp)
    lnp.combine <- pchisq(chisq, df, lower.tail = F, log.p = T)
    lgp.combine <- -lnp.combine * log10(exp(1))
    if (lgp.combine <= cutoff.lgp) {lgp.combine <- 0}

    return(lgp.combine)
}


correct.score <- function(df.score) {
    df.score[df.score == '.'] <- NA
    df.score.num <- apply(df.score, 1, as.numeric)
    if (is.null(dim(df.score.num))) {
        df.score.num <- as.data.frame(df.score.num)
        num.cols <- 1
    } else {
        df.score.num <- as.data.frame(t(df.score.num))
        num.cols <- dim(df.score.num)[2]
    }
    
    names <- names(df.score.num)
    i = 1
    for (col in names) {
        df.sub <- as.data.frame(df.score.num[,col])
        df.sub[,'peak_id'] <- row.names(df.score.num)
        
        df.sub.na <- df.sub[is.na(df.sub[,1]),]
        df.na.0 <- as.data.frame(rep(0, dim(df.sub.na)[1]))
        df.na.0$peak_id <- df.sub.na[,'peak_id']
        names(df.na.0) <- c(col,'peak_id')
        
        df.sub.nona <- df.sub[!is.na(df.sub[,1]),]
        df.sort <- df.sub.nona[order(df.sub.nona[,1]),]
        func.quantile <- ecdf(df.sort[,1])
        df.sort.out <- as.data.frame(func.quantile(df.sort[,1]))
        df.sort.out$peak_id <- df.sort[,'peak_id']
        names(df.sort.out) <- c(col,'peak_id')
        
        df.quantile <- rbind(df.sort.out, df.na.0)
        
        if (i == 1) {
            df.scores.quantile <- df.quantile
        } else {
            df.scores.quantile <- 
                merge.data.frame(df.scores.quantile, df.quantile, 
                                 by = 'peak_id')
        }
        i = i + 1
    }
    
    row.names(df.scores.quantile) <- df.scores.quantile$peak_id
    df.scores.quantile$peak_id <- NULL
    
    df.max.quantile <- apply(df.scores.quantile, 1, max)
    df.out <- as.data.frame(df.max.quantile)
    
    return(df.max.quantile)
    
}


Adjust.pValue <- function(path.in, path.out, peak.num, file.num, ref.col.num) {
    df.bed <- read.delim(path.in, sep = '\t', stringsAsFactors = F, header = F)
    df.pvalue <- data.frame()
    col.pvalue <- c()
    col.score <- c()
    ref.col <- as.numeric(ref.col.num) - 4
    for (i in 1:as.numeric(file.num)) {
        col.pvalue <- c(col.pvalue, 
                        paste0('V', as.character(ref.col + 3 + 4*i)))
        col.score <- c(col.score, 
                       paste0('V', as.character(ref.col + 2 + 4*i)))
        col.use <- c(paste0('V', as.character(ref.col + 1 + 4*i)),
                     paste0('V', as.character(ref.col + 3 + 4*i)))
        sub.df.pvalue <- df.bed[, col.use]
        names(sub.df.pvalue) <- c('Label', 'P.value')
        df.pvalue <- rbind(df.pvalue, sub.df.pvalue)
    }
    # calculate cutoff
    df.pvalue <- unique(df.pvalue[df.pvalue$Label != '.',])
    # print(df.pvalue[1:5,])
    list.lgp <- as.numeric(df.pvalue$P.value)
    vec.pvalue <- 10 ^ -(list.lgp)
    vec.qvalue <- p.adjust(vec.pvalue, 'BH', n = as.numeric(peak.num))
    filter.pvalue <- vec.pvalue[vec.qvalue < 0.05]
    cutoff.lgp <- -log10(max(filter.pvalue))
    # combine p value
    df.col.pvalue <- df.bed[, col.pvalue]
    # print(cutoff.lgp)
    vec.combine.lgp <- 
        unlist(alply(.data = df.col.pvalue, .margins = 1, .fun = fisher.combine, 
                      cutoff.lgp = cutoff.lgp))
    df.bed$p.combine <- vec.combine.lgp
    # combine score
    df.score <- as.data.frame(df.bed[, col.score])
    row.names(df.score) <- df.bed$V4
    df.score.correct <- as.data.frame(correct.score(df.score))
    names(df.score.correct) <- c('score.combine')
    df.score.correct$V4 <- row.names(df.score.correct)
    df.bed <- merge(df.bed, df.score.correct, by = 'V4')
    df.bed[df.bed$p.combine == 0, 'score.combine'] <- 0

    df.out <- df.bed[, c(paste0('V', as.character(1:ref.col.num)), 
                         'score.combine', 'p.combine')]
    write.table(df.out, path.out, sep = '\t', quote = F, row.names = F,
                col.names = F)
}

args <- commandArgs(T)
Adjust.pValue(args[1], args[2], args[3], args[4], args[5])


