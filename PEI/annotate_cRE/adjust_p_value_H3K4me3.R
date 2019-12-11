fisher.combine <- function(vec.lgp, cutoff.lgp) {
    vec.lgp[vec.lgp == '.'] <- '0'
    vec.lgp.num <- as.numeric(vec.lgp)
    vec.lgp.num[vec.lgp <= cutoff.lgp] <- 0
    vec.lnp <- vec.lgp.num * log(10)
    chisq <- sum(2 * vec.lnp)
    df <- 2 * length(vec.lgp)
    lnp.combine <- pchisq(chisq, df, lower.tail = F, log.p = T)
    lgp.combine <- -lnp.combine * log10(exp(1))
    if (lgp.combine <= cutoff.lgp) {lgp.combine <- 0}

    return(lgp.combine)
}

Adjust.pValue <- function(path.in, path.out, peak.num, file.num) {
    df.bed <- read.delim(path.in, sep = '\t', stringsAsFactors = F, header = F)
    df.pvalue <- data.frame()
    col.pvalue <- c()
    for (i in 1:file.num) {
        col.pvalue <- c(col.pvalue, paste0('V', as.character(5 + 3*i)))
        col.use <- c(paste0('V', as.character(4 + 3*i)),
                     paste0('V', as.character(5 + 3*i)))
        sub.df.pvalue <- df.bed[, col.use]
        names(sub.df.pvalue) <- c('Label', 'P.value')
        df.pvalue <- rbind(df.pvalue, sub.df.pvalue)
    }
    # calculate cutoff
    df.pvalue <- unique(df.pvalue[df.pvalue$Label != '.',])
    list.lgp <- as.numeric(df.pvalue$P.value)
    vec.pvalue <- 10 ^ -(list.lgp)
    vec.qvalue <- p.adjust(vec.pvalue, 'BH', n = as.numeric(peak.num))
    filter.pvalue <- vec.pvalue[vec.qvalue < 0.05]
    cutoff.lgp <- -log10(max(filter.pvalue))
    # combine p value
    df.pvalue <- df.bed[, col.pvalue]
    vec.combine.lgp <- 
        unlist(alply(.data = df.pvalue, .margins = 1, .fun = fisher.combine, 
                      cutoff.lgp = cutoff.lgp))
    df.bed$p.combine <- vec.combine.lgp
    df.out <- df.bed[, c('V1', 'V2', 'V3', 'V4', 'V5', 'p.combine')]

    write.table(df.out, path.out, sep = '\t', quote = F, row.names = F,
                col.names = F)
}

args <- commandArgs(T)
Adjust.pValue(args[1], args[2], args[3], args[4])


