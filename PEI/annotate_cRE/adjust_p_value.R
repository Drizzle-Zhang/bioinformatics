split.pvalue <- function(str.vec) {
    num.vec <- as.numeric(str.vec)
    vec.pvalue <- 10 ^ -(num.vec)
    
    return(vec.pvalue)
}

cut.integration <- function(vec.lgp, cutoff) {
    vec.lgp <- vec.lgp[vec.lgp >= cutoff]
    len.vec <- length(vec.lgp)
    if (len.vec == 0) {
        lgp.combine <- -1
    }
    if (len.vec == 1) {
        lgp.combine <- vec.lgp
    }
    if (len.vec > 1) {
        vec.lnp <- vec.lgp * log(10)
        chisq <- sum(2 * vec.lnp)
        df <- 2 * len.vec
        lnp.combine <- pchisq(chisq, df, lower.tail = F, log.p = T)
        lgp.combine <- -lnp.combine * log10(exp(1))
    } 
    
    return(lgp.combine)
}

Adjust.pValue <- function(path.in, path.out, peak_num) {
    df.bed <- read.delim(path.in, sep = '\t', stringsAsFactors = F, header = F)
    list.split <- strsplit(df.bed[,'V5'], ',')
    # calculate cutoff
    list.pvalue <- lapply(list.split, split.pvalue)
    unlist.pvalue <- unlist(list.pvalue)
    vec.qvalue <- p.adjust(unlist.pvalue, 'BH', n = as.numeric(peak_num))
    filter.pvalue <- unlist.pvalue[vec.qvalue < 0.05]
    cutoff.lgp <- -log10(max(filter.pvalue))
    # combine p value
    list.lgp <- lapply(list.split, as.numeric)
    vec.combine.lgp <- 
        unlist(lapply(list.lgp, cut.integration, cutoff = cutoff.lgp))
    df.bed$V6 <- vec.combine.lgp
    df.out <- df.bed[df.bed$V6 != -1, c('V1', 'V2', 'V3', 'V6')]

    write.table(df.out, path.out, sep = '\t', quote = F, row.names = F,
                col.names = F)
}

args <- commandArgs(T)
Adjust.pValue(args[1], args[2], args[3])


