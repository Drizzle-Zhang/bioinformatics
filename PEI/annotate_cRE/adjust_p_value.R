library(metap)

split.pvalue <- function(str.vec) {
    num.vec <- as.numeric(str.vec)
    vec.pvalue <- 10 ^ -(num.vec)
    
    return(vec.pvalue)
}

cut.integration <- function(vec.pvalue, cutoff) {
    vec.pvalue <- vec.pvalue[vec.pvalue <= cutoff]
    if (length(vec.pvalue) == 0) {
        p.combine <- -1
    }
    if (length(vec.pvalue) == 1) {
        p.combine <- vec.pvalue
    }
    if (length(vec.pvalue) > 1) {
        p.combine <- sumlog(vec.pvalue)$p
    } 
    
    return(p.combine)
}

Adjust.pValue <- function(path.in, path.out, peak_num) {
    df.bed <- read.delim(path.in, sep = '\t', stringsAsFactors = F, header = F)
    list.split <- strsplit(df.bed[,'V5'], ',')
    list.pvalue <- lapply(list.split, split.pvalue)
    unlist.pvalue <- unlist(list.pvalue)
    vec.qvalue <- p.adjust(unlist.pvalue, 'BH', n = as.numeric(peak_num))
    filter.pvalue <- unlist.pvalue[vec.qvalue < 0.05]
    cutoff.pvalue <- max(filter.pvalue)
    vec.combine.p <- 
        unlist(lapply(list.pvalue, cut.integration, cutoff = cutoff.pvalue))
    df.bed$V6 <- vec.combine.p
    df.out <- df.bed[df.bed$V6 != -1, c('V1', 'V2', 'V3', 'V6')]
    df.out$V6 <- -log10(df.out$V6)
    df.out$V6[df.out$V6 == Inf] <- 500
    
    write.table(df.out, path.out, sep = '\t', quote = F, row.names = F,
                col.names = F)
}

args <- commandArgs(T)
Adjust.pValue(args[1], args[2], args[3])


