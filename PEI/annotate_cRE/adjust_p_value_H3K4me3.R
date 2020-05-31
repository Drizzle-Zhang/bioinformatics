library(plyr)
library(VIM)
library(Hmisc)
library(car)

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
    
    if (num.cols <= 1) {
        df.score.num[,'peak_id'] <- row.names(df.score.num)
        
        df.score.num.na <- df.score.num[is.na(df.score.num[,1]),]
        df.na.0 <- as.data.frame(rep(0, dim(df.score.num.na)[1]))
        row.names(df.na.0) <- df.score.num.na[,'peak_id']
        names(df.na.0) <- c('score.combine')
        
        df.score.num.nona <- df.score.num[!is.na(df.score.num[,1]),]
        df.sort <- df.score.num.nona[order(df.score.num.nona[,1]),]
        func.quantile <- ecdf(df.sort[,1])
        df.sort.out <- as.data.frame(func.quantile(df.sort[,1]))
        row.names(df.sort.out) <- df.sort[,'peak_id']
        names(df.sort.out) <- c('score.combine')
        
        df.quantile <- rbind(df.sort.out, df.na.0)
        
    } else {
        names <- names(df.score.num)
        len.cols <- apply(df.score.num, 2, function(x){length(x[!is.na(x)])})
        
        # define length
        # len.median <- median(len.cols)
        names(len.cols) <- names
        len.cols <- sort(len.cols)
        # idx.ref <- ceiling(length(len.cols)/2)
        idx.ref <- length(len.cols)
        # define reference file
        name.ref <- names(len.cols[idx.ref])
        
        # normalization
        i = 1
        for (col in names) {
            df.sub <- df.score.num[,col]
            box <- summary(powerTransform(df.sub))
            index <- box$result[4]
            # if (index > 0 & col == name.ref) {
            #     print(col)
            #     print(index)
            #     return()
            # }
            # print(index)
            if (col == name.ref) {
                index.ref <- index
            }
            df.sub <- df.sub^index
            if (i == 1) {
                df.norm <- df.sub
            } else {
                df.norm <- cbind(df.norm, df.sub)
            }
            i = i + 1
        }
        # print(index.ref)
        dimnames(df.norm)[[1]] <- row.names(df.score.num)
        dimnames(df.norm)[[2]] <- names
        
        cols.correct <- setdiff(names, name.ref)
        df.out <- as.data.frame(df.norm[,name.ref])
        names(df.out) <- name.ref
        for (col in cols.correct) {
            df.sub <- df.norm[,c(name.ref, col)]
            df.sub.omitna <- as.data.frame(na.omit(df.sub))
            boxtidwell <- boxTidwell(as.formula(paste0(name.ref, ' ~ ', col)), 
                                     data = df.sub.omitna)
            df.sub.omitna$col.boxtidwell <- 
                (df.sub.omitna[,col])^boxtidwell$result[1]
            fit <- lm(as.formula(paste0(name.ref, ' ~ col.boxtidwell')), 
                      data = df.sub.omitna)
            df.sub.omitna$hatvalue <- hatvalues(fit)
            df.sub.omitna$rstudent <- rstudent(fit)
            cutoff.hatvalue <- 3*mean(df.sub.omitna$hatvalue)
            df.sub.omitna.del <- df.sub.omitna[
                (df.sub.omitna$hatvalue < cutoff.hatvalue) & 
                    (df.sub.omitna$rstudent < 2) & (df.sub.omitna$rstudent > -2),]
            fit <- lm(as.formula(paste0(name.ref, ' ~ col.boxtidwell')), 
                      data = df.sub.omitna.del)
            # print(summary(fit)$r.squared)
            
            sub.out <- 
                fit$coefficients[2]*(df.norm[,col]^boxtidwell$result[1]) + 
                fit$coefficients[1]
            sub.out <- as.data.frame(sub.out)
            names(sub.out) <- col
            df.out <- cbind(df.out, sub.out)
            
        }
        
        if (index.ref < 0) {
            df.impute <- df.out
            for (col in names) {
                df.impute[,col] <- impute(df.out[,col], fun = max)
            }
            df.mean <- rowMeans(df.impute, na.rm = T)
            # split NA and non-NA
            num.na <- max(df.mean)
            
            df.na <- df.mean[df.mean == num.na]
            df.na.0 <- as.data.frame(rep(0, length(df.na)))
            row.names(df.na.0) <- names(df.na)
            names(df.na.0) <- c('score.combine')
            
            df.nona <- df.mean[df.mean < num.na]
            df.sort <- sort(df.nona)
            df.sort <- as.data.frame(df.sort)
            func.quantile <- ecdf(df.sort[,1])
            df.sort.out <- 
                as.data.frame(1 - func.quantile(df.sort[,1]) + 1/dim(df.sort)[1])
            row.names(df.sort.out) <- row.names(df.sort)
            names(df.sort.out) <- c('score.combine')
            
            df.quantile <- rbind(df.sort.out, df.na.0)
            
        }
        if (index.ref > 0) {
            df.impute <- df.out
            for (col in names) {
                df.impute[,col] <- impute(df.out[,col], fun = min)
            }
            df.mean <- rowMeans(df.impute, na.rm = T)
            # split NA and non-NA
            num.na <- min(df.mean)
            
            df.na <- df.mean[df.mean == num.na]
            df.na.0 <- as.data.frame(rep(0, length(df.na)))
            row.names(df.na.0) <- names(df.na)
            names(df.na.0) <- c('score.combine')
            
            df.nona <- df.mean[df.mean > num.na]
            df.sort <- sort(df.nona, decreasing = T)
            df.sort <- as.data.frame(df.sort)
            func.quantile <- ecdf(df.sort[,1])
            df.sort.out <- as.data.frame(func.quantile(df.sort[,1]))
            row.names(df.sort.out) <- row.names(df.sort)
            names(df.sort.out) <- c('score.combine')
            
            df.quantile <- rbind(df.sort.out, df.na.0)
        }
    }
    
    return(df.quantile)
}


Adjust.pValue <- function(path.in, path.out, peak.num, file.num) {
    df.bed <- read.delim(path.in, sep = '\t', stringsAsFactors = F, header = F)
    df.pvalue <- data.frame()
    col.pvalue <- c()
    col.score <- c()
    for (i in 1:file.num) {
        col.pvalue <- c(col.pvalue, paste0('V', as.character(6 + 4*i)))
        col.score <- c(col.score, paste0('V', as.character(5 + 4*i)))
        col.use <- c(paste0('V', as.character(4 + 4*i)),
                     paste0('V', as.character(6 + 4*i)))
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
    df.score.correct <- correct.score(df.score)
    df.score.correct$V4 <- row.names(df.score.correct)
    df.bed <- merge(df.bed, df.score.correct, by = 'V4')
    df.bed[df.bed$p.combine == 0, 'score.combine'] <- 0

    df.out <- df.bed[, c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'score.combine', 
                         'p.combine')]
    write.table(df.out, path.out, sep = '\t', quote = F, row.names = F,
                col.names = F)
}

args <- commandArgs(T)
Adjust.pValue(args[1], args[2], args[3], args[4])


