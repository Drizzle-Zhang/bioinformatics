library(VIM)
library(Hmisc)
library(car)

# file.in = '/home/drizzle_zhang/driver_mutation/cRE_plot/model_test/score.txt'
# file.out = 
#     '/home/drizzle_zhang/driver_mutation/cRE_plot/model_test/correct_score.txt'
correct.dhs.score <- function(file.in, file.out) {
    df.scores <- read.delim(file.in, sep = '\t', header = T,
                            row.names = 'peak_id', stringsAsFactors = F)
    
    names <- names(df.scores)
    len.cols <- apply(df.scores, 2, function(x){length(x[!is.na(x)])})
    
    # define length
    len.median <- median(len.cols)
    names(len.cols) <- names
    len.cols <- sort(len.cols)
    idx.ref <- ceiling(length(len.cols)/2)
    # define reference file
    name.ref <- names(len.cols[idx.ref])
    
    # normalization
    i = 1
    for (col in names) {
        df.sub <- df.scores[,col]
        box <- summary(powerTransform(df.sub))
        index <- box$result[4]
        if (index > 0) {
            print(col)
            print(index)
            return()
        }
        # print(index)
        df.sub <- df.sub^index
        if (i == 1) {
            df.norm <- df.sub
        } else {
            df.norm <- cbind(df.norm, df.sub)
        }
        i = i + 1
    }
    dimnames(df.norm)[[1]] <- row.names(df.scores)
    dimnames(df.norm)[[2]] <- names
    
    cols.correct <- setdiff(names, name.ref)
    df.out <- as.data.frame(df.norm[,name.ref])
    names(df.out) <- name.ref
    for (col in cols.correct) {
        df.sub <- df.norm[,c(name.ref, col)]
        df.sub.omitna <- as.data.frame(na.omit(df.sub))
        boxtidwell <- boxTidwell(as.formula(paste0(name.ref, ' ~ ', col)), 
                                 data = df.sub.omitna)
        df.sub.omitna$col.boxtidwell <- (df.sub.omitna[,col])^boxtidwell$result[1]
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
        
        sub.out <- fit$coefficients[2]*(df.norm[,col]^boxtidwell$result[1]) + 
            fit$coefficients[1]
        sub.out <- as.data.frame(sub.out)
        names(sub.out) <- col
        df.out <- cbind(df.out, sub.out)
        
    }
    
    df.impute <- df.out
    for (col in names) {
        df.impute[,col] <- impute(df.out[,col], fun = max)
    }
    df.mean <- rowMeans(df.impute, na.rm = T)
    df.sort <- sort(df.mean)
    df.sort <- as.data.frame(df.sort[1:len.median])
    func.quantile <- ecdf(df.sort[1:len.median])
    df.sort.out <- as.data.frame(1 - func.quantile(df.sort) + 1/len.median)
    row.names(df.sort.out) <- names(df.sort)
    
    write.table(df.sort.out, file.out, sep = '\t', quote = F, row.names = T,
                col.names = F)
    
}

args <- commandArgs(T)
correct.dhs.score(args[1], args[2])
