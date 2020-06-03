
# file.in = '/home/drizzle_zhang/driver_mutation/cRE_plot/model_test/score.txt'
# file.out =
#     '/home/drizzle_zhang/driver_mutation/cRE_plot/model_test/correct_score.txt'

correct.dhs.score <- function(file.in, file.out) {
    df.scores <- read.delim(file.in, sep = '\t', header = T,
                            row.names = 'peak_id', stringsAsFactors = F)
    df.scores[df.scores == '.'] <- NA
    df.score.num <- apply(df.scores, 1, as.numeric)
    if (is.null(dim(df.score.num))) {
        df.score.num <- as.data.frame(df.score.num)
        num.cols <- 1
    } else {
        df.score.num <- as.data.frame(t(df.score.num))
        num.cols <- dim(df.score.num)[2]
    }
    
    names(df.score.num) <- names(df.scores)
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

    write.table(df.out, file.out, sep = '\t', quote = F, row.names = T,
                col.names = F)
    
}

args <- commandArgs(T)
correct.dhs.score(args[1], args[2])
