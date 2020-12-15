# former
file_ko_1 <- '~/microbiome/result_1/9.PICRUSt/ko_predict/ko_predictions.txt'
df_ko_1 <- read.delim(file_ko_1, row.names = 1)
df_ko_1 <- df_ko_1[, 1:505]

# later
file_ko_2 <- '~/microbiome/result/9.PICRUSt/ko_predict/ko_predictions.txt'
df_ko_2 <- read.delim(file_ko_2, row.names = 1)
df_ko_2 <- df_ko_2[, 1:743]

colnames1 <- colnames(df_ko_1)
colnames2 <- colnames(df_ko_2)
colnames2.sel <- setdiff(colnames2, colnames1)

df_ko_combine <- cbind(df_ko_1, df_ko_2[, colnames2.sel])
file_combine <- '~/microbiome/result/9.PICRUSt/ko_predict/ko_predictions_combine.txt'
write.table(df_ko_combine, file_combine, sep = '\t', quote = F, 
            row.names = T, col.names = T)
