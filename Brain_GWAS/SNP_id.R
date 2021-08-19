library(biomaRt)
library(openxlsx)

mart.snp <- useMart(host="http://grch37.ensembl.org",  # 这里用b37版本
                    biomart="ENSEMBL_MART_SNP",
                    dataset="hsapiens_snp")

getAnno <- function(rs = "rs3043732", mart = mart.snp) {
    results <- getBM(attributes =
                         c("chr_name","chrom_start","chrom_end", "refsnp_id","allele"),
                     filters = "snp_filter", values = rs, mart = mart)
    return(results)
}

# get brain SNPs
df.snp <- read.xlsx('/mdshare/node9/zy/Brain_GWAS/brain_disease.xlsx')
df.snp.simple <- df.snp[, c("DISEASE/TRAIT", "SNPS", "CONTEXT", "CHR_ID", "CHR_POS")]
df.snp.simple <- df.snp.simple[!duplicated(df.snp.simple),]
df.snp.simple <- na.omit(df.snp.simple)
df.snp.simple$numchr <- unlist(lapply(df.snp.simple$CHR_ID, function(x) {nchar(x)}))
df.snp.simple <- df.snp.simple[df.snp.simple$numchr < 3,]
snps <- unique(df.snp.simple$SNPS)

# get SNP positions
df.SNP.pos <- getAnno(snps)
df.SNP.pos <- df.SNP.pos[, c("chr_name", "chrom_start", "chrom_end", "refsnp_id", "allele")]
df.SNP.pos$numchr <- unlist(lapply(df.SNP.pos$chr_name, function(x) {nchar(x)}))
df.SNP.pos <- df.SNP.pos[df.SNP.pos$numchr < 3,]
df.SNP.pos$chr_name <- paste0('chr', df.SNP.pos$chr_name)
# df.final.snp <- merge(df.SNP.pos, df.snp.simple, by.x = 'refsnp_id', by.y = 'SNPS')
df.SNP.pos$chrom_end <- df.SNP.pos$chrom_end + 1
df.SNP.pos <- df.SNP.pos[, c("chr_name", "chrom_start", "chrom_end", "refsnp_id", "allele")]
file.bed <- '/mdshare/node9/zy/Brain_GWAS/Brain_SNP_hg19.bed'
write.table(df.SNP.pos, file.bed, sep = '\t',
            quote = F, row.names = F, col.names = F)

# +- 1kb
df.SNP.pos <- getAnno(snps)
df.SNP.pos <- df.SNP.pos[, c("chr_name", "chrom_start", "chrom_end", "refsnp_id", "allele")]
df.SNP.pos$numchr <- unlist(lapply(df.SNP.pos$chr_name, function(x) {nchar(x)}))
df.SNP.pos <- df.SNP.pos[df.SNP.pos$numchr < 3,]
df.SNP.pos$chr_name <- paste0('chr', df.SNP.pos$chr_name)
df.SNP.pos$chrom_start <- df.SNP.pos$chrom_start - 1000
df.SNP.pos$chrom_end <- df.SNP.pos$chrom_end + 1000
df.SNP.pos <- df.SNP.pos[, c("chr_name", "chrom_start", "chrom_end", "refsnp_id", "allele")]
file.bed <- '/mdshare/node9/zy/Brain_GWAS/Brain_SNP_hg19_1kb.bed'
write.table(df.SNP.pos, file.bed, sep = '\t',
            quote = F, row.names = F, col.names = F)


### hg38
# mart.snp <- useMart(host="http://asia.ensembl.org",  # 这里用b37版本
#                     biomart="ENSEMBL_MART_SNP",
#                     dataset="hsapiens_snp")

# get brain SNPs
df.snp <- read.xlsx('/mdshare/node9/zy/Brain_GWAS/brain_disease.xlsx')
df.snp.simple <- df.snp[, c("DISEASE/TRAIT", "SNPS", "CONTEXT", "CHR_ID", "CHR_POS")]
df.snp.simple <- df.snp.simple[!duplicated(df.snp.simple),]
df.snp.simple <- na.omit(df.snp.simple)
df.snp.simple$numchr <- unlist(lapply(df.snp.simple$CHR_ID, function(x) {nchar(x)}))
df.snp.simple <- df.snp.simple[df.snp.simple$numchr < 3,]
df.SNP.pos <- data.frame(chr_name = df.snp.simple$CHR_ID,
                         chrom_start = df.snp.simple$CHR_POS,
                         chrom_end = as.numeric(df.snp.simple$CHR_POS) + 1,
                         refsnp_id = df.snp.simple$SNPS)
df.SNP.pos$chr_name <- paste0('chr', df.SNP.pos$chr_name)
file.bed <- '/mdshare/node9/zy/Brain_GWAS/Brain_SNP_hg38.bed'
write.table(df.SNP.pos, file.bed, sep = '\t',
            quote = F, row.names = F, col.names = F)

# +- 1kb
df.SNP.pos <- data.frame(chr_name = df.snp.simple$CHR_ID,
                         chrom_start = as.numeric(df.snp.simple$CHR_POS) - 1000,
                         chrom_end = as.numeric(df.snp.simple$CHR_POS) + 1000,
                         refsnp_id = df.snp.simple$SNPS)
df.SNP.pos$chr_name <- paste0('chr', df.SNP.pos$chr_name)
file.bed <- '/mdshare/node9/zy/Brain_GWAS/Brain_SNP_hg38_1kb.bed'
write.table(df.SNP.pos, file.bed, sep = '\t',
            quote = F, row.names = F, col.names = F)
