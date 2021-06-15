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
df.snp <- read.xlsx('/mdshare/node9/zy/GWAS/brain_disease.xlsx')
df.snp.simple <- df.snp[, c("DISEASE/TRAIT", "SNPS", "CONTEXT")]
df.snp.simple <- df.snp.simple[!duplicated(df.snp.simple),]
snps <- unique(df.snp.simple$SNPS)

# get SNP positions
df.SNP.pos <- getAnno(snps)
df.SNP.pos <- df.SNP.pos[, c("chr_name","chrom_start","chrom_end", "refsnp_id","allele")]
file.bed <- '/mdshare/node9/zy/GWAS/Brain_SNP.bed'
write.table(df.SNP.pos, file.bed, sep = '\t',
            quote = F, row.names = F, col.names = F)
