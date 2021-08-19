library(biomaRt)

mart.snp <- useMart(host="http://grch37.ensembl.org",  # 这里用b37版本
                    biomart="ENSEMBL_MART_SNP",
                    dataset="hsapiens_snp")

getAnno <- function(rs = "rs3043732", mart = mart.snp) {
    results <- getBM(attributes =
                         c("chr_name","chrom_start","chrom_end", "refsnp_id","allele"),
                     filters = "snp_filter", values = rs, mart = mart)
    return(results)
}

file_disease <- '/mdshare/node9/zy/Reg_brain/DisGeNET/disease_sel.tsv'
df_disease <- read.delim(file_disease)

file_snp <- '/mdshare/node9/zy/Reg_brain/DisGeNET/curated_variant_disease_associations.tsv'
df_snp <- read.delim(file_snp)
df_snp_simple <- df_snp[, c("snpId", "diseaseId")]

path_disease_SNP <- '/mdshare/node9/zy/Reg_brain/disease_SNP/'
# function to generate snp files
generate_SNP_file <- function(df_snp_simple, path_disease_SNP, disease) {
    snps <- unique(df_snp_simple[df_snp_simple$diseaseId == disease, 'snpId'])
    if (length(snps) == 0) {return()}
    df.SNP.pos <- getAnno(snps)
    df.SNP.pos <- df.SNP.pos[, c("chr_name", "chrom_start", "chrom_end", "refsnp_id", "allele")]
    df.SNP.pos$numchr <- unlist(lapply(df.SNP.pos$chr_name, function(x) {nchar(x)}))
    df.SNP.pos <- df.SNP.pos[df.SNP.pos$numchr < 3,]
    df.SNP.pos$chr_name <- paste0('chr', df.SNP.pos$chr_name)
    # df.final.snp <- merge(df.SNP.pos, df.snp.simple, by.x = 'refsnp_id', by.y = 'SNPS')
    df.SNP.pos$chrom_end <- df.SNP.pos$chrom_end + 1
    df.SNP.pos <- df.SNP.pos[, c("chr_name", "chrom_start", "chrom_end", "refsnp_id", "allele")]
    df.SNP.pos$Disease <- rep(disease, nrow(df.SNP.pos))
    file.bed <- paste0(path_disease_SNP, disease, '.bed')
    write.table(df.SNP.pos, file.bed, sep = '\t',
                quote = F, row.names = F, col.names = F)
    return(file.bed)
}

nervous_diseases <- df_disease$diseaseId
for (disease in nervous_diseases) {
    generate_SNP_file(df_snp_simple, path_disease_SNP, disease)
}


# +- 1kb
path_disease_SNP <- '/mdshare/node9/zy/Reg_brain/disease_SNP_1kb/'
# function to generate snp files
generate_SNP_file <- function(df_snp_simple, path_disease_SNP, disease) {
    snps <- unique(df_snp_simple[df_snp_simple$diseaseId == disease, 'snpId'])
    if (length(snps) == 0) {return()}
    df.SNP.pos <- getAnno(snps)
    df.SNP.pos <- df.SNP.pos[, c("chr_name", "chrom_start", "chrom_end", "refsnp_id", "allele")]
    df.SNP.pos$numchr <- unlist(lapply(df.SNP.pos$chr_name, function(x) {nchar(x)}))
    df.SNP.pos <- df.SNP.pos[df.SNP.pos$numchr < 3,]
    df.SNP.pos$chr_name <- paste0('chr', df.SNP.pos$chr_name)
    df.SNP.pos$chrom_start <- df.SNP.pos$chrom_start - 1000
    df.SNP.pos$chrom_end <- df.SNP.pos$chrom_end + 1000
    df.SNP.pos <- df.SNP.pos[, c("chr_name", "chrom_start", "chrom_end", "refsnp_id", "allele")]
    df.SNP.pos$Disease <- rep(disease, nrow(df.SNP.pos))
    file.bed <- paste0(path_disease_SNP, disease, '.bed')
    write.table(df.SNP.pos, file.bed, sep = '\t',
                quote = F, row.names = F, col.names = F)
    return(file.bed)
}

nervous_diseases <- df_disease$diseaseId
for (disease in nervous_diseases) {
    generate_SNP_file(df_snp_simple, path_disease_SNP, disease)
}

