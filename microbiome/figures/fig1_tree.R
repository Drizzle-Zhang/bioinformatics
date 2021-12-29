library(RColorBrewer)
library(foreach)
library(doParallel)
library(stringr)
registerDoParallel(4)

# parameters
# filter by abundance
abundance <- 0.0001
# top num diff OTU
topN <- 50
# path
pathout <- '/home/drizzle_zhang/microbiome/result/6.Phylogenetic/graphlan'

# read OTUs matrix
file.OTUs <- '/home/drizzle_zhang/microbiome/result/2.OTUs/OTUs_stat/OTUs_tax_even.txt'
df.OTUs <- read.delim(file.OTUs, row.names = 1)
df.OTUs <- df.OTUs[, 1:dim(df.OTUs)[2] - 1]

# read taxonomy file
file.tax <- '/home/drizzle_zhang/microbiome/result/2.OTUs/OTUs_stat/taxonomy.txt'
df.taxonomy <- read.delim(file.tax, row.names = 1, stringsAsFactors = F)

# filter
# normalization
df.OTUs.norm <- as.data.frame(t(t(df.OTUs)/colSums(df.OTUs,na.rm = T)*100))
idx = order(rowMeans(df.OTUs.norm), decreasing = T)
df.OTUs.norm = df.OTUs.norm[idx,]
# filter by abundance
idx = rowMeans(df.OTUs.norm) > abundance
filtered_otutab = df.OTUs.norm[idx,]
# filter by top num
# filtered_otutab = head(df.OTUs.norm, topN)
# add means
filtered_otutab = round(cbind(rowMeans(filtered_otutab), filtered_otutab), digits = 4)
colnames(filtered_otutab)[1] = "Mean"
# filter taxonomy
idx = rownames(filtered_otutab) %in% rownames(df.taxonomy)
filtered_otutab = filtered_otutab[idx,]
filtered_taxonomy = df.taxonomy[rownames(filtered_otutab),]

# save filter results
setwd(pathout)
saveRDS(filtered_otutab, file = paste0('OTUs_filtered_', abundance, '.Rdata'))
saveRDS(filtered_taxonomy, file = paste0('taxonomy_filtered_', abundance, '.Rdata'))

# read Rdata
setwd(pathout)
df.OTUs <- readRDS(paste0('OTUs_filtered_', abundance, '.Rdata'))
df.taxonomy <- readRDS(paste0('taxonomy_filtered_', abundance, '.Rdata'))

# difference analysis
# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)

# delete unclassified families
idx.del <- df.taxonomy[, 'family'] %in% c('Unassigned', 'f__uncultured')
df.tax.filter <- df.taxonomy[!idx.del, ]
idx.del <- df.tax.filter[, 'genus'] %in% c('Unassigned', 'g__uncultured', 
                                           'g__uncultured_bacterium', 'Ambiguous_taxa')
df.tax.filter <- df.tax.filter[!idx.del, ]

# dose
list.vec.dose <- list(c(0, 1, 2, 3), c(0, 1), c(0, 2), c(0, 3))
for (vec.dose in list.vec.dose) {
    # dose
    # vec.dose <- c(0, 1, 2, 3)
    str_dose <- paste0(as.character(vec.dose), collapse = '')
    path.dose <- paste(pathout, str_dose, sep = '/')
    if (!file.exists(path.dose)) {
        dir.create(path.dose)
    }
    # gender
    vec.gender <- c('male', 'female')
    df.tree.merge <- data.frame(stringsAsFactors = F)
    for (gender in vec.gender) {
        if (gender == 'female') {
            sel.time <- c(1,  5,  9, 17, 21, 25, 29, 33, 41, 49, 60, 68, 84)
        } else {
            sel.time <- c(1, 5, 9, 17, 21, 25, 29, 33, 41)
        }
        df.meta.gender <- df.meta[df.meta$Gender == gender, ]
        # cutoff
        type.cutoff <- 'diff'
        
        # creat folder
        path.diff <- 
            paste0(path.dose, '/', gender, '_', type.cutoff, '_', abundance)
        if (!file.exists(path.diff)) {
            dir.create(path.diff)
        }
        series.time <- unique(df.meta$Time)

        # function of difference analysis
        diff.OTU <- function(df.meta.gender, df.OTUs, df.tax.filter, path.diff, vec.dose, 
                             type.cutoff, sub.time) {
            library(edgeR)
            sel.meta <- df.meta.gender[df.meta.gender$Time == sub.time,]
            sel.meta <- sel.meta[sel.meta$Dose %in% vec.dose,]
            use.sample <- sel.meta$SampleName
            
            sub.OTUs <- df.OTUs[, use.sample]
            
            factor.group <- as.factor(sel.meta$Group)
            d = DGEList(counts=sub.OTUs, group=factor.group)
            d = calcNormFactors(d)
            
            # 生成实验设计矩阵
            design.mat = model.matrix(~ 0 + d$samples$group)
            dimnames(design.mat)[[2]] <- levels(factor.group)
            d2 = estimateGLMCommonDisp(d, design.mat)
            d2 = estimateGLMTagwiseDisp(d2, design.mat)
            fit = glmFit(d2, design.mat)
            
            # 设置比较组
            BvsA <- makeContrasts(
                contrasts = paste(levels(factor.group), collapse = '-'),
                levels=design.mat)
            # 组间比较,统计Fold change, Pvalue
            lrt = glmLRT(fit,contrast=BvsA)
            # FDR检验，控制假阳性率小于5%
            # de_lrt = decideTestsDGE(lrt, adjust.method="fdr", p.value=0.05)
            
            # 导出计算结果
            res.edgeR=lrt$table
            # res.edgeR$sig.edger=de_lrt
            # vec.sig <- rep(0, dim(res.edgeR)[1])
            # vec.sig[(res.edgeR$logFC > 1) & (res.edgeR$PValue < 0.05)] <- 1
            # vec.sig[(res.edgeR$logFC < -1) & (res.edgeR$PValue < 0.05)] <- -1
            # res.edgeR$sig <- vec.sig
            # if (type.cutoff == 'fdr') {
            #     enriched = row.names(subset(res.edgeR, sig.edger==1))
            #     depleted = row.names(subset(res.edgeR, sig.edger==-1))
            # } else {
            #     enriched = row.names(subset(res.edgeR, sig==1))
            #     depleted = row.names(subset(res.edgeR, sig==-1))
            # }
            
            # write results
            # bool.fc <- rep(0, dim(res.edgeR)[1])
            # bool.fc[abs(res.edgeR$logFC) > 1.5] <- 1
            # res.edgeR$bool.fc <- bool.fc
            file.res <- paste0(path.diff, "/OUT_edgeR_", sub.time, 
                               paste0(as.character(vec.dose), collapse = ''), 
                               ".txt")
            res.edgeR$qvalue <- p.adjust(res.edgeR$PValue, method = 'fdr')
            # add taxonomy
            res.edgeR <- merge(res.edgeR, df.tax.filter, by = 'row.names')
            res.edgeR <- res.edgeR[order(res.edgeR$logFC), ]
            write.table(res.edgeR, file = file.res, quote = F, sep = '\t')
            
            return(file.res)
        }
        
        files.res <- foreach(sub.time = series.time, .combine = rbind) %dopar% 
            diff.OTU(df.meta.gender, df.OTUs, df.tax.filter, 
                     path.diff, vec.dose, type.cutoff, sub.time)
        
        # select significantly different OTUs
        for (i in 1:length(files.res)) {
            file <- files.res[i]
            df.res <- read.delim(file)
            df.res$log10pval <- df.res$logFC / abs(df.res$logFC) * log10(df.res$PValue)
            df.res$log10pval[is.na(df.res$log10pval)] <- 0
            df.res$log10pval[abs(df.res$log10pval) < 1] <- 0
            df.res <- df.res[, c('Row.names', 'log10pval')]
            names(df.res) <- c('OTU_id', series.time[i])
            if (i == 1) {
                df.combine <- df.res
            } else {
                df.combine <- merge(df.combine, df.res, by = 'OTU_id')
            }
        }
        
        df.combine$Sum <- rowSums(df.combine[, as.character(sel.time)])
        df.combine$absSum <- abs(df.combine$Sum)
        
        # OTUs.diff <- df.combine[(df.combine$Sum > 1) | (df.combine$Sum < -1), 'OTU_id']
        df.combine <- df.combine[order(df.combine$absSum, decreasing = T),]
        num.otu <- nrow(df.combine[df.combine$absSum > 0,])
        OTUs.diff <- df.combine[1:min(topN, num.otu), 'OTU_id']
        
        # output tree file
        df.tree <- df.tax.filter[as.character(OTUs.diff), c('phylum', 'class', 'order', 'family', 'genus', 'species')]
        df.tree$species <- str_replace_all(df.tree$species, fixed('.'), '_')
        df.tree$OTU_ID <- row.names(df.tree)
        write.table(df.tree, file = paste0(path.diff, '/tree_backbone_', topN, '.txt'), sep = '.', 
                    col.names = F, row.names = F, quote = F)
        df.tree.merge <- rbind(df.tree.merge, df.tree)
    }
    
    # merge 2 tree files
    length(unique(df.tree.merge$OTU_ID))
    df.tree.unique <- df.tree.merge[!duplicated(df.tree.merge$OTU_ID),]
    file.tree <- 
        paste0(path.dose, '/tree_backbone_', topN, '.txt')
    write.table(df.tree.unique, file = file.tree, sep = '.', 
                col.names = F, row.names = F, quote = F)
}
