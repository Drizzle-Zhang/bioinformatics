library(RColorBrewer)
library(foreach)
library(doParallel)
registerDoParallel(4)

# parameters
# filter by abundance
abundance <- 0.0001
# filter by top num
topN <- 500

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
setwd('/home/drizzle_zhang/microbiome/result/6.Phylogenetic/graphlan')
saveRDS(filtered_otutab, file = paste0('OTUs_filtered_', abundance, '.Rdata'))
saveRDS(filtered_taxonomy, file = paste0('taxonomy_filtered_', abundance, '.Rdata'))

# read Rdata
topN <- 500
abundance <- 0.0001
setwd('/home/drizzle_zhang/microbiome/result/6.Phylogenetic/graphlan')
df.OTUs <- readRDS(paste0('OTUs_filtered_', abundance, '.Rdata'))
df.taxonomy <- readRDS(paste0('taxonomy_filtered_', abundance, '.Rdata'))

# difference analysis
# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.out.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)

# delete unclassified families
idx.del <- df.taxonomy[, 'family'] %in% c('Unassigned', 'f__uncultured')
df.tax.filter <- df.taxonomy[!idx.del, ]
idx.del <- df.tax.filter[, 'genus'] %in% c('Unassigned', 'g__uncultured', 'g__uncultured_bacterium')
df.tax.filter <- df.tax.filter[!idx.del, ]

# gender
# gender = 'female'
# sel.time <- c(1,  5,  9, 17, 21, 25, 29, 33, 41, 49, 60, 68, 84)
gender = 'male'
sel.time <- c(1,  5,  9, 17, 21, 25, 29, 33, 41)
df.meta.gender <- df.meta[df.meta$Gender == gender, ]
# dose
# vec.dose <- c(0, 1, 2, 3)
vec.dose <- c(0, 3)
# cutoff
type.cutoff <- 'diff'

# creat folder
path.diff <- 
    paste0('/home/drizzle_zhang/microbiome/result/6.Phylogenetic/graphlan/',
           gender, '_', type.cutoff, '_', abundance)
if (!file.exists(path.diff)) {
    dir.create(path.diff)
}
# series.time <- unique(df.meta$Time)
series.time <- c(-1, 1,  5,  9, 17, 21, 25, 29, 33, 41, 49, 60, 68, 84)

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
OTUs.diff <- df.combine[1:100, 'OTU_id']

# output tree file
df.tree <- df.tax.filter[as.character(OTUs.diff), c('phylum', 'class', 'order', 'family', 'genus')]
df.tree$OTU_ID <- row.names(df.tree)
write.table(df.tree, file = 'tree_backbone.txt', sep = '.', 
            col.names = F, row.names = F, quote = F)

# set color for each phylum
table.phylum <- sort(table(df.tree[, 'phylum']), decreasing = T)
phylums <- names(table.phylum)
df.color <- data.frame(stringsAsFactors = F)
set.color <- brewer.pal(8, "Set2")
for (i in 1:length(phylums)) {
    if (i > 6 | table.phylum[phylums[i]] < 3) {
        sub.color <- '#D3D3D3'
    } else {
        sub.color <- set.color[i]
    }
    df.color <- rbind(df.color, 
                      data.frame(phylum = phylums[i], color = sub.color,
                                 stringsAsFactors = F))
}
df.color <- merge(unique(df.tree[, c('phylum', 'genus')]),
                  df.color, by = 'phylum')

# generate phylum annotation file
phylum_color <- data.frame(stringsAsFactors = F)
for (sub.phylum in phylums) {
    sub.anno <- data.frame(stringsAsFactors = F)
    
    sub.anno[1, 1] <- sub.phylum
    sub.anno[1, 2] <- 'annotation_background_color'
    sub.anno[1, 3] <- unique(df.color[df.color$phylum == sub.phylum, 'color'])
    
    sub.anno[2, 1] <- paste0(sub.phylum, '*')
    sub.anno[2, 2] <- 'clade_marker_color'
    sub.anno[2, 3] <- unique(df.color[df.color$phylum == sub.phylum, 'color'])
    
    phylum_color <- rbind(phylum_color, sub.anno)
    
}

# generate genus annotation file
families <- unique(df.tree[, 'genus'])
genus_color <- data.frame(stringsAsFactors = F)
for (sub.genus in families) {
    # name.genus <- substr(sub.genus, 4, nchar(sub.genus))
    sub.anno <- data.frame(stringsAsFactors = F)
    sub.anno[1, 1] <- sub.genus
    sub.anno[1, 2] <- 'annotation'
    sub.anno[1, 3] <- '*'
    # rotate text 90 degree
    sub.anno[2, 1] <- sub.genus
    sub.anno[2, 2] <- 'annotation_rotation'
    sub.anno[2, 3] <- '90'
    # set background color
    sub.anno[3, 1] <- sub.genus
    sub.anno[3, 2] <- 'annotation_background_color'
    sub.anno[3, 3] <- df.color[df.color$genus == sub.genus, 'color']
    
    genus_color <- rbind(genus_color, sub.anno)
    
}

label_color <- rbind(phylum_color, genus_color)
write.table(label_color, paste0('tree_label_color_', topN, '.txt'), 
            sep = '\t', quote = F, row.names = F, col.names = F, na = '')

# generate ring annotation file
limit <- 3
# color.pos <- '#F8766D'
# color.neg <- '#00BFC4'
color.pos <- '#EE4000'
color.neg <- '#20B2AA'
df.ring <- data.frame(stringsAsFactors = F)
df.ring.global <- data.frame(stringsAsFactors = F)
for (i in 1:length(series.time)) {
    sub.time <- series.time[i]
    
    time.anno <- data.frame(stringsAsFactors = F)
    time.anno[1, 1] <- 'ring_internal_separator_thickness'
    time.anno[1, 2] <- i
    time.anno[1, 3] <- 0.1
    time.anno[2, 1] <- 'ring_width'
    time.anno[2, 2] <- i
    time.anno[2, 3] <- 1
    time.anno[3, 1] <- 'ring_height'
    time.anno[3, 2] <- i
    time.anno[3, 3] <- 1
    time.anno[4, 1] <- 'ring_label'
    time.anno[4, 2] <- i
    time.anno[4, 3] <- sub.time
    time.anno[5, 1] <- 'ring_label_font_size'
    time.anno[5, 2] <- i
    time.anno[5, 3] <- 10
    df.ring.global <- rbind(df.ring.global, time.anno)
    
    file.res <- paste0(path.diff, "/OUT_edgeR_", sub.time, 
                       paste0(as.character(vec.dose), collapse = ''), 
                       ".txt")
    res.edgeR <- read.delim(file.res, row.names = 1)
    res.edgeR <- res.edgeR[res.edgeR$Row.names %in% OTUs.diff,]
    res.edgeR$logPval <- log10(res.edgeR$PValue) * 
        (res.edgeR$logFC / abs(res.edgeR$logFC))
    res.edgeR$logPval[is.na(res.edgeR$logPval)] <- 0
    for (j in 1:dim(res.edgeR)[1]) {
        sub.anno <- data.frame(stringsAsFactors = F)
        OTU_id <- res.edgeR[j, 'Row.names']
        sub.logPval <- res.edgeR[j, "logPval"]
        sub.anno[1, 1] <- OTU_id
        sub.anno[1, 2] <- 'ring_color'
        sub.anno[1, 3] <- i
        if (sub.logPval > 0) {
            sub.anno[1, 4] <- color.pos
        } else {
            sub.anno[1, 4] <- color.neg
        }
        sub.anno[2, 1] <- OTU_id
        sub.anno[2, 2] <- 'ring_alpha'
        sub.anno[2, 3] <- i
        sub.anno[2, 4] <- min(1, abs(sub.logPval/limit))
        df.ring <- rbind(df.ring, sub.anno)
    }
}

write.table(df.ring.global, paste0(path.diff, '/tree_ring_global_', topN, '.txt'), 
            sep = '\t', quote = F, row.names = F, col.names = F, na = '')
write.table(df.ring, paste0(path.diff, '/tree_ring_', topN, '.txt'), 
            sep = '\t', quote = F, row.names = F, col.names = F, na = '')


# shell 
system(paste0('cd /home/drizzle_zhang/microbiome/result/8.Phylogenetic/graphlan'))

cd /home/drizzle_zhang/microbiome/result/6.Phylogenetic/graphlan/female_fdr_500
cd /home/drizzle_zhang/microbiome/result/6.Phylogenetic/graphlan/male_fdr_500
rm -rf track*
cat ../cfg/global.cfg ../tree_label_color_500.txt > track0
cat tree_ring_global_500.txt tree_ring_500.txt > track1
cat track* > graphlan_annotate.txt
graphlan_annotate.py --annot graphlan_annotate.txt ../tree_backbone.txt graphlan.xml
# 绘图，size决定图片大小，越大字越小
graphlan.py graphlan.xml graphlan500_tree.pdf --size 8


