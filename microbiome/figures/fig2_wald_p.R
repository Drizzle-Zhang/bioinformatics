library("AER")
# parameters
# filter by abundance
abundance <- 0.0001
# path
pathout <- '/home/drizzle_zhang/microbiome/result/6.Phylogenetic/graphlan'

# read Rdata
setwd(pathout)
df.OTUs <- readRDS(paste0('OTUs_filtered_', abundance, '.Rdata'))
df.taxonomy <- readRDS(paste0('taxonomy_filtered_', abundance, '.Rdata'))

# difference analysis
# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.out.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)

# delete unclassified families
idx.del <- df.taxonomy[, 'family'] %in% c('Unassigned', 'f__uncultured')
df.tax.filter <- df.taxonomy[!idx.del, ]
idx.del <- df.tax.filter[, 'genus'] %in% c('Unassigned', 'g__uncultured', 
                                           'g__uncultured_bacterium', 'Ambiguous_taxa')
df.tax.filter <- df.tax.filter[!idx.del, ]

# dose
vec.dose <- c(0, 1, 2, 3)
str_dose <- paste0(as.character(vec.dose), collapse = '')
path.dose <- paste(pathout, str_dose, sep = '/')

# OTU
sel.OTU <- c('13', '55', '122', '73', '128', '50190', '143', '100', 
             '38', '512', '205')
series.time <- c(-1, 1,  5,  9, 17, 21, 25, 29, 33, 41, 49, 60, 68, 84)

# gender
vec.gender <- c('male', 'female')
for (gender in vec.gender) {
    print(gender)
    df.meta.gender <- df.meta[df.meta$Gender == gender, ]
    sel.meta <- df.meta.gender[df.meta.gender$Dose %in% vec.dose,]
    sel.meta <- sel.meta[sel.meta$Time %in% series.time,]
    use.sample <- sel.meta$SampleName
    rownames(sel.meta) <- use.sample
    
    sub.OTUs <- t(df.OTUs[sel.OTU, use.sample])
    col.OTUs <- paste('OTU', colnames(sub.OTUs), sep = '')
    colnames(sub.OTUs) <- col.OTUs
    sub.OTUs <- cbind(sub.OTUs, sel.meta[, c('Time', 'Group')])
    
    # wald
    for (OTU in col.OTUs) {
        sub.OTU <- sub.OTUs[, c(OTU, 'Time', 'Group')]
        sub.OTU$Time <- factor(sub.OTU$Time)
        sub.OTU$Group <- factor(sub.OTU$Group)
        colnames(sub.OTU) <- c('OTU', 'Time', 'Group')
        sub.OTU <- as.data.frame(sub.OTU)
        print(OTU)
        fm1 <- lm(OTU ~ Time, data = sub.OTU)
        fm2 <- lm(OTU ~ Time + Group, data = sub.OTU)
        print(waldtest(fm2, fm1))
        sub.fit <- aov(OTU ~ Time*Group, data = sub.OTU)
        # print(summary(sub.fit)[[1]]$`Pr(>F)`)
        print(summary(sub.fit))
    }
}
