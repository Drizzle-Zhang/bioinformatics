setwd('~/driver_mutation/cRE_plot/DHS')
library(ggplot2)

df.overlap.exp <- read.delim('./overlap_exp.txt', sep = '\t', 
                             stringsAsFactors = F)
df.overlap.exp$Level <- rep('Experiment', dim(df.overlap.exp)[1])

df.overlap.term <- read.delim('./overlap_term.txt', sep = '\t', 
                              stringsAsFactors = F)
df.overlap.term$Level <- rep('Term', dim(df.overlap.term)[1])

df.overlap.suborgan <- read.delim('./overlap.suborgan.txt', sep = '\t', 
                                  stringsAsFactors = F)
df.overlap.suborgan$Level <- rep('Suborgan', dim(df.overlap.suborgan)[1])

df.overlap.organ <- read.delim('./overlap_organ.txt', sep = '\t', 
                               stringsAsFactors = F)
df.overlap.organ <- df.overlap.organ[
    !(df.overlap.organ$Combination %in% df.overlap.suborgan$Combination),]
df.overlap.organ$Level <- rep('Organ', dim(df.overlap.organ)[1])

df.overlap.all <- read.delim('./all_organs_overlap.txt', sep = '\t', 
                             stringsAsFactors = F)
df.overlap.all$Level <- rep('All', dim(df.overlap.all)[1])

col_select <- c('Jaccard.distance', 'Level')
df.plot <- rbind(
    df.overlap.exp[,col_select], df.overlap.term[,col_select], 
    df.overlap.suborgan[,col_select], df.overlap.organ[,col_select], 
    df.overlap.all[,col_select]
)


