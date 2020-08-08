library(ape)
library(ggtree)

# read OTUs matrix
file.OTUs <- '/home/drizzle_zhang/microbiome/result/2.OTUs/OTUs_stat/OTUs_tax_even.csv'
df.OTUs <- read.csv(file.OTUs, row.names = 1)
df.OTUs <- df.OTUs[, 1:dim(df.OTUs)[2] - 1]
# choose topN OTUs
df.OTUs.mean <- rowSums(df.OTUs)
topN <- 500
df.OTUs.topN <- df.OTUs[order(df.OTUs.mean, decreasing = T)[1:topN],]
OTUs.topN <- row.names(df.OTUs.topN)

# tree including all OTUs
file.tree <- '/home/drizzle_zhang/microbiome/result/2.OTUs/rep_phylo.tre'
tree <- read.tree(file = file.tree)
tree.topN <- tree
index.use <- which((tree$edge[,1] %in% OTUs.topN) & (tree$edge[,2] %in% OTUs.topN))
tree.topN$edge <- tree$edge[index.use,]
tree.topN$edge.length <- tree$edge.length[index.use]
tree.topN$

file.tree.top100 <- '/home/drizzle_zhang/microbiome/result/8.Phylogenetic/top100.tre'
tree.top100 <- read.tree(file = file.tree.top100)

plot.tree <- ggtree(tree.top100, layout = 'circular')

