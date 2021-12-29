source('/home/zy/my_git/scRef/main/scRef.v5.R')

# calculate overlap
overlap <- function(list.genes1, list.genes2) {
    for (i in 1:length(ref.names)) {
        cell1 <- names(list.genes1)[i]
        cell2 <- names(list.genes2)[i]
        print(cell1)
        print(cell2)
        overlap <- length(intersect(list.genes1[[cell1]], list.genes2[[cell2]]))
        print(overlap)
    }
}

# input file
setwd('/home/zy/scRef/try_data')
file.ref <- './scRef/Reference/MouseBrain_Bulk_Zhang2014/Reference_expression.txt'
exp_ref_mat <- read.table(file.ref, header=T, row.names=1, sep='\t', check.name=F)

# combat auto
out.markers <- find.markers.auto(exp_ref_mat, type = 'fpkm', topN = 100)
auto.cell.genes <- out.markers[['list.cell.genes']]
# genes.ref <- dimnames(out.markers[['exp_ref_mat']])[[1]]

# combat manual match
ref.names <- c("Astrocyte", "Neuron", "Oligodendrocyte precursor cell",
               "Oligodendrocyte", "Myelinating oligodendrocyte",
               "Microglia", "Endothelial cell")
names(exp_ref_mat) <- ref.names
out.markers <- find.markers(exp_ref_mat, type = 'fpkm', topN = 100)
match.cell.genes <- out.markers[['list.cell.genes']]

# combat manual partial match
ref.names <- c("astrocyte", "Neuron", "Oligodendrocyte precursor cell",
               "Newly Formed Oligodendrocyte", "Myelinating oligodendrocyte",
               "Microglia", "Endothelial cell")
names(exp_ref_mat) <- ref.names
out.markers <- find.markers(exp_ref_mat, type = 'fpkm', topN = 100)
partial.cell.genes <- out.markers[['list.cell.genes']]

# HK ratio
ref.names <- c("Astrocyte", "Neuron", "Oligodendrocyte precursor cell",
               "Oligodendrocyte", "Myelinating oligodendrocyte",
               "Microglia", "Endothelial cell")
names(exp_ref_mat) <- ref.names
out.markers <- find.markers.HKratio(exp_ref_mat, topN = 100)
ratio.cell.genes <- out.markers[['list.cell.genes']]

# overlap compare
overlap(match.cell.genes, partial.cell.genes)
overlap(match.cell.genes, ratio.cell.genes)


