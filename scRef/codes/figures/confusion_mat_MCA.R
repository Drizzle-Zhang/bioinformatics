library(ggplot2)
library(RColorBrewer)
library(gridExtra)

# heatmap
plot.heatmap <- function(ref.dataset, dataset, true.tags, path.res, path.fig, 
                         method, names.sc, names.ref) {
    rda.scRef <- paste0(path.res, ref.dataset, '_', dataset, '_', method, '.Rdata')
    pred.scRef <- readRDS(rda.scRef)
    pred.scRef[!(pred.scRef %in% names.ref)] <- 'Unassigned'
    mytable <- table(true.tags, pred.scRef)
    supp.names <- setdiff(names.ref, colnames(mytable))
    old.colnames <- colnames(mytable)
    mytable <- cbind(mytable, matrix(rep(0, length(supp.names)*length(names.sc)), 
                                     nrow = length(names.sc), ncol = length(supp.names)))
    colnames(mytable) <- c(old.colnames, supp.names)
    mydata <- data.frame(stringsAsFactors = F)
    table.true <- table(true.tags)
    for (label1 in rownames(mytable)) {
        row.sum <- table.true[label1]
        for (label2 in colnames(mytable)) {
            mydata <- rbind(mydata, data.frame(origin = label1, annotation = label2, 
                                               count = mytable[label1, label2], 
                                               prop = mytable[label1, label2]/row.sum))
        }
    }
    mydata$origin <- factor(mydata$origin, levels = names.sc)
    mydata$annotation <- factor(mydata$annotation, levels = names.ref)
    
    plot.heatmap <- 
        ggplot(data = mydata, aes(x = origin, y = annotation)) + 
        geom_tile(aes(fill = prop)) + 
        scale_fill_gradient2(low = "#FFF5EE", mid = '#EE7700', high = "#B22222", midpoint = 0.5) + 
        labs(fill = 'Proportion', title = method) + 
        # geom_text(aes(label = round(prop, 2)), family = "Arial", size = 2.5) +
        # labs(fill = 'Proportion', title = 'scMAGIC') +
        theme_bw() +
        theme(
            title = element_text(size = 12, color = "black", family = 'Arial'),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.title = element_blank(),
            plot.title = element_text(hjust = 0.5),
            axis.text.y = element_text(size = 10, color = "black", family = 'Arial'),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, 
                                       size = 10, color = "black", family = 'Arial')
        )
    ggsave(filename = paste0('heatmap_', ref.dataset, '_', dataset, '_', method, '.png'), 
           path = path.fig, plot = plot.heatmap,
           units = 'cm', height = 10, width = 15)
    return(plot.heatmap)
}

path.res <- '/home/zy/scRef/Benchmark/mouse_brain/'
path.fig <- '/home/zy/scRef/figure/MCA/'
if (!file.exists(path.fig)) {
    dir.create(path.fig)
}

ref.dataset <- 'MCA'
dataset <- 'Campbell'
OUT <- readRDS(paste0(path.res, dataset, '.Rdata'))
label_sc <- OUT$label
true.tags <- label_sc[,1]
names.ref <- c('Neuron', 'Myelinating oligodendrocyte', 'Oligodendrocyte precursor cell', 
               'Microglia', 'Macrophage', 'Astrocyte', 'Hypothalamic ependymal cell',
               'Unassigned', 'Astroglial cell', 'Granulocyte', 'Schwann cell')
names.sc <- c('Neurons', 'Oligodendrocytes', 'OPC', 'PVMs & Microglia', 
              'Astrocytes', 'Ependymocytes', 
              'VLMCs','Endothelial cells',  'Mural cells', 'Pars tuberalis', 'Tanycytes')

method <- 'scMAGIC'
plot.heatmap(ref.dataset, dataset, true.tags, path.res, path.fig, method, names.sc, names.ref)

method <- 'sciBet'
plot.heatmap(ref.dataset, dataset, true.tags, path.res, path.fig, method, names.sc, names.ref)

method <- 'scPred'
plot.heatmap(ref.dataset, dataset, true.tags, path.res, path.fig, method, names.sc, names.ref)

# method <- 'singleCellNet'
# plot.heatmap(ref.dataset, dataset, true.tags, path.res, path.fig, method, names.sc, names.ref)

method <- 'scClassify'
plot.heatmap(ref.dataset, dataset, true.tags, path.res, path.fig, method, names.sc, names.ref)

