library(ggplot2)
library(gridExtra)

ref.dataset <- c('Tasic', 'MCA')
ref.GSE_id <- c('GSE71585', 'GSE108097')
ref.legend <- c('Tasic, Primary visual cortex, FACS+SMARTer', 
                'Han, Brain, Microwell-seq')

sc.dataset <- c('Tasic2018', 'Campbell', 'HochgernerA', 'Mizrak')
sc.GSE_id <- c('GSE115746', 'GSE93374', 'GSE95315', 'GSE109447')
sc.legends <- c('Tasic, Neocortex, SMART-Seq v4',
               'Campbell, hypothalamic Arc-ME, Drop-seq',
               'Hochgerner, Denatate gyrus, 10X',
               'Mizrak, V-SVZ, Drop-seq')
# sc.legends <- c('Tasic (SMART-Seq v4)', 'Campbell (Drop-seq)', 
#                 'Hochgerner (10X)', 'Mizrak (Drop-seq)')

path.res <- '/home/zy/scRef/Benchmark/mouse_brain/'
path.fig <- '/home/zy/scRef/figure/mouse_brain'

# Accuracy
for (i in 1:length(ref.dataset)) {
    ref.data <- ref.dataset[i]
    ref.GSE <- ref.GSE_id[i]
    df.plot <- data.frame(stringsAsFactors = F)
    for (j in 1:length(sc.dataset)) {
        sc.data <- sc.dataset[j]
        sc.GSE <- sc.GSE_id[j]
        sc.legend <- sc.legends[j]
        file.res <- paste0(path.res, 'results_', ref.data, '_', sc.data, '.txt')
        sub.res <- read.delim(file = file.res, row.names = 1, stringsAsFactors = F)
        sub.res <- sub.res[sub.res$term == 'Accuracy',]
        acc.scRef <- sub.res[sub.res$method == 'scMAGIC', 'value']
        sub.plot <- data.frame(other = sub.res[sub.res$method != 'scMAGIC', 'value'], 
                               scRef = rep(acc.scRef, nrow(sub.res)-1), 
                               method = sub.res[sub.res$method != 'scMAGIC', 'method'],
                               dataset = rep(sc.legend, nrow(sub.res)-1))
        df.plot <- rbind(df.plot, sub.plot)
    }
    plot.scatter <- 
        ggplot(data = df.plot, aes(x = other, y = scRef, color = method, shape = dataset)) + 
        geom_point(size = 3) + 
        geom_abline(intercept = 0, slope = 1, linetype = 2) + 
        ylim(0, 1) + 
        labs(x = 'Accuracy(Other method)', y = 'Accuracy(scMAGIC)', 
             shape = 'Query dataset', color = 'Method') + 
        theme(panel.background = element_rect(fill = 'transparent', color = 'gray'))
    ggsave(filename = paste0('scatter_', ref.data, '.png'), 
           path = path.fig, plot = plot.scatter,
           units = 'cm', height = 12, width = 18)
}

# Macro F1
for (i in 1:length(ref.dataset)) {
    ref.data <- ref.dataset[i]
    ref.GSE <- ref.GSE_id[i]
    df.plot <- data.frame(stringsAsFactors = F)
    for (j in 1:length(sc.dataset)) {
        sc.data <- sc.dataset[j]
        sc.GSE <- sc.GSE_id[j]
        sc.legend <- sc.legends[j]
        file.res <- paste0(path.res, 'results_', ref.data, '_', sc.data, '.txt')
        sub.res <- read.delim(file = file.res, row.names = 1, stringsAsFactors = F)
        sub.res <- sub.res[sub.res$term == 'Macro F1',]
        acc.scRef <- sub.res[sub.res$method == 'scMAGIC', 'value']
        sub.plot <- data.frame(other = sub.res[sub.res$method != 'scMAGIC', 'value'], 
                               scRef = rep(acc.scRef, nrow(sub.res)-1), 
                               method = sub.res[sub.res$method != 'scMAGIC', 'method'],
                               dataset = rep(sc.legend, nrow(sub.res)-1))
        df.plot <- rbind(df.plot, sub.plot)
    }
    plot.scatter <- 
        ggplot(data = df.plot, aes(x = other, y = scRef, color = method, shape = dataset)) + 
        geom_point(size = 3) + 
        geom_abline(intercept = 0, slope = 1, linetype = 2) + 
        ylim(0, 1) + xlim(0, 1) + 
        labs(x = 'Macro F1(Other method)', y = 'Macro F1(scMAGIC)', 
             shape = 'Query dataset', color = 'Method') + 
        theme(panel.background = element_rect(fill = 'transparent', color = 'gray'))
    ggsave(filename = paste0('scatter_MacroF1_', ref.data, '.png'), 
           path = path.fig, plot = plot.scatter,
           units = 'cm', height = 12, width = 18)
}

# Accuracy (rm)
for (i in 1:length(ref.dataset)) {
    ref.data <- ref.dataset[i]
    ref.GSE <- ref.GSE_id[i]
    df.plot <- data.frame(stringsAsFactors = F)
    for (j in 1:length(sc.dataset)) {
        sc.data <- sc.dataset[j]
        sc.GSE <- sc.GSE_id[j]
        sc.legend <- sc.legends[j]
        file.res <- paste0(path.res, 'results_', ref.data, '_', sc.data, '.txt')
        sub.res <- read.delim(file = file.res, row.names = 1, stringsAsFactors = F)
        sub.res <- sub.res[sub.res$term == 'Accuracy (remove unassigned)',]
        acc.scRef <- sub.res[sub.res$method == 'scMAGIC', 'value']
        sub.plot <- data.frame(other = sub.res[sub.res$method != 'scMAGIC', 'value'], 
                               scRef = rep(acc.scRef, nrow(sub.res)-1), 
                               method = sub.res[sub.res$method != 'scMAGIC', 'method'],
                               dataset = rep(sc.legend, nrow(sub.res)-1))
        df.plot <- rbind(df.plot, sub.plot)
    }
    plot.scatter <- 
        ggplot(data = df.plot, aes(x = other, y = scRef, color = method, shape = dataset)) + 
        geom_point(size = 3) + 
        geom_abline(intercept = 0, slope = 1, linetype = 2) + 
        ylim(0, 1) + 
        labs(x = 'Accuracy(Other method, remove unassigned cells)', 
             y = 'Accuracy(scMAGIC, remove unassigned cells)', 
             shape = 'Query dataset', color = 'Method') + 
        theme(panel.background = element_rect(fill = 'transparent', color = 'gray'))
    ggsave(filename = paste0('scatter_Accuracy_rm_', ref.data, '.png'), 
           path = path.fig, plot = plot.scatter,
           units = 'cm', height = 12, width = 18)
}


# heatmap
plot.heatmap <- function(ref.dataset, dataset, true.tags, path.res, path.fig, method, names.sc, names.ref) {
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
        scale_fill_continuous(low = "#FFFAFA", high = "#CD5C5C") + 
        labs(fill = 'Proportion', title = method) +
        # labs(fill = 'Proportion', title = 'scMAGIC') +
        theme_bw() +
        theme(
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.title = element_blank(),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
        ) + 
        geom_text(aes(label = round(prop, 2)), family = "Arial", size = 2.5)
    ggsave(filename = paste0('heatmap_', ref.dataset, '_', dataset, '_', method, '.png'), 
           path = path.fig, plot = plot.heatmap,
           units = 'cm', height = 12, width = 22)
}

ref.dataset <- 'Tasic'
dataset <- 'Campbell'
OUT <- readRDS(paste0(path.res, dataset, '.Rdata'))
label_sc <- OUT$label
true.tags <- label_sc[,1]
names.ref <- c('Neuron', 'Oligodendrocyte', 'Oligodendrocyte Precursor Cell', 
               'Microglia', 'Endothelial Cell', 'Astrocyte', 'Unassigned')
names.sc <- c('Neurons', 'Oligodendrocytes', 'OPC', 'PVMs & Microglia', 
              'Endothelial cells', 'Astrocytes', 
              'Ependymocytes', 'VLMCs', 'Mural cells', 'Pars tuberalis', 'Tanycytes')

ref.dataset <- 'MCA'
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

method <- 'singleCellNet'
plot.heatmap(ref.dataset, dataset, true.tags, path.res, path.fig, method, names.sc, names.ref)

method <- 'scClassify'
plot.heatmap(ref.dataset, dataset, true.tags, path.res, path.fig, method, names.sc, names.ref)

