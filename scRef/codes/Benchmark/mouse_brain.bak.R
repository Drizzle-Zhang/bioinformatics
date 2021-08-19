# import python package: sklearn.metrics
library(reticulate)
use_python('/home/zy/tools/anaconda3/bin/python3', required = T)
# py_config()
py_module_available('sklearn')
metrics <- import('sklearn.metrics')

# function of data preparation
prepare.data <- function(file.data.unlabeled, file.label.unlabeled, 
                         del.label = c('miss')) {
    library(stringr)
    data.unlabeled <- read.delim(file.data.unlabeled, row.names=1)
    data.unlabeled <- floor(data.unlabeled)
    names(data.unlabeled) <- str_replace_all(names(data.unlabeled), '_', '.')
    names(data.unlabeled) <- str_replace_all(names(data.unlabeled), '-', '.')
    # read label file
    file.label.unlabeled <- file.label.unlabeled
    label.unlabeled <- read.delim(file.label.unlabeled, row.names=1)
    row.names(label.unlabeled) <- str_replace_all(row.names(label.unlabeled), '_', '.')
    row.names(label.unlabeled) <- str_replace_all(row.names(label.unlabeled), '-', '.')
    col.name1 <- names(data.unlabeled)[1]
    if (substring(col.name1, 1, 1) == 'X') {
        row.names(label.unlabeled) <- paste0('X', row.names(label.unlabeled))
    }
    # filter data
    use.cols <- row.names(label.unlabeled)[!label.unlabeled[,1] %in% del.label]
    data.filter <- data.unlabeled[,use.cols]
    label.filter <- data.frame(label.unlabeled[use.cols,], row.names = use.cols)
    
    OUT <- list()
    OUT$data.filter <- data.filter
    OUT$label.filter <- label.filter
    return(OUT)
    
}

# evaluation
simple.evaluation <- function(true.tag, )

source('/home/zy/my_git/scRef/main/scRef.v10.R')

############# regard counts data as reference
path.input <- '/home/zy/scRef/summary/'
path.output <- '/home/zy/scRef/Benchmark/mouse_brain/'
dataset <- 'Tasic'
file.data.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat_cluster_merged.txt')
# OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, del.label = c('Unclassified'))
# saveRDS(OUT, file = paste0(path.output, dataset, '.Rdata'))
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
data.filter <- OUT$data.filter
label.filter <- OUT$label.filter
label.in <- data.frame(cell_id = row.names(label.filter), tag = label.filter$label.unlabeled.use.cols...)
exp.Tasic.sum <- .generate_ref(data.filter, label.in, M='SUM')
exp_ref_mat <- exp.Tasic.sum

############### import unlabeled data
############### Habib
dataset <- 'Habib'
file.data.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat_cluster_original.txt')
# OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, del.label = c('miss'))
# saveRDS(OUT, file = paste0(path.output, dataset, '.Rdata'))
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_sc_mat <- OUT$data.filter
label.origin <- OUT$label.filter

ref.names <- colnames(exp_ref_mat)
# list of cell names
all.cell <- unique(label.origin[,1])
sc.name <- c("Astrocyte", "EndothelialCells",
             "microglia", "Neurons", "Oligodend", "OPC")
unknow.cell <- setdiff(all.cell, sc.name)
df.cell.names <- data.frame(ref.name = ref.names, sc.name = sc.name, idx = 1:length(sc.name))


#################
### scRef
source('/home/zy/my_git/scRef/main/scRef.v10.R')
setwd('~/my_git/scRef')
result.scref <- SCREF(exp_sc_mat, exp_ref_mat, type_ref = 'count', 
                      cluster.speed = T, cluster.cell = 10,
                      min_cell = 10, CPU = 8)

library(ggplot2)
library("scales")
library(mclust)
sub.astrocyte <- df.tags[df.tags$scRef.tag.12 == 'Astrocyte', ]
sub.astrocyte <- df.tags1[df.tags1$scRef.tag == 'Astrocyte', ]
ggplot(sub.astrocyte, aes(x = log10Pval)) + geom_histogram(binwidth = 1)
model.astrocyte <- densityMclust(sub.astrocyte$log10Pval)
summary(model.astrocyte, parameters = T)
sub.astrocyte$cluster <- model.astrocyte$classification
Thickness <- sub.astrocyte$log10Pval
dens <- model.astrocyte
x <- seq(min(Thickness)-diff(range(Thickness))/10,max(Thickness)+diff(range(Thickness))/10, length = 2000)
cdens <- predict(dens, x, what = "cdens")
cdens <- t(apply(cdens, 1, function(d) d*dens$parameters$pro))
for (i in 1:dim(cdens)[2]) {
    if (i == 1) {
        df.dens <- data.frame(score = x, density = cdens[, 1], cluster = rep('1', 2000))
    } else {
        df.dens <- rbind(df.dens, data.frame(score = x, density = cdens[, i], cluster = rep(as.character(i), 2000)))
    }
}
# ggplot(df.dens, aes(x = score, y = density, color = cluster)) + geom_line()
plot.2 <- ggplot() + 
    geom_histogram(data = sub.astrocyte, aes(x = log10Pval), fill = 'gray', alpha = 1, binwidth = 1) + 
    geom_line(data = df.dens, aes(x = score, y = rescale(density, c(0, 250)), color = cluster)) + 
    geom_vline(xintercept = min(Thickness[dens$classification == '6'])) + 
    labs(x = 'Confidence Score', y = 'Count') + 
    scale_y_continuous(breaks=pretty_breaks(5),sec.axis = sec_axis( ~rescale(.,c(0,0.25)),name = "Density"))+
    theme(
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'transparent', color = 'gray'),
        legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)
        # legend.text = element_text(size = 15)
    )
ggsave(plot = plot.2, path = '/home/zy/scRef/figure', filename = 'GMM3.png',
       units = 'cm', height = 10, width = 12)

sub.neuron <- df.tags[df.tags$scRef.tag.12 == 'Neuron', ]
sub.neuron <- df.tags1[df.tags1$scRef.tag == 'Neuron', ]
ggplot(sub.neuron, aes(x = log10Pval)) + geom_histogram(binwidth = 1)
model.neuron <- densityMclust(sub.neuron$log10Pval, G=2)
summary(model.neuron, parameters = T)

sub.Endo <- df.tags1[df.tags1$scRef.tag == 'Endothelial Cell', ]
ggplot(sub.Endo, aes(x = log10Pval)) + geom_histogram(binwidth = 1)
model.Endo <- densityMclust(sub.Endo$log10Pval)
summary(model.Endo, parameters = T)
Thickness <- sub.Endo$log10Pval
dens <- model.Endo
x <- seq(min(Thickness)-diff(range(Thickness))/10,max(Thickness)+diff(range(Thickness))/10, length = 2000)
cdens <- predict(dens, x, what = "cdens")
cdens <- t(apply(cdens, 1, function(d) d*dens$parameters$pro))
for (i in 1:dim(cdens)[2]) {
    if (i == 1) {
        df.dens <- data.frame(score = x, density = cdens[, 1], cluster = rep('1', 2000))
    } else {
        df.dens <- rbind(df.dens, data.frame(score = x, density = cdens[, i], cluster = rep(as.character(i), 2000)))
    }
}
# ggplot(df.dens, aes(x = score, y = density, color = cluster)) + geom_line()
plot.2 <- ggplot() + 
    geom_histogram(data = sub.Endo, aes(x = log10Pval), fill = 'gray', alpha = 1, binwidth = 1) + 
    geom_line(data = df.dens, aes(x = score, y = rescale(density, c(0, 25)), color = cluster)) + 
    geom_vline(xintercept = min(Thickness[dens$classification == '2'])) + 
    labs(x = 'Confidence Score', y = 'Count') + 
    scale_y_continuous(breaks=pretty_breaks(5),sec.axis = sec_axis( ~rescale(.,c(0,0.25)),name = "Density"))+
    theme(
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'transparent', color = 'gray'),
        legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)
        # legend.text = element_text(size = 15)
    )
ggsave(plot = plot.2, path = '/home/zy/scRef/figure', filename = 'GMM5.png',
       units = 'cm', height = 10, width = 12)

sub.Microglia <- df.tags[df.tags$scRef.tag.12 == 'Microglia', ]
sub.Microglia <- df.tags1[df.tags1$scRef.tag == 'Microglia', ]
ggplot(sub.Microglia, aes(x = log10Pval)) + geom_histogram(binwidth = 1)
model.Microglia <- densityMclust(sub.Microglia$log10Pval, G=6)
summary(model.Microglia, parameters = T)

ref.oligo <- df.tags1[select.tag1[select.tag1$tag == 'Oligodend', 'cell_id'], ]
ggplot(ref.oligo, aes(x = log10Pval)) + geom_histogram(binwidth = 1)
sub.Oligo <- df.tags1[df.tags1$scRef.tag == 'Myelinating oligodendrocyte', ]
sub.Oligo <- df.tags[df.tags$scRef.tag.12 == 'Myelinating oligodendrocyte', ]
sub.Oligo <- df.tags1[df.tags1$scRef.tag == 'Oligodendrocyte', ]
sub.Oligo <- df.tags[df.tags$scRef.tag == 'Oligodendrocyte', ]
ggplot(sub.Oligo, aes(x = log10Pval)) + geom_histogram(binwidth = 1)
model.Oligo <- densityMclust(sub.Oligo$log10Pval)
summary(model.Oligo, parameters = T)
Thickness <- sub.Oligo$log10Pval
dens <- model.Oligo
x <- seq(min(Thickness)-diff(range(Thickness))/10,max(Thickness)+diff(range(Thickness))/10, length = 2000)
cdens <- predict(dens, x, what = "cdens")
cdens <- t(apply(cdens, 1, function(d) d*dens$parameters$pro))
for (i in 1:dim(cdens)[2]) {
    if (i == 1) {
        df.dens <- data.frame(score = x, density = cdens[, 1], cluster = rep('1', 2000))
    } else {
        df.dens <- rbind(df.dens, data.frame(score = x, density = cdens[, i], cluster = rep(as.character(i), 2000)))
    }
}
# ggplot(df.dens, aes(x = score, y = density, color = cluster)) + geom_line()
plot.1 <- ggplot() + 
    geom_histogram(data = sub.Oligo, aes(x = log10Pval), fill = 'gray', alpha = 1, binwidth = 3) + 
    geom_line(data = df.dens, aes(x = score, y = rescale(density, c(0, 50)), color = cluster)) + 
    geom_vline(xintercept = min(Thickness[dens$classification == '4'])) + 
    labs(x = 'Confidence Score', y = 'Count') + 
    scale_y_continuous(breaks=pretty_breaks(5),sec.axis = sec_axis( ~rescale(.,c(0,0.6)),name = "Density"))+
    theme(
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'transparent', color = 'gray'),
        legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)
        # legend.text = element_text(size = 15)
    )
ggsave(plot = plot.1, path = '/home/zy/scRef/figure', filename = 'GMM4.png',
       units = 'cm', height = 10, width = 12)

sub.OPC <- df.tags1[df.tags1$scRef.tag == 'Oligodendrocyte Precursor Cell', ]
sub.OPC <- df.tags1[df.tags1$scRef.tag == 'Oligodendrocyte precursor cell', ]
sub.OPC <- df.tags[df.tags$scRef.tag.12 == 'Oligodendrocyte precursor cell', ]
ggplot(sub.OPC, aes(x = log10Pval)) + geom_histogram(binwidth = 1)
model.OPC <- densityMclust(sub.OPC$log10Pval)
summary(model.OPC, parameters = T)
Thickness <- sub.OPC$log10Pval
dens <- model.OPC
x <- seq(min(Thickness)-diff(range(Thickness))/10,max(Thickness)+diff(range(Thickness))/10, length = 2000)
cdens <- predict(dens, x, what = "cdens")
cdens <- t(apply(cdens, 1, function(d) d*dens$parameters$pro))
for (i in 1:dim(cdens)[2]) {
    if (i == 1) {
        df.dens <- data.frame(score = x, density = cdens[, 1], cluster = rep('1', 2000))
    } else {
        df.dens <- rbind(df.dens, data.frame(score = x, density = cdens[, i], cluster = rep(as.character(i), 2000)))
    }
}
# ggplot(df.dens, aes(x = score, y = density, color = cluster)) + geom_line()
plot.1 <- ggplot() + 
    geom_histogram(data = sub.OPC, aes(x = log10Pval), fill = 'gray', alpha = 1, binwidth = 3) + 
    geom_line(data = df.dens, aes(x = score, y = rescale(density, c(0, 40)), color = cluster)) + 
    geom_vline(xintercept = min(Thickness[dens$classification == '3'])) + 
    labs(x = 'Confidence Score', y = 'Count') + 
    scale_y_continuous(breaks=pretty_breaks(5),sec.axis = sec_axis( ~rescale(.,c(0,0.50)),name = "Density"))+
    theme(
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'transparent', color = 'gray'),
        legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)
        # legend.text = element_text(size = 15)
    )
ggsave(plot = plot.1, path = '/home/zy/scRef/figure', filename = 'GMM5.png',
       units = 'cm', height = 10, width = 12)

sub.ependymal <- df.tags1[df.tags1$scRef.tag == 'Hypothalamic ependymal cell', ]
sub.ependymal <- df.tags[df.tags$scRef.tag.12 == 'Hypothalamic ependymal cell', ]
ggplot(sub.ependymal, aes(x = log10Pval)) + geom_histogram(binwidth = 1)
model.ependymal <- densityMclust(sub.ependymal$log10Pval)
summary(model.ependymal, parameters = T)
min(sub.ependymal$log10Pval[model.ependymal$classification == '4'])

meta.tag <- merge(result.scref$final.out, label.origin, by = 'row.names')
row.names(meta.tag) <- meta.tag$Row.names
meta.tag$Row.names <- NULL
names(meta.tag) <- c('scRef.tag', 'log10Pval', 'ori.tag')

### evaluation
true.tag = meta.tag$ori.tag
scRef.tag = meta.tag$scRef.tag

# uniform tags
for (j in 1:dim(df.cell.names)[1]) {
    scRef.tag[scRef.tag == df.cell.names[j, 'ref.name']] <- 
        df.cell.names[j, 'sc.name']
}
meta.tag$scRef.tag <- scRef.tag

# default cutoff
true.tag[true.tag %in% unknow.cell] <- 'Unassigned'
our.tag <- meta.tag$scRef.tag
metrics$f1_score(true.tag, our.tag, average = 'weighted')
metrics$f1_score(true.tag, our.tag, average = 'macro')
metrics$accuracy_score(true.tag, our.tag)

