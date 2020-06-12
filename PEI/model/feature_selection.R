library(ggplot2)
library(ggsignif)

setwd('/local/zy/PEI/mid_data/cell_line/model_input/GM12878')
setwd('/lustre/tianlab/zhangyu/PEI/mid_data_correct/cell_line/model_input/GM12878')
# setwd('/home/drizzle_zhang/driver_mutation/cRE_plot/model_test')
file.in <- './selection.txt'
df.mat <- read.delim(file.in, stringsAsFactors = F)
df.mat$label <- as.character(df.mat$label)

unique(df.mat$label)
compaired <- list(c("0", "1"))

aggre <- aggregate(df.mat$DHS_DHS, by = list(df.mat$label), FUN = median)
aggre$x[2] - aggre$x[1]
my.plot <- 
    ggplot(df.mat, aes(label, DHS_DHS, fill = label)) +
    geom_boxplot(width = 0.5) +
    theme(
        plot.title = element_text(size = 25),
        axis.text.x = element_text(size = 15, angle = 0),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 23),
        axis.title.y = element_text(size = 23)
    ) + labs(x = 'Label', y = 'DHS_DHS') + geom_signif(
        comparisons = compaired,
        step_increase = 0.1,
        map_signif_level = F,
        test = 'wilcox.test'
    )
ggsave('./DHS_DHS.pdf', my.plot)


aggre <- aggregate(df.mat$expression_DHS, by = list(df.mat$label), FUN = median)
aggre$x[2] - aggre$x[1]
my.plot <- 
    ggplot(df.mat, aes(label, expression_DHS, fill = label)) +
    geom_boxplot(width = 0.5) +
    theme(
        plot.title = element_text(size = 25),
        axis.text.x = element_text(size = 15, angle = 0),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 23),
        axis.title.y = element_text(size = 23)
    ) + labs(x = 'Label', y = 'expression_DHS') + geom_signif(
        comparisons = compaired,
        step_increase = 0.1,
        map_signif_level = F,
        test = 'wilcox.test'
    )
ggsave('./expression_DHS.pdf', my.plot)

# z-score
# assay
file.assay <- './Assay_comparation_1_0.txt'
df.compare <- read.delim(file.assay, stringsAsFactors = F)
ggplot(df.compare, aes(x = feature, y = diff_median, fill = correlation)) + 
    geom_boxplot() + 
    labs(x = 'Feature', y = 'Diff of median', fill = 'Methods of correlation')

df.compare$pval[df.compare$pval == 0] <- 10^-300
df.compare$logp <- -log10(df.compare$pval)
ggplot(df.compare, aes(x = feature, y = logp, fill = correlation)) + 
    geom_boxplot() + 
    labs(x = 'Feature', y = '-log(p-value)', fill = 'Methods of correlation')

# Cell
file.cell <- './Cell_comparation_1_0.txt'
df.compare <- read.delim(file.cell, stringsAsFactors = F)
ggplot(df.compare, aes(x = feature, y = diff_median, 
                       color = correlation, shape = label)) + 
    geom_point(size = 4) + 
    labs(x = 'Feature', y = 'Diff of median', color = 'Methods of correlation',
         shape = 'Cell line')

df.compare$pval[df.compare$pval == 0] <- 10^-300
df.compare$logp <- -log10(df.compare$pval)
ggplot(df.compare, aes(x = feature, y = logp, 
                       fill = correlation, shape = label)) + 
    geom_boxplot() + 
    labs(x = 'Feature', y = '-log(p-value)', fill = 'Methods of correlation')


# quantile
# assay
file.assay <- './Assay_comparation_1_0.txt'
df.compare <- read.delim(file.assay, stringsAsFactors = F)
ggplot(df.compare, aes(x = feature, y = diff_median, fill = correlation)) + 
    geom_boxplot() + 
    labs(x = 'Feature', y = 'Diff of median', fill = 'Methods of correlation')

df.compare$pval[df.compare$pval == 0] <- 10^-300
df.compare$logp <- -log10(df.compare$pval)
ggplot(df.compare, aes(x = feature, y = logp, fill = correlation)) + 
    geom_boxplot() + 
    labs(x = 'Feature', y = '-log(p-value)', fill = 'Methods of correlation')

# Cell
file.cell <- './Cell_comparation_1_0_quantile.txt'
df.compare <- read.delim(file.cell, stringsAsFactors = F)
ggplot(df.compare, aes(x = feature, y = diff_median, 
                       color = correlation, shape = label)) + 
    geom_point(size = 4) + 
    labs(x = 'Feature', y = 'Diff of median', color = 'Methods of correlation',
         shape = 'Cell line')

df.compare$pval[df.compare$pval == 0] <- 10^-300
df.compare$logp <- -log10(df.compare$pval)
ggplot(df.compare, aes(x = feature, y = logp, 
                       color = correlation, shape = label)) + 
    geom_point(size = 4) + 
    labs(x = 'Feature', y = '-log(p-value)', color = 'Methods of correlation',
         shape = 'Cell line')








