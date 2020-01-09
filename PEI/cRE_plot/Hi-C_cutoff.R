setwd('/home/drizzle_zhang/driver_mutation/cRE_plot')
library(ggplot2)


### cutoff
# pcHi-C
df.pcHiC <- read.delim('./HiC/pcHiCresult.txt', 
                       sep = '\t', stringsAsFactors = F)
df.num.pcHiC <- data.frame()
col.names <- c("count", "num_egene", "num_egene_pairs", "num_overlap_egene",
               "num_pairs", "overlap_num", "total")
for (col.name in col.names) {
    df.sub <- df.pcHiC[, c('label', col.name)]
    names(df.sub) <- c('Cutoff', 'Number')
    df.sub$Label <- rep(col.name, dim(df.pcHiC)[1])
    df.num.pcHiC <- rbind(df.num.pcHiC, df.sub)
}

# proportion
effective.interation <- df.pcHiC$ratio
overlap.egene <- df.pcHiC$overlap_num/df.pcHiC$num_egene_pairs
overlap.eqtl <- df.pcHiC$overlap_num/5747
df.ratio.pcHiC <- data.frame(
    Cutoff = df.pcHiC$label, Proportion = effective.interation,
    Label = rep('effective.interation', dim(df.pcHiC)[1]))
df.ratio.pcHiC <- rbind(df.ratio.pcHiC, data.frame(
    Cutoff = df.pcHiC$label, Proportion = overlap.egene,
    Label = rep('overlap.egene', dim(df.pcHiC)[1])))
df.ratio.pcHiC <- rbind(df.ratio.pcHiC, data.frame(
    Cutoff = df.pcHiC$label, Proportion = overlap.eqtl,
    Label = rep('overlap.eqtl', dim(df.pcHiC)[1])))

pcHiC.plot <-
    ggplot(
        df.num.pcHiC,
        aes(
            x = Cutoff,
            y = Number,
            color = Label,
            shape = Label
        )
    ) +
    geom_line() +
    geom_point(size = 1) +
    labs(x = 'Cutoff', y = 'Number', title = "Promoter Capture Hi-C") +
    scale_colour_discrete(
        name = 'Legend',
        breaks = c('total', 'count', 'num_pairs', 'num_egene_pairs',
                   'num_egene', 'overlap_num', 'num_overlap_egene'),
        labels = c(
            'pcHi-C interations', 'pcHi-C interations annotated by cREs', 
            'Gene-cRE pairs', 'eGene-cRE pairs', 
            'Expressed genes in eGene-cRE pairs', 
            'eGene-cRE pairs overlapped with eQTL', 
            'Expressed genes in eGene-cRE pairs overlapped with eQTL'
        )
    ) +
    scale_shape_manual(
        values = c(0, 1, 2, 4, 20, 5, 7),
        name = 'Legend',
        breaks = c('total', 'count', 'num_pairs', 'num_egene_pairs',
                   'num_egene', 'overlap_num', 'num_overlap_egene'),
        labels = c(
            'pcHi-C interations', 'pcHi-C interations annotated by cREs', 
            'Gene-cRE pairs', 'eGene-cRE pairs', 
            'Expressed genes in eGene-cRE pairs', 
            'eGene-cRE pairs overlapped with eQTL', 
            'Expressed genes in eGene-cRE pairs overlapped with eQTL'
        )
    )

ggsave(
    plot = pcHiC.plot, path = './HiC', filename = "pcHi-C.png",
    units = 'cm', width = 25, height = 15)


df.num.pcHiC$Number <- log10(df.num.pcHiC$Number)

pcHiC.plot <-
    ggplot(
        df.num.pcHiC,
        aes(
            x = Cutoff,
            y = Number,
            color = Label,
            shape = Label
        )
    ) +
    geom_line() +
    geom_point(size = 1) +
    labs(x = 'Cutoff', y = 'log10(Number)', title = "Promoter Capture Hi-C") +
    scale_colour_discrete(
        name = 'Legend',
        breaks = c('total', 'count', 'num_pairs', 'num_egene_pairs',
                   'num_egene', 'overlap_num', 'num_overlap_egene'),
        labels = c(
            'pcHi-C interations', 'pcHi-C interations annotated by cREs', 
            'Gene-cRE pairs', 'eGene-cRE pairs', 
            'Expressed genes in eGene-cRE pairs', 
            'eGene-cRE pairs overlapped with eQTL', 
            'Expressed genes in eGene-cRE pairs overlapped with eQTL'
        )
    ) +
    scale_shape_manual(
        values = c(0, 1, 2, 4, 20, 5, 7),
        name = 'Legend',
        breaks = c('total', 'count', 'num_pairs', 'num_egene_pairs',
                   'num_egene', 'overlap_num', 'num_overlap_egene'),
        labels = c(
            'pcHi-C interations', 'pcHi-C interations annotated by cREs', 
            'Gene-cRE pairs', 'eGene-cRE pairs', 
            'Expressed genes in eGene-cRE pairs', 
            'eGene-cRE pairs overlapped with eQTL', 
            'Expressed genes in eGene-cRE pairs overlapped with eQTL'
        )
    )

ggsave(
    plot = pcHiC.plot, path = './HiC', filename = "pcHi-C_log.png",
    units = 'cm', width = 25, height = 15)


pcHiC.plot <-
    ggplot(
        df.ratio.pcHiC,
        aes(
            x = Cutoff,
            y = Proportion,
            color = Label,
            shape = Label
        )
    ) +
    geom_line() +
    geom_point(size = 1) +
    labs(x = 'Cutoff', y = 'Proportion', title = "Promoter Capture Hi-C") +
    scale_colour_discrete(
        name = 'Legend',
        breaks = c('effective.interation', 'overlap.egene', 'overlap.eqtl'),
        labels = c(
            'pcHi-C interations annotated by cREs / pcHi-C interations', 
            'eGene-cRE pairs overlapped with eQTL / eGene-cRE pairs from pcHi-C', 
            'eGene-cRE pairs overlapped with eQTL / eGene-cRE pairs from eQTL'
        )
    ) +
    scale_shape_manual(
        values = c(0, 1, 2),
        name = 'Legend',
        breaks = c('effective.interation', 'overlap.egene', 'overlap.eqtl'),
        labels = c(
            'pcHi-C interations annotated by cREs / pcHi-C interations', 
            'eGene-cRE pairs overlapped with eQTL / eGene-cRE pairs from pcHi-C', 
            'eGene-cRE pairs overlapped with eQTL / eGene-cRE pairs from eQTL'
        )
    )


ggsave(
    plot = pcHiC.plot, path = './HiC', filename = "pcHi-C_proportion.png",
    units = 'cm', width = 25, height = 15)



#################################################################
# 3DIV
df.3DIV <- read.delim('./HiC/result.3DIV.txt', 
                       sep = '\t', stringsAsFactors = F)
df.num.3DIV <- data.frame()
col.names <- c("count", "num_egene", "num_egene_pairs", "num_overlap_egene",
               "num_pairs", "overlap_num", "total")
for (col.name in col.names) {
    df.sub <- df.3DIV[, c('label', col.name)]
    names(df.sub) <- c('Cutoff', 'Number')
    df.sub$Label <- rep(col.name, dim(df.3DIV)[1])
    df.num.3DIV <- rbind(df.num.3DIV, df.sub)
}

# proportion
effective.interation <- df.3DIV$ratio
overlap.egene <- df.3DIV$overlap_num/df.3DIV$num_egene_pairs
overlap.eqtl <- df.3DIV$overlap_num/5747
df.ratio.3DIV <- data.frame(
    Cutoff = df.3DIV$label, Proportion = effective.interation,
    Label = rep('effective.interation', dim(df.3DIV)[1]))
df.ratio.3DIV <- rbind(df.ratio.3DIV, data.frame(
    Cutoff = df.3DIV$label, Proportion = overlap.egene,
    Label = rep('overlap.egene', dim(df.3DIV)[1])))
df.ratio.3DIV <- rbind(df.ratio.3DIV, data.frame(
    Cutoff = df.3DIV$label, Proportion = overlap.eqtl,
    Label = rep('overlap.eqtl', dim(df.3DIV)[1])))

plot.3DIV <-
    ggplot(
        df.num.3DIV,
        aes(
            x = Cutoff,
            y = Number,
            color = Label,
            shape = Label
        )
    ) +
    geom_line() +
    geom_point(size = 1) +
    labs(x = 'Cutoff', y = 'Number', title = "Hi-C from 3DIV") +
    scale_colour_discrete(
        name = 'Legend',
        breaks = c('total', 'count', 'num_pairs', 'num_egene_pairs',
                   'num_egene', 'overlap_num', 'num_overlap_egene'),
        labels = c(
            'pcHi-C interations', 'pcHi-C interations annotated by cREs', 
            'Gene-cRE pairs', 'eGene-cRE pairs', 
            'Expressed genes in eGene-cRE pairs', 
            'eGene-cRE pairs overlapped with eQTL', 
            'Expressed genes in eGene-cRE pairs overlapped with eQTL'
        )
    ) +
    scale_shape_manual(
        values = c(0, 1, 2, 4, 20, 5, 7),
        name = 'Legend',
        breaks = c('total', 'count', 'num_pairs', 'num_egene_pairs',
                   'num_egene', 'overlap_num', 'num_overlap_egene'),
        labels = c(
            'pcHi-C interations', 'pcHi-C interations annotated by cREs', 
            'Gene-cRE pairs', 'eGene-cRE pairs', 
            'Expressed genes in eGene-cRE pairs', 
            'eGene-cRE pairs overlapped with eQTL', 
            'Expressed genes in eGene-cRE pairs overlapped with eQTL'
        )
    )

ggsave(
    plot = plot.3DIV, path = './HiC', filename = "HiC3DIV.png",
    units = 'cm', width = 25, height = 15)


df.num.3DIV$Number <- log10(df.num.3DIV$Number)

plot.3DIV <-
    ggplot(
        df.num.3DIV,
        aes(
            x = Cutoff,
            y = Number,
            color = Label,
            shape = Label
        )
    ) +
    geom_line() +
    geom_point(size = 1) +
    labs(x = 'Cutoff', y = 'log10(Number)', title = "Hi-C from 3DIV") +
    scale_colour_discrete(
        name = 'Legend',
        breaks = c('total', 'count', 'num_pairs', 'num_egene_pairs',
                   'num_egene', 'overlap_num', 'num_overlap_egene'),
        labels = c(
            'pcHi-C interations', 'pcHi-C interations annotated by cREs', 
            'Gene-cRE pairs', 'eGene-cRE pairs', 
            'Expressed genes in eGene-cRE pairs', 
            'eGene-cRE pairs overlapped with eQTL', 
            'Expressed genes in eGene-cRE pairs overlapped with eQTL'
        )
    ) +
    scale_shape_manual(
        values = c(0, 1, 2, 4, 20, 5, 7),
        name = 'Legend',
        breaks = c('total', 'count', 'num_pairs', 'num_egene_pairs',
                   'num_egene', 'overlap_num', 'num_overlap_egene'),
        labels = c(
            'pcHi-C interations', 'pcHi-C interations annotated by cREs', 
            'Gene-cRE pairs', 'eGene-cRE pairs', 
            'Expressed genes in eGene-cRE pairs', 
            'eGene-cRE pairs overlapped with eQTL', 
            'Expressed genes in eGene-cRE pairs overlapped with eQTL'
        )
    )

ggsave(
    plot = plot.3DIV, path = './HiC', filename = "log.3DIV.png",
    units = 'cm', width = 25, height = 15)


plot.3DIV <-
    ggplot(
        df.ratio.3DIV,
        aes(
            x = Cutoff,
            y = Proportion,
            color = Label,
            shape = Label
        )
    ) +
    geom_line() +
    geom_point(size = 1) +
    labs(x = 'Cutoff', y = 'Proportion', title = "Hi-C from 3DIV") +
    scale_colour_discrete(
        name = 'Legend',
        breaks = c('effective.interation', 'overlap.egene', 'overlap.eqtl'),
        labels = c(
            'pcHi-C interations annotated by cREs / pcHi-C interations', 
            'eGene-cRE pairs overlapped with eQTL / eGene-cRE pairs from pcHi-C', 
            'eGene-cRE pairs overlapped with eQTL / eGene-cRE pairs from eQTL'
        )
    ) +
    scale_shape_manual(
        values = c(0, 1, 2),
        name = 'Legend',
        breaks = c('effective.interation', 'overlap.egene', 'overlap.eqtl'),
        labels = c(
            'pcHi-C interations annotated by cREs / pcHi-C interations', 
            'eGene-cRE pairs overlapped with eQTL / eGene-cRE pairs from pcHi-C', 
            'eGene-cRE pairs overlapped with eQTL / eGene-cRE pairs from eQTL'
        )
    )


ggsave(
    plot = plot.3DIV, path = './HiC', filename = "proportion.3DIV.png",
    units = 'cm', width = 25, height = 15)

