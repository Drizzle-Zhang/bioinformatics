setwd('/home/drizzle_zhang/driver_mutation/cRE_plot')
library(ggplot2)


### flank cutoff
# DHS
df.flank.DHS <- read.delim('./DHS/flank_count_exp.txt', 
                           sep = '\t', stringsAsFactors = F)

DHS.plot <-
    ggplot(
        df.flank.DHS,
        aes(
            x = flank_percent,
            y = num_rows,
            color = label,
            shape = label
        )
    ) +
    geom_line() +
    geom_point(size = 1) +
    labs(x = 'Extended Flank Percent', y = 'DHS Numbers',
         title = "Curve of DHS number with increasing of flank region(Experiments Level)") +
    scale_colour_discrete(
        name = 'Legend'
    ) +
    scale_shape_manual(
        values = c(0, 1, 2, 4, 20, 5),
        name = 'Legend'
    )
# theme(panel.background = element_rect(fill = "transparent",colour = NA),
#       plot.background = element_rect(fill = "transparent",colour = NA),
#       panel.grid.major = element_line(colour = 'grey'),
#       panel.grid.minor = element_line(colour = 'grey'),
#       # axis.line = element_line(colour = "black")
#       )

ggsave(
    plot = DHS.plot, path = './DHS', filename = "flank_DHS_exp.png",
    units = 'cm', width = 25, height = 15)
df.flank.DHS[df.flank.DHS$flank_percent == 0, ]
df.flank.DHS[df.flank.DHS$flank_percent == 0.56, ]


# H3K4me3
df.flank.H3K4me3 <- read.delim('./H3K4me3/flank_count_exp.txt', 
                               sep = '\t', stringsAsFactors = F)

H3K4me3.plot <-
    ggplot(
        df.flank.H3K4me3,
        aes(
            x = flank_percent,
            y = num_rows,
            color = label,
            shape = label
        )
    ) +
    geom_line() +
    geom_point(size = 1) +
    labs(x = 'Extended Flank Percent', y = 'H3K4me3 Peak Numbers',
         title = "Curve of H3K4me3 Peak number with increasing of flank region(Experiments level)") +
    scale_colour_discrete(
        name = 'Legend'
    ) +
    scale_shape_manual(
        values = c(0, 1, 2, 4, 20, 5),
        name = 'Legend'
    )
ggsave(
    plot = H3K4me3.plot, path = './H3K4me3', filename = "flank_H3K4me3_exp.png",
    units = 'cm', width = 25, height = 15)
df.flank.H3K4me3[df.flank.H3K4me3$flank_percent == 0, ]
df.flank.H3K4me3[df.flank.H3K4me3$flank_percent == 0.12, ]


# H3K27ac
df.flank.H3K27ac <- read.delim('./H3K27ac/flank_count_exp.txt', 
                               sep = '\t', stringsAsFactors = F)

H3K27ac.plot <-
    ggplot(
        df.flank.H3K27ac,
        aes(
            x = flank_percent,
            y = num_rows,
            color = label,
            shape = label
        )
    ) +
    geom_line() +
    geom_point(size = 1) +
    labs(x = 'Extended Flank Percent', y = 'H3K27ac Peak Numbers',
         title = "Curve of H3K27ac Peak number with increasing of flank region(Experiments Level)") +
    scale_colour_discrete(
        name = 'Legend'
    ) +
    scale_shape_manual(
        values = c(0, 1, 2, 4, 20, 5),
        name = 'Legend')
ggsave(
    plot = H3K27ac.plot, path = './H3K27ac', filename = "flank_H3K27ac_exp.png",
    units = 'cm', width = 25, height = 15)
df.flank.H3K27ac[df.flank.H3K27ac$flank_percent == 0, ]
df.flank.H3K27ac[df.flank.H3K27ac$flank_percent == 0.4, ]
