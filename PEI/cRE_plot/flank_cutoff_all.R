setwd('/home/drizzle_zhang/driver_mutation/cRE_plot')
library(ggplot2)


### flank cutoff
# DHS
df.flank.DHS <- read.delim('./DHS/flank_count_all.txt', 
                           sep = '\t', stringsAsFactors = F)

DHS.plot <-
    ggplot(
        df.flank.DHS,
        aes(
            x = flank_percent,
            y = num_rows
        )
    ) +
    geom_line() +
    geom_point(size = 1) +
    labs(x = 'Extended Flank Percent', y = 'DHS Numbers',
         title = "Curve of DHS number with increasing of flank region(All data)")

ggsave(
    plot = DHS.plot, path = './DHS', filename = "flank_DHS_all.png",
    units = 'cm', width = 25, height = 15)
df.flank.DHS[df.flank.DHS$flank_percent == 0, ]
df.flank.DHS[df.flank.DHS$flank_percent == 0.7000000000000001, ]
df.flank.DHS[df.flank.DHS$flank_percent == 0.8, ]


