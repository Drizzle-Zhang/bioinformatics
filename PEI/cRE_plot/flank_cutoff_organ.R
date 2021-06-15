setwd('/home/drizzle_zhang/driver_mutation/cRE_plot')
library(ggplot2)


### flank cutoff
# DHS
df.flank.DHS <- read.delim('./DHS/flank_count_organ.txt', 
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
         title = "Curve of DHS number with increasing of flank region(Organs Level)") +
    scale_colour_discrete(
        name = 'Legend',
        breaks = c('adrenal_gland', 'brain', 'intestine', 'musculature_of_body',
                   'testis', 'heart'),
        labels = c(
            'Adrenal gland', 'Brain', 'Intestine', 'Musculature of body',
            'Testis', 'Heart'
        )
    ) +
    scale_shape_manual(
        values = c(0, 1, 2, 4, 20, 5),
        name = 'Legend',
        breaks = c('adrenal_gland', 'brain', 'intestine', 'musculature_of_body',
                   'testis', 'heart'),
        labels = c(
            'Adrenal gland', 'Brain', 'Intestine', 'Musculature of body',
            'Testis', 'Heart'
        )
    )
# theme(panel.background = element_rect(fill = "transparent",colour = NA),
#       plot.background = element_rect(fill = "transparent",colour = NA),
#       panel.grid.major = element_line(colour = 'grey'),
#       panel.grid.minor = element_line(colour = 'grey'),
#       # axis.line = element_line(colour = "black")
#       )

ggsave(
    plot = DHS.plot, path = './DHS', filename = "flank_DHS_organ.png",
    units = 'cm', width = 25, height = 15)
df.flank.DHS[df.flank.DHS$flank_percent == 0, ]
df.flank.DHS[df.flank.DHS$flank_percent == 0.6, ]