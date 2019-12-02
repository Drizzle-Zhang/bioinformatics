setwd('/home/drizzle_zhang/driver_mutation/cRE_plot')
library(ggplot2)


### flank cutoff
# DHS
df.flank.DHS <- read.delim('./DHS/flank_count.txt', 
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
         title = "Curve of DHS number with increasing of flank region(Terms Level)") +
    scale_colour_discrete(
        name = 'Legend',
        breaks = c('brain_embryonic_brain', 'frontal-cortex_adult_brain',
                   'sigmoid-colon_adult_intestine', 'limb_embryonic_limb',
                   'lung_embryonic_lung', 'cerebellar-cortex_adult_brain'),
        labels = c(
            'Embryonic brain',
            'Adult brain-frontal cortex',
            'Adult intestine-sigmoid colon',
            'Embryonic limb',
            'Embryonic lung',
            'Adult brain-cerebellar cortex'
        )
    ) +
    scale_shape_manual(
        values = c(0, 1, 2, 4, 20, 5),
        name = 'Legend',
        breaks = c('brain_embryonic_brain', 'frontal-cortex_adult_brain',
                   'sigmoid-colon_adult_intestine', 'limb_embryonic_limb',
                   'lung_embryonic_lung', 'cerebellar-cortex_adult_brain'),
        labels = c(
            'Embryonic brain',
            'Adult brain-frontal cortex',
            'Adult intestine-sigmoid colon',
            'Embryonic limb',
            'Embryonic lung',
            'Adult brain-cerebellar cortex'
        )
    )
    # theme(panel.background = element_rect(fill = "transparent",colour = NA),
    #       plot.background = element_rect(fill = "transparent",colour = NA),
    #       panel.grid.major = element_line(colour = 'grey'),
    #       panel.grid.minor = element_line(colour = 'grey'),
    #       # axis.line = element_line(colour = "black")
    #       )

ggsave(
    plot = DHS.plot, path = './DHS', filename = "flank_DHS_term.png",
    units = 'cm', width = 25, height = 15)
df.flank.DHS[df.flank.DHS$flank_percent == 0, ]
df.flank.DHS[df.flank.DHS$flank_percent == 0.5, ]


# H3K4me3
df.flank.H3K4me3 <- read.delim('./H3K4me3/flank_count.txt', 
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
         title = "Curve of H3K4me3 Peak number with increasing of flank region(Terms level)") +
    scale_colour_discrete(
        name = 'Legend',
        breaks = c('stomach_adult_stomach', 'layer-of-hippocampus_adult_brain',
                   'sigmoid-colon_adult_intestine', 
                   'adrenal-gland_adult_adrenal-gland',
                   'placenta_embryonic_extraembryonic-component',
                   'heart-left-ventricle_adult_heart'),
        labels = c(
            'Adult stomach',
            'Adult brain-layer of hippocampus',
            'Adult intestine-sigmoid colon',
            'Adult adrenal gland',
            'Embryonic extraembryonic component-placenta',
            'Adult Heart-heart left ventricle'
        )
    ) +
    scale_shape_manual(
        values = c(0, 1, 2, 4, 20, 5),
        name = 'Legend',
        breaks = c('stomach_adult_stomach', 'layer-of-hippocampus_adult_brain',
                   'sigmoid-colon_adult_intestine', 
                   'adrenal-gland_adult_adrenal-gland',
                   'placenta_embryonic_extraembryonic-component',
                   'heart-left-ventricle_adult_heart'),
        labels = c(
            'Adult stomach',
            'Adult brain-layer of hippocampus',
            'Adult intestine-sigmoid colon',
            'Adult adrenal gland',
            'Embryonic extraembryonic component-placenta',
            'Adult Heart-heart left ventricle'
        )
    )
ggsave(
    plot = H3K4me3.plot, path = './H3K4me3', filename = "flank_H3K4me3_term.png",
    units = 'cm', width = 25, height = 15)
df.flank.H3K4me3[df.flank.H3K4me3$flank_percent == 0, ]
df.flank.H3K4me3[df.flank.H3K4me3$flank_percent == 0.2, ]


# H3K27ac
df.flank.H3K27ac <- read.delim('./H3K27ac/flank_count.txt', 
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
         title = "Curve of H3K27ac Peak number with increasing of flank region(Terms Level)") +
    scale_colour_discrete(
        name = 'Legend',
        breaks = c('stomach_adult_stomach', 'layer-of-hippocampus_adult_brain',
                   'sigmoid-colon_adult_intestine', 
                   'adrenal-gland_adult_adrenal-gland',
                   'placenta_embryonic_extraembryonic-component',
                   'heart-left-ventricle_adult_heart'),
        labels = c(
            'Adult stomach',
            'Adult brain-layer of hippocampus',
            'Adult intestine-sigmoid colon',
            'Adult adrenal gland',
            'Embryonic extraembryonic component-placenta',
            'Adult Heart-heart left ventricle'
        )
    ) +
    scale_shape_manual(
        values = c(0, 1, 2, 4, 20, 5),
        name = 'Legend',
        breaks = c('stomach_adult_stomach', 'layer-of-hippocampus_adult_brain',
                   'sigmoid-colon_adult_intestine', 
                   'adrenal-gland_adult_adrenal-gland',
                   'placenta_embryonic_extraembryonic-component',
                   'heart-left-ventricle_adult_heart'),
        labels = c(
            'Adult stomach',
            'Adult brain-layer of hippocampus',
            'Adult intestine-sigmoid colon',
            'Adult adrenal gland',
            'Embryonic extraembryonic component-placenta',
            'Adult Heart-heart left ventricle'
        )
    )
ggsave(
    plot = H3K27ac.plot, path = './H3K27ac', filename = "flank_H3K27ac_term.png",
    units = 'cm', width = 25, height = 15)
df.flank.H3K27ac[df.flank.H3K27ac$flank_percent == 0, ]
df.flank.H3K27ac[df.flank.H3K27ac$flank_percent == 0.5, ]
