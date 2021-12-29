library(ggpubr)
library(stringr)


df_organ <- read.delim2('/home/drizzle_zhang/microbiome/experiment_data/OrganWeight.txt',
                        stringsAsFactors = F)
df_organ$RelativeWeight <- as.numeric(df_organ$RelativeWeight)

df_organ$Dosage <- str_replace_all(df_organ$Dosage, fixed('0'), '0.0')
df_organ$Dosage <- str_replace_all(df_organ$Dosage, fixed('6'), '6.0')
df_organ$Dosage <- str_replace_all(df_organ$Dosage, fixed('0.0.375'), '0.375')

mycomparisons <- list(c("0.0", "0.375"), c("0.0", "1.5"), c("0.0", "6.0"))
plot_organ <- 
    ggboxplot(df_organ, x = "Dosage", y = "RelativeWeight", color = "Dosage", 
              add = c("dotplot", "boxplot"), add.params = list(fill = 'transparent'),
              facet.by = c("Gender", "Organ"), outlier.shape = NA) + 
    stat_compare_means(comparisons = mycomparisons, size = 2.5, step.increase = 0.17)
plot_organ <-
    facet(plot_organ, facet.by = c("Gender", "Organ")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.11))) + 
    scale_color_manual(breaks = c('0.0', '0.375', '1.5', '6.0'),
                       values = c('dimgrey', 'firebrick', '#33A02C', '#1F78B4'), 
                       labels = c('0.0', '0.375', '1.5', '6.0')) +
    labs(x = 'Time (week)', y = 'Absolute Weights (mg)', color = 'Dosage (mg/kg IAA)') +
    coord_cartesian(ylim = c(0.006, 0.08)) + 
    theme_bw() + 
    theme(panel.background = element_rect(color = 'black',
                                          fill = 'transparent'),
          strip.background = element_rect(
              color = 'black', fill = 'transparent'),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(colour = 'black'),
          axis.ticks.x = element_blank(),
          legend.title = element_text(size = 9, family = 'Arial'),
          legend.text = element_text(size = 8, family = 'Arial'),
          panel.grid = element_blank())

file.organ <- '/home/drizzle_zhang/microbiome/result/Figs/Organ.png'
ggsave(filename = file.organ, plot = plot_organ, units = 'cm', height = 10, width = 14)
