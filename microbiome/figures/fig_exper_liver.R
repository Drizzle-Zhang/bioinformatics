library(ggpubr)
library(Rmisc)
library(cowplot)

df_liver <- read.delim2('/home/drizzle_zhang/microbiome/experiment_data/LiverFunction.txt',
                        stringsAsFactors = F)
df_liver$varience <- as.character(df_liver$varience)
df_liver$Gender <- unlist(lapply(df_liver$varience, 
                                  function(x) {strsplit(x, split = '-', fixed = T)[[1]][2]}))
df_liver$Index <- unlist(lapply(df_liver$varience, 
                                 function(x) {strsplit(x, split = '-', fixed = T)[[1]][1]}))
df_liver$value <- as.numeric(df_liver$value)
df_liver$Group <- toupper(df_liver$group)

mycomparisons <- list(c("A", "B"), c("A", "C"), c("A", "D"))
plot_liver <- 
    ggbarplot(df_liver, x = "Group", y = "value", color = "Group", 
          add = c("mean_se", "dotplot"), 
          palette = "jco", facet.by = c("Index", "Gender")) + 
    stat_compare_means(comparisons = mycomparisons, size = 2.5, step.increase = 0.17)
plot_liver <-
    facet(plot_liver, facet.by = c("Index", "Gender"), scales = 'free')+
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
    scale_color_manual(breaks = c('A', 'B', 'C', 'D'),
                       values = c('dimgrey', 'firebrick', '#33A02C', '#1F78B4'), 
                       labels = c('0.0', '0.375', '1.5', '6.0')) +
    labs(x = 'Time (week)', y = 'Concentrations in Serum (U/L)', color = 'Dosage (mg/kg IAA)') +
    theme_bw() + 
    theme(panel.background = element_rect(color = 'black',
                                          fill = 'transparent'),
          strip.background = element_rect(
              color = 'black', fill = 'transparent'),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(colour = 'black'),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank())


file.liver <- '/home/drizzle_zhang/microbiome/result/Figs/Liver.png'
ggsave(filename = file.liver, plot = plot_liver, units = 'cm', height = 15, width = 15)
