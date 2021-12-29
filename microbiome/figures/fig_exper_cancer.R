library(ggpubr)

df_cancer <- 
    data.frame(Gender = c(rep('Male', 8), rep('Female', 8)), 
               Cancer = c(rep('Thyroid Follicular Cell\nAdenoma or Carcinoma', 4),
                          rep('Pituitary Adenoma', 4), 
                          rep('Thyroid Follicular Cell\nAdenoma or Carcinoma', 4),
                          rep('Pituitary Adenoma', 4)),
               Dosage = rep(c('0.0', '0.375', '1.5', '6.0'), 4),
               Rate = c(4/37, 1/36, 4/35, 6/36, 11/37, 15/36, 15/35, 19/36,
                        1/40, 1/39, 4/38, 8/38, 10/40, 23/39, 27/38, 30/38),
               stringsAsFactors = F)

plot_cancer <-  
    ggbarplot(df_cancer, x = "Dosage", y = "Rate", color = "Dosage", 
              facet.by = c("Gender", "Cancer"))
plot_cancer <-
    facet(plot_cancer, facet.by = c("Gender", "Cancer"), scales = 'free') +
    # scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
    scale_color_manual(breaks = c('0.0', '0.375', '1.5', '6.0'),
                       values = c('dimgrey', 'firebrick', '#33A02C', '#1F78B4'), 
                       labels = c('0.0', '0.375', '1.5', '6.0')) +
    labs(x = '', y = 'Incidence rate of neoplasms lesions', 
         color = 'Dosage (mg/kg IAA)') +
    coord_cartesian(ylim = c(0, 0.9)) + 
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

file.cancer <- '/home/drizzle_zhang/microbiome/result/Figs/CancerRate.png'
ggsave(filename = file.cancer, plot = plot_cancer, units = 'cm', 
       height = 10.4, width = 14)

df_52 <- 
    data.frame(Gender = c(rep('Male', 8), rep('Female', 8)), 
               Cancer = c(rep('Pituitary Adenoma', 4), rep('No pituitary Adenoma', 4), 
                          rep('Pituitary Adenoma', 4), rep('No pituitary Adenoma', 4)),
               Dosage = rep(c('0.0', '0.375', '1.5', '6.0'), 4),
               Num = c(0, 1, 2, 1, 5, 4, 3, 4,
                        0, 0, 0, 0, 5, 5, 5, 5),
               stringsAsFactors = F)
df_52$Cancer <- factor(df_52$Cancer, levels = c('Pituitary Adenoma', 'No pituitary Adenoma'))
plot_52 <-  
    ggbarplot(df_52, x = "Dosage", y = "Num", color = "Dosage", 
              facet.by = c("Gender", "Cancer"))
plot_52 <-
    facet(plot_52, facet.by = c("Gender", "Cancer"), scales = 'free') +
    # scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
    scale_color_manual(breaks = c('0.0', '0.375', '1.5', '6.0'),
                       values = c('dimgrey', 'firebrick', '#33A02C', '#1F78B4'), 
                       labels = c('0.0', '0.375', '1.5', '6.0')) +
    labs(x = '', y = 'Number of rats', 
         color = 'Dosage (mg/kg IAA)') +
    coord_cartesian(ylim = c(0, 5.5)) + 
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
file.52 <- '/home/drizzle_zhang/microbiome/result/Figs/Cancer52.png'
ggsave(filename = file.52, plot = plot_52, units = 'cm', 
       height = 10, width = 14)



