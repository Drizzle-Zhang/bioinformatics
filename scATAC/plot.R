library(ggplot2)

df_plot <- data.frame(dataset = c('PBMC', 'Forebrain'), Accuracy = c(0.89, 0.91))

plot_bar <- 
    ggplot(df_plot, aes(x = dataset, y = Accuracy)) + 
    geom_bar(stat = 'identity', color = 'dimgray', fill = 'dimgray', 
             position=position_dodge(0.6)) + 
    coord_cartesian(ylim = c(0, 1)) + 
    # coord_flip() + 
    labs(x = 'Dataset', y = 'Accuracy') +
    theme_bw() +
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 11, color = 'black', family = 'Arial'),
          axis.text.y = element_text(size = 8, color = 'black', family = 'Arial'),
          axis.title = element_text(size = 12, family = 'Arial'))
ggsave(filename = 'acc_plot.png',
       path = '/mdshare/node9/zy/scATAC', plot = plot_bar,
       units = 'cm', height = 10, width = 7)
