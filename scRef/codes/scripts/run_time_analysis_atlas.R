library(ggplot2)
path.res <- '/mdshare/node9/zy/scRef/figure/atlas_anno/'

df.time <- data.frame(method = c(rep('scMAGIC', 4), rep('scHCL/scMCA', 4)),
                      dataset = rep(c('MCA -> Mouse Neocortex', 'MCA -> Mouse Duodenum', 
                                      'HCL -> Human Pancreas', 'HCL -> Human PBMC'), 2),
                      time = c(10.46, 8.72, 11.42, 7.98,
                               133.19, 79.86, 120.18, 39.72))

plot.bar <-
    ggplot(df.time, aes(x = dataset, y = time, color = method, fill = method)) +
    geom_bar(position = 'dodge', stat = 'identity') +
    scale_color_manual(breaks = c('scMAGIC', 'scHCL/scMCA'), values = c("#BC80BD", "#80B1D3")) +
    scale_fill_manual(breaks = c('scMAGIC', 'scHCL/scMCA'), values = c("#BC80BD", "#80B1D3")) +
    labs(title = "", y = 'Runtime (min)', x = '', color = '', fill = '') +
    coord_flip() +
    theme_bw() +
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 8, color = 'black', family = 'Arial'),
          axis.text.y = element_text(size = 10, color = 'black', family = 'Arial'),
          axis.title = element_text(size = 10, family = 'Arial'),
          legend.text = element_text(size = 10, family = 'Arial'),
          legend.position = 'right',
          # legend.key.size = unit(1, 'cm'),
          legend.key = element_blank())
ggsave(filename = 'RunTime_atlas.png',
       path = path.res, plot = plot.bar,
       units = 'cm', height = 6, width = 14)

