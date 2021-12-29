library(ggplot2)
# color bar
df_plot <- data.frame(pval = c(0, 1, 2, 3, 4), TF = as.character(c(0, 1, 2, 3, 4)))
df_plot$NUL <- rep('1', nrow(df_plot))

plot_bar_up <- 
    ggplot(data = df_plot, aes(x = TF, y = NUL, fill = pval)) + 
    geom_tile() + 
    scale_fill_gradient(low = 'transparent', high = '#8B0000', breaks = c(0, 2, 4)) + 
    labs(fill = '') + 
    theme(legend.title = element_text(size = 6, color = "black"),
          legend.text = element_text(size = 9, color = "black")) 
file.up <- '/home/drizzle_zhang/microbiome/result/Figs/ColorBar_up.png'
ggsave(plot = plot_bar_up, filename = file.up,
       height = 5, width = 5, units = 'cm')

plot_bar_down <- 
    ggplot(data = df_plot, aes(x = TF, y = NUL, fill = pval)) + 
    geom_tile() + 
    scale_fill_gradient(low = 'transparent', high = '#00008B', breaks = c(0, 2, 4)) + 
    labs(fill = '') + 
    theme(legend.title = element_text(size = 6, color = "black"),
          legend.text = element_text(size = 9, color = "black")) 
file.down <- '/home/drizzle_zhang/microbiome/result/Figs/ColorBar_down.png'
ggsave(plot = plot_bar_down, filename = file.down,
       height = 5, width = 5, units = 'cm')
