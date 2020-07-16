library(ggplot2)

df.evaluate <- 
    data.frame(test.method = c(rep('Neuron merged', 16), rep('Neuron removed', 16)),
               dataset = c(rep('Zeisel', 4), rep('Tasic', 4),
                           rep('Habib', 4), rep('Hochgerner', 4),
                           rep('Zeisel', 4), rep('Tasic', 4),
                           rep('Habib', 4), rep('Hochgerner', 4)),
               method = c('Control', 'Wilcox test', 'T test', 'Permutation test',
                          'Control', 'Wilcox test', 'T test', 'Permutation test',
                          'Control', 'Wilcox test', 'T test', 'Permutation test',
                          'Control', 'Wilcox test', 'T test', 'Permutation test'),
               macro.f1 = c(0.8924, 0.9479, 0.9492, 0.9492,
                            0.9668, 0.9954, 0.9953, 0.9942,
                            0.9879, 0.9922, 0.9925, 0.9957,
                            0.8318, 0.6700, 0.4581, 0.6209,
                            0.8035, 0.9268, 0.9281, 0.9281, 
                            0.6338, 0.9550, 0.9537, 0.9412,
                            0.9139, 0.9678, 0.9692, 0.9777,
                            0.3621, 0.8883, 0.8047, 0.9001))

ggplot(df.evaluate,
       aes(x = dataset, y = macro.f1, fill = method)) +
    geom_bar(position = 'dodge', stat = 'identity') +
    facet_wrap(~ test.method, scales = 'free', ncol = 2) +
    labs(title = "", y = 'Weighted Macro F1', x = '', fill = 'Method') + 
    coord_cartesian(ylim = c(0.25, 1)) + 
    theme(axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1))   
