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

#############################################
library(ggplot2)
vec.dataset <- c('Habib', 'Tasic', 'Zeisel')
setwd('/home/drizzle_zhang/scRef/try_data/evaluation_del_cell')
for (dataset in vec.dataset) {
    rds <- readRDS(paste0('plot_', dataset, '_scRef2.Rdata'))
    df.plot <- rds$plot
    # modify plot dataframe
    df.uplimit <- df.plot[df.plot$V3 == 'UpLimit',]
    merge.uplimit <- df.uplimit[df.uplimit$V2 == 'Neuron merged', 'V4']
    remove.uplimit <- df.uplimit[df.uplimit$V2 == 'Neuron removed', 'V4']
    df.plot.modify <- 
        df.plot[df.plot$V3 != 'UpLimit' & df.plot$V7 == 'Macro F1',]
    for (cell in unique(df.plot.modify$V3)) {
        df.plot.modify <- rbind(df.plot.modify,
                                data.frame(V1 = dataset, V2 = 'Neuron merged',
                                           V3 = cell, V4 = merge.uplimit,
                                           V5 = 'Up Limit', V6 = '-1', 
                                           V7 = 'Macro F1', V8 = 'scRef'))
        df.plot.modify <- rbind(df.plot.modify,
                                data.frame(V1 = dataset, V2 = 'Neuron removed',
                                           V3 = cell, V4 = remove.uplimit,
                                           V5 = 'Up Limit', V6 = '-1', 
                                           V7 = 'Macro F1', V8 = 'scRef'))
    }
    for (i in row.names(df.plot.modify)) {
        if (df.plot.modify[i, 'V5'] == 'control') {
            df.plot.modify[i, 'V9'] <- 'scRef'
        }
        if (df.plot.modify[i, 'V5'] == 'Up Limit') {
            df.plot.modify[i, 'V9'] <- 'Upper limit'
        }
        if (df.plot.modify[i, 'V5'] == 'best cutoff') {
            df.plot.modify[i, 'V9'] <- 'scRef plus(best cutoff)'
        }
        if (df.plot.modify[i, 'V5'] == 'default cutoff') {
            df.plot.modify[i, 'V9'] <- 'scRef plus(default cutoff)'
        }
    }
    
    plot.f1 <- ggplot(df.plot.modify,
                      aes(x = V3, y = V4, fill = V9)) +
        geom_bar(position = 'dodge', stat = 'identity') +
        facet_wrap(~ V2, scales = 'free', ncol = 2) +
        labs(title = dataset, y = 'Weighted Macro F1', 
             x = 'Deleted Cell', fill = 'Method') + 
        coord_cartesian(ylim = c(0.7, 1)) + 
        theme(axis.text = element_text(size = 10),
              axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(filename = paste0('./performance_', dataset, '.png'),
           plot = plot.f1)
    
    df.correct <- df.plot[df.plot$V7 == 'Rate of error correction',]
    plot.correct <- ggplot(df.correct,
                           aes(x = V3, y = V4, fill = V5)) +
        geom_bar(position = 'dodge', stat = 'identity') +
        facet_wrap(~ V2, scales = 'free', ncol = 2) +
        labs(title = dataset, y = 'F1 score', 
             x = 'Deleted Cell', fill = 'Method') + 
        coord_cartesian(ylim = c(0.7, 1)) + 
        theme(axis.text = element_text(size = 10),
              axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(filename = paste0('./correction_', dataset, '.png'),
           plot = plot.correct)
    
    # change cutoff
    list.metrics <- rds$metrics
    df.metrics <- data.frame(stringsAsFactors = F)
    for (cell in unique(df.correct$V3)) {
        sub.list <- list.metrics[[cell]]
        df.merge <- sub.list[['merge']]
        df.merge$cell <- rep(cell, dim(df.merge)[1])
        df.merge$type <- rep('Neuron merged', dim(df.merge)[1])
        df.remove <- sub.list[['remove']]
        df.remove$cell <- rep(cell, dim(df.remove)[1])
        df.remove$type <- rep('Neuron removed', dim(df.remove)[1])
        df.metrics <- rbind(df.metrics, df.merge)
        df.metrics <- rbind(df.metrics, df.remove)
    }
    
    plot.metrics <- ggplot(df.metrics,
                           aes(x = cutoff, y = weighted.f1, color = cell)) +
        geom_line() +
        facet_wrap(~ type, scales = 'free', ncol = 2) +
        labs(title = dataset, y = 'Weighted Macro F1', 
             x = 'Cutoff', color = 'Deleted cell') + 
        coord_cartesian(ylim = c(0.6, 1))
    ggsave(filename = paste0('./metrics_', dataset, '.png'),
           plot = plot.metrics)
    
}

