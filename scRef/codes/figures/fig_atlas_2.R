library(ggplot2)

path.res <- '/mdshare/node9/zy/scRef/figure/atlas_anno/'
vec.ref <- rep(c('MCA', 'MCA', 'HCL', 'HCL'), 2)
vec.dataset <- rep(c('Tasic2018', 'Haber_Duodenum', 'panc8_indrop', 'hPBMC_10Xv2'), 2)
vec.title <- rep(c('Mouse Neocortex', 'Mouse Duodenum', 
                   'Human Pancreas', 'Human PBMC'), 2)
vec.methods <- c(rep('scMAGIC', 4), rep('scHCL', 4))

df.plot <- data.frame(stringsAsFactors = F)
for (i in 1:length(vec.dataset)) {
    reference <- vec.ref[i]
    dataset <- vec.dataset[i]
    title <- vec.title[i]
    method <- vec.methods[i]
    file.res.scMAGIC <- paste0(path.res, 'RES_', reference,'_', dataset, '_', method, '.Rdata')
    res.scMAGIC <- readRDS(file.res.scMAGIC)
    if (reference == 'MCA' & method == 'scHCL') {method <- 'scMCA'}
    df.sub <- data.frame(term = 'Accuracy', method = method,
                         value = res.scMAGIC$accuracy, stringsAsFactors = F)
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Balanced accuracy', method = method,
                               value = res.scMAGIC$balanced.accuracy, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'labeled Accuracy', method = method,
                               value = res.scMAGIC$accuracy.rm.unassigned, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'labeled Balanced accuracy', method = method,
                               value = res.scMAGIC$balanced.accuracy.rm.unassigned, stringsAsFactors = F))
    df.sub$dataset <- rep(title, nrow(df.sub))
    df.plot <- rbind(df.plot, df.sub)
}
df.plot$dataset <- factor(df.plot$dataset, levels = c('Mouse Neocortex', 'Mouse Duodenum', 
                                                      'Human Pancreas', 'Human PBMC'))
df.plot$term <- factor(df.plot$term, 
                       levels = rev(c('labeled Balanced accuracy', 'labeled Accuracy',
                                  'Balanced accuracy', 'Accuracy')))
df.plot$method <- factor(df.plot$method, levels = rev(c('scMCA', 'scHCL', 'scMAGIC')))
# library(RColorBrewer)
# brewer.pal(11,"Paired")
plot.bar <-
    ggplot(df.plot, aes(x = method, y = value, color = term, fill = term)) +
    geom_bar(stat = 'identity', position=position_dodge(0.7), width = 0.7) +
    facet_grid( ~ dataset, scales = 'free', space = 'free') + 
    scale_color_manual(breaks = c('Accuracy', 'Balanced accuracy', 'labeled Accuracy',
                                  'labeled Balanced accuracy'),
                       values = c("#80B1D3", "#5AB4AC", "#BC80BD", "#FFA07A")) +
    scale_fill_manual(breaks = c('Accuracy', 'Balanced accuracy', 'labeled Accuracy',
                                 'labeled Balanced accuracy'),
                      values = c("#80B1D3", "#5AB4AC", "#BC80BD", "#FFA07A")) +
    labs(title = "", y = '', x = '', color = '', fill = '') +
    theme_bw() +
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          strip.text.x = element_text(
              size = 10, color = "black", family = 'Arial', angle = 0), 
          strip.background = element_rect(
              color = "transparent", fill = "transparent"),
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 10, color = 'black', family = 'Arial'),
          axis.text.y = element_text(size = 8, color = 'black', family = 'Arial'),
          axis.title = element_text(size = 10, family = 'Arial'),
          legend.text = element_text(size = 10, family = 'Arial'),
          legend.position = 'right',
          # legend.key.size = unit(1, 'cm'),
          legend.key = element_blank())
ggsave(filename = 'Balanced_accuracy_atlas_methods.png',
       path = path.res, plot = plot.bar,
       units = 'cm', height = 10, width = 21)



