# scatter
# GM12878 DHS_DHS
file.in <- '/lustre/tianlab/zhangyu/PEI/mid_data_correct/cell_line/model_input/GM12878/GM12878_feature_label.txt'
df.feature.label <- read.delim(file.in, stringsAsFactors = F)
df.feature.label.plot <- 
    df.feature.label[sample.int(dim(df.feature.label)[1], round(dim(df.feature.label)[1]/10)),]

fig.out <- '/lustre/tianlab/zhangyu/PEI/mid_data_correct/cell_line/model_input/GM12878/GM12878_scatter_DHS_DHS.pdf'
ggplot.scatter <- ggplot(
    df.feature.label.plot, aes(x = signal_DHS_DHS, y = DHS_DHS, color = label)) + 
    geom_point(size = 0.1)
ggsave(fig.out, ggplot.scatter)

cutoff.signal_DHS_DHS <- quantile(df.feature.label$signal_DHS_DHS, 0.8)
cutoff.DHS_DHS <- quantile(df.feature.label$DHS_DHS, 0.8)
