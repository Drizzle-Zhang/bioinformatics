library(ggplot2)
setwd('/home/drizzle_zhang/driver_mutation/cRE_plot/model_test')

df_feature <- data.frame(
    Feature = factor(c('score_dhs_enhancer', 'score_h3k4me3_enhancer', 'pval_h3k4me3_enhancer', 
                'score_h3k27ac_enhancer', 'pval_h3k27ac_enhancer', 'distance', 
                'score_dhs_promoter', 'score_h3k4me3_promoter', 'pval_h3k4me3_promoter', 
                'score_h3k27ac_promoter', 'pval_h3k27ac_promoter', 'gene_expression', 
                'pcHi-C_ng2019', '3DIV', 'Thurman'), 
                levels = c('score_dhs_enhancer', 'score_h3k4me3_enhancer', 'pval_h3k4me3_enhancer', 
                           'score_h3k27ac_enhancer', 'pval_h3k27ac_enhancer', 'distance', 
                           'score_dhs_promoter', 'score_h3k4me3_promoter', 'pval_h3k4me3_promoter', 
                           'score_h3k27ac_promoter', 'pval_h3k27ac_promoter', 'gene_expression', 
                           'pcHi-C_ng2019', '3DIV', 'Thurman')),
    Importance = c(0.06509043, 0.06040399, 0.05754984, 0.06348088, 0.06302784, 
                   0.12191288, 0.10307148, 0.09548811, 0.09602043, 0.08009987, 
                   0.07834308, 0.10562864, 0.00443098, 0.00075129, 0.00470028))

ggplot(df_feature, aes(x = Feature, y = Importance)) + 
    geom_bar(stat = "identity") + coord_flip() + 
    theme(axis.ticks = element_blank())+
    theme(axis.title = element_text(size = 17)) +   
    theme(axis.text = element_text(size = 12)) +   
    theme(panel.background=element_blank())   ## ȥ?����???ɫ


df_feature_corr <- data.frame(
    Feature = factor(c('score_dhs_enhancer', 'score_h3k4me3_enhancer', 'pval_h3k4me3_enhancer', 
                       'score_h3k27ac_enhancer', 'pval_h3k27ac_enhancer', 'distance', 
                       'score_dhs_promoter', 'score_h3k4me3_promoter', 'pval_h3k4me3_promoter', 
                       'score_h3k27ac_promoter', 'pval_h3k27ac_promoter', 'gene_expression', 
                       'pcHi-C_ng2019', '3DIV', 'Thurman'), 
                     levels = c('score_dhs_enhancer', 'score_h3k4me3_enhancer', 'pval_h3k4me3_enhancer', 
                                'score_h3k27ac_enhancer', 'pval_h3k27ac_enhancer', 'distance', 
                                'score_dhs_promoter', 'score_h3k4me3_promoter', 'pval_h3k4me3_promoter', 
                                'score_h3k27ac_promoter', 'pval_h3k27ac_promoter', 'gene_expression', 
                                'pcHi-C_ng2019', '3DIV', 'Thurman')),
    Importance = c(-0.008468, 0.019101, 0.057032, 0.046944, 0.018738, 
                   -0.658697, -0.086066, -0.023175, -0.039543, -0.030743, 
                   0.006081, -0.015617, 0.073823, -0.021836, 0.135527))

ggplot(df_feature_corr, aes(x = Feature, y = Importance)) + 
    geom_bar(stat = "identity") + coord_flip() + 
    theme(axis.ticks = element_blank())+
    theme(axis.title = element_text(size = 17)) +   
    theme(axis.text = element_text(size = 12)) +   
    theme(panel.background=element_blank())   ## ȥ?����???ɫ


df_feature_corr <- data.frame(
    Feature = factor(c('score_dhs_enhancer', 
                'score_h3k4me3_enhancer', 'pval_h3k4me3_enhancer',
                'score_h3k27ac_enhancer', 'pval_h3k27ac_enhancer',
                'distance', 'score_dhs_promoter',
                'score_h3k4me3_promoter', 'pval_h3k4me3_promoter',
                'score_h3k27ac_promoter', 'pval_h3k27ac_promoter',
                'gene_expression', 'score_dhs_insulator', 'score_ctcf_insulator',
                'pcHi-C_ng2019', '3DIV', 'Thurman'), 
                levels = c('score_dhs_enhancer', 
                           'score_h3k4me3_enhancer', 'pval_h3k4me3_enhancer',
                           'score_h3k27ac_enhancer', 'pval_h3k27ac_enhancer',
                           'gene_expression', 'score_dhs_promoter',
                           'score_h3k4me3_promoter', 'pval_h3k4me3_promoter',
                           'score_h3k27ac_promoter', 'pval_h3k27ac_promoter',
                           'distance', 
                           'score_dhs_insulator', 'score_ctcf_insulator',
                           'pcHi-C_ng2019', '3DIV', 'Thurman')),
    Importance = c(0.04808697, 0.04425646, 0.04335354, 0.05527308, 0.04855947, 
                   0.10078233, 0.0793875, 0.0761882, 0.07374874, 0.06673095, 
                   0.05808765, 0.08013773, 0.07231681, 0.07006732, 0.06248565, 
                   0.00378482, 0.01675277))

ggplot(df_feature_corr, aes(x = Feature, y = Importance)) + 
    geom_bar(stat = "identity") + coord_flip() + 
    theme(axis.ticks = element_blank())+
    theme(axis.title = element_text(size = 17)) +   
    theme(axis.text = element_text(size = 12)) +   
    theme(panel.background=element_blank())   ## ȥ?����???ɫ



# compare among different models
df.accuracy <- data.frame(
    Model = c(rep('Model(distance matched)', 5),
              rep('Model(without distance matched)', 5)),
    Evaluator = c('F1 Score', 'Precision', 'Recall', 'AUC of ROC', 'AUC of PRC',
                  'F1 Score', 'Precision', 'Recall', 'AUC of ROC', 'AUC of PRC'),
    Value = c(0.72543, 0.70133, 0.75125, 0.79244, 0.78579,
              0.87345, 0.89694, 0.85120, 0.94325, 0.94112)
)

df.accuracy <- data.frame(
    Model = c(rep('Model(distance matched)', 5),
              rep('Model(without distance matched, no CTCF)', 5),
              rep('Model(without distance matched)', 5)),
    Evaluator = c('F1 Score', 'Precision', 'Recall', 'AUC of ROC', 'AUC of PRC',
                  'F1 Score', 'Precision', 'Recall', 'AUC of ROC', 'AUC of PRC',
                  'F1 Score', 'Precision', 'Recall', 'AUC of ROC', 'AUC of PRC'),
    Value = c(0.72543, 0.70133, 0.75125, 0.79244, 0.78579,
              0.64133, 0.64327, 0.63947, 0.69948, 0.69776,
              0.87345, 0.89694, 0.85120, 0.94325, 0.94112)
)

ggplot(df.accuracy,
       aes(x = Model, y = Value, fill = Evaluator)) +
    geom_bar(position = 'dodge', stat = 'identity') +
    labs(title = "", y = '', x = '') + 
    theme(axis.text = element_text(size = 12)) +   
    coord_cartesian(ylim = c(0.60, 1))
# scale_fill_discrete(
#     name = 'Evaluation methods',
#     breaks = c('pcReg', 'sil', 'kBET', 'mixent',
#                'ARI', 'NMI', 'ldaReg', 'sd'),
#     labels = c(
#         'PC Regression',
#         'Silhouettes',
#         'kBET',
#         'Entropy of batch mixing',
#         'Adjusted rand index',
#         'Mormalized mutual information',
#         'LDA Regression',
#         'Standard Deviation'
#     )
# )


df_feature_corr <- data.frame(
    Feature = factor(
        c('DHS_DHS', 'H3K4me3_DHS', 'DHS_H3K27ac',
          'H3K4me3_H3K27ac', 'score_ctcf_insulator'), 
        levels = c('DHS_DHS', 'H3K4me3_DHS', 'DHS_H3K27ac',
                   'H3K4me3_H3K27ac', 'score_ctcf_insulator')),
        Importance = c(0.19175378, 0.18946368, 0.21554017, 0.20636483, 
                       0.19687755))

ggplot(df_feature_corr, aes(x = Feature, y = Importance)) + 
    geom_bar(stat = "identity") + coord_flip() + 
    theme(axis.ticks = element_blank())+
    theme(axis.title = element_text(size = 17)) +   
    theme(axis.text = element_text(size = 12)) +   
    theme(panel.background=element_blank())   


df_feature_corr <- data.frame(
    Feature = factor(
        c('DHS_DHS', 'H3K4me3_DHS', 'DHS_H3K27ac',
          'H3K4me3_H3K27ac', 'score_ctcf_insulator'), 
        levels = c('DHS_DHS', 'H3K4me3_DHS', 'DHS_H3K27ac',
                   'H3K4me3_H3K27ac', 'score_ctcf_insulator')),
    Importance = c(0.19019998, 0.18159677, 0.20561532, 0.19609595, 
                   0.22649198))

ggplot(df_feature_corr, aes(x = Feature, y = Importance)) + 
    geom_bar(stat = "identity") + coord_flip() + 
    theme(axis.ticks = element_blank())+
    theme(axis.title = element_text(size = 17)) +   
    theme(axis.text = element_text(size = 12)) +   
    theme(panel.background=element_blank())   


df_feature_corr <- data.frame(
    Dataset = factor(
        c('Cross_validation', 'H1', 'IMR-90'), 
        levels = c('Cross_validation', 'H1', 'IMR-90')),
    AUROC = c(0.9422056376836075, 0.8452468496196388, 0.8492239553508365))

ggplot(df_feature_corr, aes(x = Dataset, y = AUROC)) + 
    geom_bar(stat = "identity") + coord_flip() + 
    theme(axis.ticks = element_blank()) +
    theme(axis.title = element_text(size = 17)) +   
    theme(axis.text = element_text(size = 12)) +   
    coord_cartesian(ylim = c(0.50, 1)) +   
    theme(panel.background = element_blank()) 

################################# 
df.accuracy <- data.frame(
    TestMethod = c(rep('Random split', 6),
                   rep('Chrom split', 6)),
    Model = c(rep('Full model', 2),
              rep('Full model - CTCF', 2),
              rep('Only CTCF', 2),
              rep('Full model', 2),
              rep('Full model - CTCF', 2),
              rep('Only CTCF', 2)),
    Evaluator = c('AUC of PRC', 'AUC of ROC',
                  'AUC of PRC', 'AUC of ROC',
                  'AUC of PRC', 'AUC of ROC',
                  'AUC of PRC', 'AUC of ROC',
                  'AUC of PRC', 'AUC of ROC',
                  'AUC of PRC', 'AUC of ROC'),
    Value = c(0.55136, 0.95148, 0.23839, 0.78941, 0.50187,
              0.92050, 0.45323, 0.92791, 0.13468, 0.72471,
              0.40021, 0.88268)
)

ggplot(df.accuracy,
       aes(x = Model, y = Value, fill = TestMethod)) +
    geom_bar(position = 'dodge', stat = 'identity') +
    facet_wrap(~ Evaluator, scales = 'free_x', ncol = 2) +
    labs(title = "", y = '', x = '') + 
    theme(axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1))   
    
    # coord_cartesian(ylim = c(0, 1))






