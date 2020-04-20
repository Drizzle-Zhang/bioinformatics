library(ggplot2)

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
    theme(panel.background=element_blank())   ## È¥µô±³¾°ÑÕÉ«


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
    theme(panel.background=element_blank())   ## È¥µô±³¾°ÑÕÉ«


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
    theme(panel.background=element_blank())   ## È¥µô±³¾°ÑÕÉ«



# compare among different models
df.accuracy <- data.frame(
    Model = c(rep('Baseline', 4), rep('Add CTCF and continuous pcHi-C', 4), 
              rep('Simple Model', 4)),
    Evaluator = c('F1 Score', 'Precision', 'Recall', 'AUC of ROC', 
                  'F1 Score', 'Precision', 'Recall', 'AUC of ROC', 
                  'F1 Score', 'Precision', 'Recall', 'AUC of ROC'),
    Value = c(0.92790, 0.90936, 0.94725, 0.97448,
              0.97350, 0.96844, 0.97865, 0.99171,
              0.96392, 0.95434, 0.97376, 0.99035)
)

ggplot(df.accuracy,
       aes(x = Model, y = Value, fill = Evaluator)) +
    geom_bar(position = 'dodge', stat = 'identity') +
    labs(title = "", y = '', x = '') + 
    theme(axis.text = element_text(size = 12)) +   
    coord_cartesian(ylim = c(0.85, 1))
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
        c('score_h3k27ac_enhancer', 'distance', 'score_h3k4me3_promoter',
          'gene_expression', 'score_ctcf_insulator',
          'pcHi-C_ng2019', '3DIV', 'Thurman'), 
        levels = c('score_h3k27ac_enhancer', 'score_h3k4me3_promoter',
                   'gene_expression', 'distance', 'score_ctcf_insulator',
                   'pcHi-C_ng2019', '3DIV', 'Thurman')),
        Importance = c(0.13116104, 0.15870940, 0.19263499, 0.18950501, 
                       0.16377640, 0.12306651, 0.00732349, 0.03382314))

ggplot(df_feature_corr, aes(x = Feature, y = Importance)) + 
    geom_bar(stat = "identity") + coord_flip() + 
    theme(axis.ticks = element_blank())+
    theme(axis.title = element_text(size = 17)) +   
    theme(axis.text = element_text(size = 12)) +   
    theme(panel.background=element_blank())   ## È¥µô±³¾°ÑÕÉ«
