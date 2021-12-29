#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: organize_dataset.py
# @time: 2020/3/6 21:17

from time import time
import pandas as pd
import numpy as np
import os


def build_dataset(path_cre, path_data, path_inter, path_out, dict_in):
    term = dict_in['Biosample term name']
    str_term = term.replace(' ', '_').replace('/', '+').replace("'", '--')
    pchic = dict_in['pcHi-C_ng2019']
    hic_3div = dict_in['3DIV']
    thurman = dict_in['Thurman']
    col_expressin = dict_in['Gene expression']

    file_positive = os.path.join(
        path_data, f"{str_term}/positive/cRE_pairs.txt")
    file_negative = os.path.join(
        path_data, f"{str_term}/negative/cRE_pairs.sample.1.txt")
    df_positive = pd.read_csv(file_positive, sep='\t', header=None)
    df_positive[4] = np.ones(df_positive.shape[0], dtype='int')
    df_negative = pd.read_csv(file_negative, sep='\t', header=None)
    df_negative[4] = np.zeros(df_negative.shape[0], dtype='int')
    df_data = pd.concat([df_positive, df_negative])
    df_data = df_data.loc[df_data[2].apply(
        lambda x: x in {'Protein-Promoter(Enhancer)',
                        'Other-Promoter(Enhancer)', 'Enhancer'}), :]

    genes = df_data[0].unique()
    sub_promoter = df_promoter.loc[genes, [0, 1, 2, 6]]
    file_pro = os.path.join(path_out, f"{str_term}/promoter.txt")
    file_pro_sort = os.path.join(path_out, f"{str_term}/promoter.sort.txt")
    sub_promoter.to_csv(file_pro, sep='\t', header=False, index=False,
                        encoding='utf-8')
    os.system(f"bedtools sort -i {file_pro} > {file_pro_sort}")

    file_cre = os.path.join(path_cre, f"{str_term}/cRE.txt")
    file_pro_origin = os.path.join(path_out, f"{str_term}/promoter.origin.txt")
    os.system(f"bedtools intersect -a {file_pro_sort} -b {file_cre} -wao | "
              f"cut -f 4,8,9,10,11,12,13,14,15,16 | "
              f"grep -w 'Promoter' | uniq > {file_pro_origin}")
    os.remove(file_pro)
    os.remove(file_pro_sort)
    df_origin = pd.read_csv(file_pro_origin, sep='\t', header=None)

    def _integrate_promoter(df_in):
        df_num = df_in.loc[:, [3, 5, 6, 7, 8]]

        # weighted mean
        # df_out = np.sum(df_num.mul(df_in[9], axis=0)) / np.sum(df_in[9])
        # df_out[0] = df_in[0].unique()[0]

        # max
        df_out = np.max(df_num, axis=0)

        return df_out

    df_integrate = df_origin.groupby(0).apply(_integrate_promoter)
    df_epigenome = pd.merge(df_data, df_integrate, on=0, how='left')
    sub_mat_exp = mat_exp.loc[:, col_expressin]
    df_genome = pd.merge(df_epigenome, sub_mat_exp, left_on=0,
                         right_on='Description', how='left')
    df_positive = df_genome.loc[df_genome[4] == 1, :]
    df_negative = df_genome.loc[df_genome[4] == 0, :]
    df_positive = df_positive.dropna()
    df_negative = df_negative.dropna()
    df_positive['pair'] = df_positive.apply(
        lambda x: x[0] + '_' + x[1], axis=1)
    df_negative['pair'] = df_negative.apply(
        lambda x: x[0] + '_' + x[1], axis=1)
    df_negative = df_negative.loc[
                  df_negative['pair'].apply(
                      lambda x: x not in set(df_positive['pair'].tolist())), :]
    df_negative = df_negative.sample(df_positive.shape[0], random_state=2020)
    genes_negative = set(df_negative[0].tolist())
    df_positive = df_positive.loc[
                  df_positive[0].apply(lambda x: x in genes_negative), :]
    df_negative = df_negative.sample(df_positive.shape[0], random_state=2021)

    # for col in ['5_y', '7_y']:
    #     col_neg = df_negative.loc[:, col].fillna(-10000)
    #     df_negative = df_negative.drop(col, axis=1)
    #     df_negative[col] = col_neg
    # for col in ['6_y', '8_y']:
    #     col_neg = df_negative.loc[:, col].fillna(0)
    #     df_negative = df_negative.drop(col, axis=1)
    #     df_negative[col] = col_neg
    # df_negative = df_negative.sample(df_positive.shape[0], random_state=2020)
    df_genome = pd.concat([df_positive, df_negative], sort=False)
    df_genome = df_genome.drop('pair', axis=1)
    file_tmp = os.path.join(path_out, f"{str_term}/input_file.tmp")
    df_tmp = df_genome.iloc[:, :2]
    df_tmp.columns = ['gene', 'dhs_id']
    df_promoter_tmp = df_promoter
    df_promoter_tmp.index = list(range(df_promoter_tmp.shape[0]))
    df_tmp = pd.merge(df_tmp, df_promoter_tmp, left_on='gene', right_on=6,
                      how='inner')
    df_tmp = df_tmp.loc[:, ['gene', 'dhs_id', 0, 1, 2]]
    df_tmp[3] = df_tmp['dhs_id'].apply(
        lambda x: int(x.split(':')[-1].split('-')[0]))
    df_tmp[4] = df_tmp['dhs_id'].apply(
        lambda x: int(x.split(':')[-1].split('-')[1]))
    df_tmp[5] = df_tmp.apply(
        lambda x: x[2] if x[3] > x[2] else x[4], axis=1)
    df_tmp[6] = df_tmp.apply(
        lambda x: x[3] if x[3] > x[2] else x[1], axis=1)
    df_tmp_out = df_tmp.loc[:, [0, 5, 6, 'gene', 'dhs_id']]
    df_tmp_out.to_csv(file_tmp, sep='\t', header=None, index=None)

    file_ctcf_tmp = os.path.join(path_out, f"{str_term}/CTCF.tmp")
    os.system(f"bedtools intersect -a {file_tmp} -b {file_cre} -wao | "
              f"cut -f 4,5,9,10,11,17 | "
              f"grep -w 'Insulator' | uniq > {file_ctcf_tmp}")
    df_ctcf = pd.read_csv(file_ctcf_tmp, sep='\t', header=None)

    def unique_ctcf(df_in):
        max_ctcf = np.max(df_in[5])
        df_out = df_in.loc[df_in[5] == max_ctcf, [0, 1, 4, 5]]

        return df_out

    df_uniq = df_ctcf.groupby([0, 1]).apply(unique_ctcf)
    df_uniq.index = list(range(df_uniq.shape[0]))
    df_genome_ctcf = pd.merge(df_genome, df_uniq, on=[0, 1], how='outer')
    df_genome_ctcf = df_genome_ctcf.fillna(-10)
    os.remove(file_ctcf_tmp)

    file_pchic = os.path.join(
        path_inter,
        f"{str_term}/{'/'.join(pchic.split('/')[0:2])}/pairs.gene.cRE.txt")
    df_pchic = pd.read_csv(file_pchic, sep='\t', header=None)
    df_pchic = df_pchic.drop(2, axis=1)
    file_3div = os.path.join(
        path_inter,
        f"{str_term}/{'/'.join(hic_3div.split('/')[0:2])}/pairs.gene.cRE.txt")
    df_3div = pd.read_csv(file_3div, sep='\t', header=None)
    df_3div = df_3div.drop(2, axis=1)
    file_thurman = os.path.join(
        path_inter,
        f"{str_term}/{'/'.join(thurman.split('/')[0:2])}/pairs.gene.cRE.txt")
    df_thurman = pd.read_csv(file_thurman, sep='\t', header=None)
    df_thurman = df_thurman.drop(2, axis=1)

    df_merge = pd.merge(df_genome_ctcf, df_pchic, on=[0, 1], how='left')
    df_merge = pd.merge(df_merge, df_3div, on=[0, 1], how='left')
    df_merge = pd.merge(df_merge, df_thurman, on=[0, 1], how='left')

    file_out = os.path.join(path_out, f"{str_term}/input_file.txt")
    df_merge = df_merge.fillna(0)
    df_merge = df_merge.drop_duplicates()

    df_merge.columns = \
        ['gene', 'dhs_id', 'type_enhancer', 'score_dhs_enhancer', 'label',
         'score_h3k4me3_enhancer', 'pval_h3k4me3_enhancer',
         'score_h3k27ac_enhancer', 'pval_h3k27ac_enhancer',
         'distance', 'score_dhs_promoter',
         'score_h3k4me3_promoter', 'pval_h3k4me3_promoter',
         'score_h3k27ac_promoter', 'pval_h3k27ac_promoter',
         'gene_expression', 'score_dhs_insulator', 'score_ctcf_insulator',
         'pcHi-C_ng2019', '3DIV', 'Thurman']

    # df_merge_neg = df_merge.loc[df_merge['label'] == 0, :]
    # df_merge_pos = df_merge.loc[df_merge['label'] == 1, :]
    # df_merge_neg_0 = df_merge_neg.loc[
    #                  (df_merge_neg['pcHi-C_ng2019'] == 0) &
    #                  (df_merge_neg['Thurman'] == 0) &
    #                  (df_merge_neg['3DIV'] == 0), :]
    # df_merge_neg_pc1 = df_merge_neg.loc[df_merge['pcHi-C_ng2019'] == 1, :]
    # df_merge_neg_t1 = df_merge_neg.loc[df_merge['Thurman'] == 1, :]
    # df_merge_neg_0 = df_merge_neg_0.sample(
    #     int(0.97*df_positive.shape[0]),
    #     random_state=2020)
    # df_merge_neg_pc1 = df_merge_neg_pc1.sample(
    #     int(0.02*df_positive.shape[0]), random_state=2021)
    # df_merge_neg_t1 = df_merge_neg_t1.sample(
    #     int(0.01*df_positive.shape[0]), random_state=2022)
    # df_merge = pd.concat(
    #     [df_merge_pos, df_merge_neg_0, df_merge_neg_pc1, df_merge_neg_t1],
    #     sort=False)

    df_merge.to_csv(file_out, sep='\t', index=None, na_rep='NA')

    return


if __name__ == '__main__':
    time_start = time()
    path_ref = '/local/zy/PEI/mid_data/cell_line/DHS/cRE_annotation/'
    path_dataset = '/local/zy/PEI/mid_data/cell_line/eQTL/'
    path_interaction = '/local/zy/PEI/mid_data/cell_line/gene_cre_pairs'
    file_meta = '/local/zy/PEI/origin_data/meta_file/meta_chrom_inter.txt'
    path_data_out = '/local/zy/PEI/mid_data/cell_line/model_input'
    list_inter = ['pcHi-C_ng2019', '3DIV', 'Thurman']
    file_promoter = '/local/zy/PEI/origin_data/gene/' \
                    'promoters.up2k.protein.gencode.v19.bed'
    df_promoter = pd.read_csv(file_promoter, sep='\t', header=None)
    df_promoter[6] = df_promoter[3].apply(lambda x: x.split('<-')[0])
    df_promoter = df_promoter.drop_duplicates(subset=6)
    df_promoter.index = df_promoter[6]

    file_expression = \
        '/local/zy/PEI/origin_data/GTEx/' \
        'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct'
    mat_exp = pd.read_csv(file_expression, sep='\t', index_col=1)
    mat_exp = mat_exp.drop('gene_id', axis=1)

    df_meta = pd.read_csv(file_meta, sep='\t')
    sub_meta = df_meta.to_dict('records')[0]
    build_dataset(path_ref, path_dataset, path_interaction,
                  path_data_out, sub_meta)

    time_end = time()
    print(time_end - time_start)

    # try
    file_input = \
        '/local/zy/PEI/mid_data/cell_line/model_input/GM12878/input_file.txt'
    df_merge = pd.read_csv(file_input, sep='\t')
    # df_merge.columns = \
    #     ['gene', 'dhs_id', 'type_enhancer', 'score_dhs_enhancer', 'label',
    #      'score_h3k4me3_enhancer', 'pval_h3k4me3_enhancer',
    #      'score_h3k27ac_enhancer', 'pval_h3k27ac_enhancer',
    #      'distance', 'score_dhs_promoter',
    #      'score_h3k4me3_promoter', 'pval_h3k4me3_promoter',
    #      'score_h3k27ac_promoter', 'pval_h3k27ac_promoter',
    #      'gene_expression', 'pcHi-C_ng2019', '3DIV', 'Thurman']
    # df_merge_pos = df_merge.loc[df_merge['label'] == 1, :]
    # df_merge_neg = df_merge.loc[df_merge['label'] == 0, :]
    # df_merge_pos_drop = df_merge_pos.loc[
    #                     df_merge_pos['pval_h3k4me3_promoter'] != 0, :]
    # df_merge_neg_drop = df_merge_neg.loc[
    #                     df_merge_neg['pval_h3k4me3_promoter'] != 0, :]
    # df_merge_pos_0 = df_merge_pos_drop.loc[
    #                  (df_merge_pos_drop['pcHi-C_ng2019'] == 0) &
    #                  (df_merge_pos_drop['Thurman'] == 0) &
    #                  (df_merge_pos_drop['3DIV'] == 0), :]
    # df_merge_neg_0 = df_merge_neg_drop.loc[
    #                  (df_merge_neg_drop['pcHi-C_ng2019'] == 0) &
    #                  (df_merge_neg_drop['Thurman'] == 0) &
    #                  (df_merge_neg_drop['3DIV'] == 0), :]
    # df_merge_pos_pc1 = df_merge_pos.loc[df_merge_pos['pcHi-C_ng2019'] == 1, :]
    # df_merge_neg_pc1 = df_merge_neg.loc[df_merge_neg['pcHi-C_ng2019'] == 1, :]
    # df_merge_pos_t1 = df_merge_pos.loc[df_merge_pos['Thurman'] == 1, :]
    # df_merge_neg_t1 = df_merge_neg.loc[df_merge_neg['Thurman'] == 1, :]
    # df_gene_pos = df_merge_pos_drop['gene'].unique()
    # df_gene_neg = df_merge_neg_drop['gene'].unique()
    #
    df_merge_drop = df_merge.loc[
                    df_merge['pval_h3k4me3_promoter'] != 0, :]
    df_corr = df_merge_drop.corr()
    file_corr = '/local/zy/PEI/mid_data/cell_line/model_input/GM12878/corr.txt'
    df_corr.to_csv(file_corr, sep='\t')
