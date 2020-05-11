#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: organize_dataset.py
# @time: 2020/3/6 21:17

from time import time
import pandas as pd
import numpy as np
import os


def build_dataset(dict_in):
    term = dict_in['Biosample term name']
    str_term = term.replace(' ', '_').replace('/', '+').replace("'", '--')
    path_term = os.path.join(path_feature_cell, str_term)

    file_cre = os.path.join(path_cre, f"{str_term}/cRE.txt")
    file_pre = os.path.join(path_term, "correlation.txt")
    df_pre = pd.read_csv(file_pre, sep='\t')

    file_tmp = os.path.join(path_term, "input_file.tmp")
    df_tmp = df_pre.iloc[:, ['gene', 'dhs_id']]
    df_tmp = pd.merge(df_tmp, df_promoter, left_on='gene', right_on=7,
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

    file_ctcf_tmp = os.path.join(path_term, "CTCF.tmp")
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
    df_uniq.columns = \
        ['gene', 'dhs_id', 'score_dhs_insulator', 'score_ctcf_insulator']
    df_genome_ctcf = pd.merge(df_pre, df_uniq, on=['gene', 'dhs_id'],
                              how='outer')
    df_genome_ctcf = df_genome_ctcf.fillna(-10)

    os.remove(file_tmp)
    os.remove(file_ctcf_tmp)

    file_out = os.path.join(path_term, "input_file.txt")
    df_genome_ctcf.to_csv(file_out, sep='\t', index=None, na_rep='NA')

    return


if __name__ == '__main__':
    time_start = time()
    path_cre = '/local/zy/PEI/mid_data/cell_line/DHS/cRE_annotation/'
    path_feature_cell = '/local/zy/PEI/mid_data/cell_line/model_input'
    file_promoter = '/local/zy/PEI/origin_data/gene/' \
                    'promoters.up2k.protein.gencode.v19.unique.bed'
    df_promoter = pd.read_csv(file_promoter, sep='\t', header=None)

    file_meta = '/local/zy/PEI/origin_data/meta_file/meta_model_input.txt'
    df_meta = pd.read_csv(file_meta, sep='\t')
    sub_meta = df_meta.to_dict('records')[0]
    build_dataset(sub_meta)

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
