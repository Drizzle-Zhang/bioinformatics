#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: organize_dataset.py
# @time: 2020/3/6 21:17

from time import time
import pandas as pd
import numpy as np
import os
from multiprocessing import Pool


def build_dataset(dict_in):
    term = dict_in['Biosample term name']
    str_term = term.replace(' ', '_').replace('/', '+').replace("'", '--')
    path_term = os.path.join(path_feature_cell, str_term)

    file_cre = os.path.join(path_cre_cell, f"{str_term}/cRE.CTCF.txt")
    if not os.path.isfile(file_cre):
        return
    file_pre = os.path.join(path_term, "correlation.txt")
    df_pre = pd.read_csv(file_pre, sep='\t')
    abs_distance = np.abs(df_pre['distance'])
    df_pre = df_pre.loc[abs_distance > 3000, :]

    file_tmp = os.path.join(path_term, "input_file.tmp")
    df_tmp = df_pre.loc[:, ['gene', 'dhs_id']]
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
    # df_tmp.loc[(df_tmp['gene'] == 'TADA2B') &
    #            (df_tmp['dhs_id'] == 'DHS<-chr4:7045626-7045765'), :]
    df_tmp_out = df_tmp.loc[:, [0, 5, 6, 'gene', 'dhs_id']]
    df_tmp_out.to_csv(file_tmp, sep='\t', header=None, index=None)

    file_cre_ctcf = os.path.join(path_term, "cRE.CTCF.tmp")
    os.system(f"grep -w 'Insulator' {file_cre} > {file_cre_ctcf}")
    file_ctcf_tmp = os.path.join(path_term, "CTCF.tmp")
    os.system(f"bedtools intersect -a {file_tmp} -b {file_cre_ctcf} -wao | "
              f"cut -f 4,5,9,10,11,17 | uniq > {file_ctcf_tmp}")
    df_ctcf = pd.read_csv(file_ctcf_tmp, sep='\t', header=None)

    def unique_ctcf(df_in):
        max_ctcf = np.max(df_in[5])
        # df_out = df_in.loc[df_in[5] == max_ctcf, [0, 1, 4, 5]]
        df_out = df_in.loc[df_in[5] == max_ctcf, [0, 1, 5]]

        return df_out

    df_uniq = df_ctcf.groupby([0, 1]).apply(unique_ctcf)
    df_uniq.index = list(range(df_uniq.shape[0]))
    # df_uniq.columns = \
    #     ['gene', 'dhs_id', 'score_dhs_insulator', 'score_ctcf_insulator']
    df_uniq.columns = ['gene', 'dhs_id', 'score_ctcf_insulator']
    df_genome_ctcf = pd.merge(df_pre, df_uniq, on=['gene', 'dhs_id'],
                              how='outer')
    df_genome_ctcf = df_genome_ctcf.fillna(-10)

    os.remove(file_tmp)
    os.remove(file_ctcf_tmp)

    file_out = os.path.join(path_term, "input_file.txt")
    df_genome_ctcf.to_csv(file_out, sep='\t', index=None, na_rep='NA')

    return


def build_training_set():
    file_num_cell = os.path.join(path_label, 'Cell_num_pairs.txt')
    df_num_cell = pd.read_csv(file_num_cell, sep='\t')
    cells = df_num_cell.loc[df_num_cell['Num of pairs'] > 0, 'Cell line']

    file_all = os.path.join(path_label, 'All_interactions.txt')
    df_all = pd.read_csv(file_all, sep='\t')
    set_all = set(df_all.apply(
        lambda x: f"{x['gene']}--{x['ref_dhs_id']}", axis=1))

    list_df_cell = []
    for cell in cells:
        file_label = os.path.join(path_label, f'{cell}.txt')
        file_feature = os.path.join(
            path_feature_cell, f'{cell}/input_file.txt')
        df_feature = pd.read_csv(file_feature, sep='\t')
        df_feature['abs_distance'] = np.abs(df_feature['distance'])
        df_label = pd.read_csv(file_label, sep='\t')
        df_positive = pd.merge(df_feature, df_label,
                               on=['gene', 'dhs_id', 'ref_dhs_id', 'type_cre'])
        assert df_positive.shape[0] == df_label.shape[0]
        df_positive = df_positive.loc[
                      df_positive['abs_distance'] < 2_000_000, :]
        df_positive['label'] = np.full(df_positive.shape[0], 1)

        df_feature.index = df_feature.apply(
            lambda x: f"{x['gene']}--{x['ref_dhs_id']}", axis=1)
        set_cell_all = set(df_feature.index)
        set_cell_all_neg = set_cell_all.difference(set_all)
        df_feature_neg = df_feature.loc[set_cell_all_neg, :]

        len_bin = 20000
        num_pos = df_positive.shape[0]
        num_neg = df_feature_neg.shape[0]
        list_df = []
        for i in range(0, 100):
            down_limit = i * len_bin
            up_limit = (i + 1) * len_bin
            sub_pos = df_positive.loc[
                      (df_positive['abs_distance'] > down_limit) &
                      (df_positive['abs_distance'] <= up_limit), :]
            # num_sample = int(num_neg * (sub_pos.shape[0]/num_pos))
            num_sample = sub_pos.shape[0]
            sub_neg = df_feature_neg.loc[
                      (df_feature_neg['abs_distance'] > down_limit) &
                      (df_feature_neg['abs_distance'] <= up_limit), :]
            # sub_neg_sample = sub_neg.sample(
            #     min(num_sample, sub_neg.shape[0]), random_state=123)
            sub_neg_sample = sub_neg.sample(num_sample, random_state=123)
            list_df.append(sub_neg_sample)

        df_negative = pd.concat(list_df, sort=False)
        df_negative['label'] = np.full(df_negative.shape[0], 0)
        print(cell)
        print(df_positive.shape[0])
        print(df_negative.shape[0])
        df_cell = pd.concat([df_positive, df_negative])
        df_cell = df_cell.drop('abs_distance', axis=1)
        file_cell_out = os.path.join(
            path_feature_cell, f'{cell}/training_set.txt')
        df_cell.to_csv(file_cell_out, sep='\t', index=None)
        df_cell['cell line'] = [cell for _ in range(df_cell.shape[0])]
        list_df_cell.append(df_cell)

    df_training_set = pd.concat(list_df_cell, sort=False)
    file_label = os.path.join(path_label, f'training_set.txt')
    df_training_set.to_csv(file_label, sep='\t', index=None)

    return


if __name__ == '__main__':
    time_start = time()
    path_cre_cell = '/local/zy/PEI/mid_data/cell_line/DHS/cRE_annotation/'
    path_feature_cell = '/local/zy/PEI/mid_data/cell_line/model_input'
    file_promoter = '/local/zy/PEI/origin_data/gene/' \
                    'promoters.up2k.protein.gencode.v19.unique.bed'
    df_promoter = pd.read_csv(file_promoter, sep='\t', header=None)

    df_meta_cell = pd.read_csv(
        os.path.join(path_cre_cell, 'meta.reference.tsv'), sep='\t')
    cell_dicts = df_meta_cell.to_dict('records')
    pool = Pool(40)
    pool.map(build_dataset, cell_dicts)
    pool.close()

    # sub_meta = cell_dicts[4]
    # build_dataset(sub_meta)

    path_label = \
        '/local/zy/PEI/mid_data/training_label/label_interactions_V1'
    build_training_set()

    time_end = time()
    print(time_end - time_start)

    # try
    # file_input = \
    #     '/local/zy/PEI/mid_data/cell_line/model_input/GM12878/input_file.txt'
    # df_merge = pd.read_csv(file_input, sep='\t')
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
    # df_merge_drop = df_merge.loc[
    #                 df_merge['pval_h3k4me3_promoter'] != 0, :]
    # df_corr = df_merge_drop.corr()
    # file_corr = '/local/zy/PEI/mid_data/cell_line/model_input/GM12878/corr.txt'
    # df_corr.to_csv(file_corr, sep='\t')
