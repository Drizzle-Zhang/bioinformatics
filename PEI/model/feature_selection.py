#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: feature_selection.py
# @time: 6/1/20 6:15 PM

from time import time
import pandas as pd
import numpy as np
import os
from scipy.stats import mannwhitneyu
from multiprocessing import Pool


def combine_corr_label(dict_in):
    file_corr = dict_in['file_corr']
    file_label = dict_in['file_label']
    file_fea_label = dict_in['file_fea_label']
    file_res = dict_in['file_res']
    label = dict_in['label']

    df_corr = pd.read_csv(file_corr, sep='\t')
    df_label = pd.read_csv(file_label, sep='\t')
    # only select distal enhancer
    df_corr = \
        df_corr.loc[df_corr['type_cre'] != 'Protein-Promoter(Enhancer)', ]
    if df_label.shape[0] == 0:
        return
    df_label['label'] = np.full(df_label.shape[0], 1)

    df_combine = pd.merge(
        df_corr, df_label, how='left',
        on=['gene', 'dhs_id', 'type_cre', 'ref_dhs_id'])
    df_combine = df_combine.fillna(0)
    df_combine.to_csv(file_fea_label, sep='\t', index=None)

    cols = df_combine.columns[4:-1]
    list_res = []
    for col in cols:
        df_sub = df_combine.loc[:, [col, 'label']]
        array_pos = df_sub.loc[df_sub['label'] == 1, col]
        array_neg = df_sub.loc[df_sub['label'] == 0, col]
        _, pval = mannwhitneyu(array_pos, array_neg, alternative='greater')
        diff_median = np.median(array_pos) - np.median(array_neg)
        # feature, corr = col.split('|')
        list_res.append({'feature': col, 'correlation': 'Spearman',
                         'diff_median': diff_median, 'pval': pval,
                         'label': label})
    df_res = pd.DataFrame(list_res)
    df_res.to_csv(file_res, sep='\t', index=False)

    return df_res


def test_features_cell():
    file_meta = os.path.join(path_label, 'meta_label.txt')
    df_meta = pd.read_csv(file_meta, sep='\t')
    df_meta_assay = df_meta.drop_duplicates(subset='Assay')

    list_assay = []
    for sub_dict in df_meta_assay.to_dict('records'):
        term = sub_dict['Biosample term name']
        label = sub_dict['Assay']
        file_corr = os.path.join(path_model_input, f"{term}/input_file.txt")
        file_label = os.path.join(path_label, f"{term}/{label}.txt")
        file_fea_label = os.path.join(
            path_model_input, f"{term}/{label}_feature_label.txt")
        file_res = os.path.join(
            path_model_input, f"{term}/{label}_result.txt")
        if not os.path.exists(file_corr):
            continue
        list_assay.append({
            'term': term, 'label': label, 'file_corr': file_corr,
            'file_label': file_label, 'file_fea_label': file_fea_label,
            'file_res': file_res
        })

    pool = Pool(processes=40)
    list_df = pool.map(combine_corr_label, list_assay)
    pool.close()
    df_assay = pd.concat(list_df)
    file_assay_out = os.path.join(path_label, 'Assay_comparation_1_0.txt')
    df_assay.to_csv(file_assay_out, sep='\t', index=None)

    terms = df_meta_assay['Biosample term name'].unique().tolist()
    list_term = []
    for term in terms:
        file_corr = os.path.join(path_model_input, f"{term}/input_file.txt")
        file_label = os.path.join(path_label, f"{term}/{term}.txt")
        file_fea_label = os.path.join(
            path_model_input, f"{term}/{term}_feature_label.txt")
        file_res = os.path.join(
            path_model_input, f"{term}/{term}_result.txt")
        list_term.append({
            'label': term, 'file_corr': file_corr,
            'file_label': file_label, 'file_fea_label': file_fea_label,
            'file_res': file_res
        })

    pool = Pool(processes=40)
    list_df = pool.map(combine_corr_label, list_term)
    pool.close()
    df_term = pd.concat(list_df)
    file_term_out = os.path.join(path_label, 'Cell_comparation_1_0.txt')
    df_term.to_csv(file_term_out, sep='\t', index=None)

    return


if __name__ == '__main__':
    time_start = time()
    # path_root = '/local/zy/PEI'
    path_root = '/lustre/tianlab/zhangyu/PEI'
    path_origin = path_root + '/origin_data'
    # path_mid = path_root + '/mid_data'
    path_mid = path_root + '/mid_data_correct'

    # high-resolution
    path_model_input = path_mid + '/cell_line/model_input'
    path_label = \
        path_mid + '/training_label/label_interactions'
    test_features_cell()

    # example
    file_corr_GM12878 = \
        path_mid + '/cell_line/model_input/GM12878/correlation.txt'
    file_label_GM12878 = \
        path_mid + '/training_label/label_interactions/GM12878/GM12878.txt'
    file_out_GM12878 = \
        path_mid + '/cell_line/model_input/GM12878/selection.txt'
    file_res_GM12878 = \
        path_mid + '/cell_line/model_input/GM12878/comparation_1_0.txt'
    # dict_in = {
    #     'label': term, 'file_corr': file_corr,
    #     'file_label': file_label, 'file_fea_label': file_fea_label,
    #     'file_res': file_res
    # }
    time_end = time()
    print(time_end - time_start)
