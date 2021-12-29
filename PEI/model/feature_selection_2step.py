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
from sklearn.metrics import recall_score, precision_score
from glob import glob


def combine_boundary():
    path_boundary = \
        '/lustre/tianlab/zhangyu/PEI/origin_data/TAD/' \
        'primary_cohort_TAD_boundaries'
    path_out = '/lustre/tianlab/zhangyu/PEI/mid_data_correct/try_CTCF'
    files_boundary = glob(rf"{path_boundary}/*.bed")
    # files_boundary = [os.path.join(path_boundary, file)
    #                   for file in os.listdir(path_boundary)]
    file_combine = os.path.join(path_out, 'TAD_boundaries.combine.bed')
    cat_out = os.path.join(path_out, 'TAD_boundaries.cat.out')
    os.system(f"cat {' '.join(files_boundary)} > {cat_out}")
    os.system(f"bedtools sort -i {cat_out} | uniq > {file_combine}")
    os.remove(cat_out)

    df_combine = pd.read_csv(file_combine, sep='\t', header=None)
    df_combine.columns = ['chrom', 'start', 'end']
    for idx, file in enumerate(files_boundary):
        df_sub = pd.read_csv(file, sep='\t', header=None)
        df_sub.columns = ['chrom', 'start', 'end']
        df_sub[idx] = np.full(df_sub.shape[0], 1)
        df_combine = pd.merge(df_combine, df_sub, how='left',
                              on=['chrom', 'start', 'end'])
        df_combine = df_combine.fillna(0)

    df_combine['sum'] = np.sum(df_combine.iloc[:, 3:], axis=1)
    df_boundary = df_combine.loc[
        df_combine['sum'] >= 4, ['chrom', 'start', 'end']]
    file_out = os.path.join(path_out, 'TAD_boundaries.bed')
    df_boundary.to_csv(file_out, sep='\t', header=None, index=None)

    return


def annotate_insulator():
    file_cre_in = \
        '/lustre/tianlab/zhangyu/PEI/mid_data_correct/cell_line/' \
        'DHS/cRE_annotation/GM12878/DHS_promoter_H3K4me3_H3K27ac_CTCF.txt'
    file_boundary = '/lustre/tianlab/zhangyu/PEI/mid_data_correct/' \
                    'try_CTCF/TAD_boundaries.bed'
    file_cre_intersect = '/lustre/tianlab/zhangyu/PEI/mid_data_correct/' \
                         'try_CTCF/cRE_GM12878.intersect'
    file_cre_out = '/lustre/tianlab/zhangyu/PEI/mid_data_correct/' \
                   'try_CTCF/cRE_GM12878.txt'

    os.system(f"bedtools intersect -a {file_cre_in} -b {file_boundary} -wao "
              f"> {file_cre_intersect}")
    with open(file_cre_intersect, 'r') as r_f:
        with open(file_cre_out, 'w') as w_f:
            fmt_dhs = "{chrom}\t{start}\t{end}\t{dhs_id}\t{cre}\t" \
                      "{dhs_score}\t{promoter_id}\t{score_h3k4me3}\t" \
                      "{p_h3k4me3}\t{score_h3k27ac}\t{p_h3k27ac}\t" \
                      "{score_ctcf}\n"
            for line in r_f:
                list_line = line.strip().split('\t')
                chrom = list_line[0]
                start = list_line[1]
                end = list_line[2]
                dhs_id = list_line[3]
                dhs_score = float(list_line[4])
                promoter_id = list_line[5]
                score_h3k4me3 = float(list_line[6])
                p_h3k4me3 = float(list_line[7])
                score_h3k27ac = float(list_line[8])
                p_h3k27ac = float(list_line[9])
                score_ctcf = float(list_line[10])
                overlap_tad = float(list_line[14])
                if (promoter_id != '.') & (p_h3k4me3 != 0) & \
                        (p_h3k27ac != 0):
                    cre = 'Protein-Promoter(Enhancer)'
                elif (promoter_id == '.') & (p_h3k4me3 != 0) & \
                        (p_h3k27ac != 0):
                    cre = 'Other-Promoter(Enhancer)'
                elif (promoter_id != '.') & (p_h3k4me3 != 0) & \
                        (p_h3k27ac == 0):
                    cre = 'Protein-Promoter'
                elif (p_h3k4me3 == 0) & (p_h3k27ac != 0):
                    cre = 'Enhancer'
                elif (p_h3k27ac == 0) & (score_ctcf != 0) & (overlap_tad != 0):
                    cre = 'Insulator'
                else:
                    cre = '.'
                w_f.write(fmt_dhs.format(**locals()))

    return


def test_one_cell():
    path_out = '/lustre/tianlab/zhangyu/PEI/mid_data_correct/try_CTCF'

    term = 'GM12878'

    file_pre = '/lustre/tianlab/zhangyu/PEI/mid_data_correct/cell_line/' \
               'model_input/GM12878/correlation.txt'
    file_cre = os.path.join(path_out, 'cRE_GM12878.txt')
    df_pre = pd.read_csv(file_pre, sep='\t')
    abs_distance = np.abs(df_pre['distance'])
    df_pre = df_pre.loc[abs_distance > 5000, :]

    file_promoter = path_origin + \
                    '/gene/promoters.up2k.protein.gencode.v19.unique.bed'
    df_promoter = pd.read_csv(file_promoter, sep='\t', header=None)

    file_tmp = os.path.join(path_out, "input_file.tmp")
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

    file_cre_ctcf = os.path.join(path_out, "cRE.CTCF.tmp")
    os.system(f"grep -w 'Insulator' {file_cre} > {file_cre_ctcf}")
    file_ctcf_tmp = os.path.join(path_out, "CTCF.tmp")
    os.system(f"bedtools intersect -a {file_tmp} -b {file_cre_ctcf} -wao | "
              f"cut -f 4,5,9,10,11,17 > {file_ctcf_tmp}")
    df_ctcf = pd.read_csv(file_ctcf_tmp, sep='\t', header=None, na_values='.',
                          dtype={2: 'str', 3: 'str', 4: 'float', 5: 'float'})

    def unique_ctcf(df_in):
        if df_in.shape[0] == 1:
            if np.isnan(df_in.iloc[0, 5]):
                df_in.iloc[0, 5] = 0
                df_out = df_in.loc[:, [0, 1, 5]]
            else:
                df_out = df_in.loc[:, [0, 1, 5]]
        else:
            max_ctcf = np.max(df_in[5])
            # df_out = df_in.loc[df_in[5] == max_ctcf, [0, 1, 4, 5]]
            df_out = df_in.loc[df_in[5] == max_ctcf, [0, 1, 5]]

        return df_out

    df_uniq = df_ctcf.groupby([0, 1]).apply(unique_ctcf)
    df_uniq.index = list(range(df_uniq.shape[0]))
    # df_uniq.columns = \
    #     ['gene', 'dhs_id', 'score_dhs_insulator', 'score_ctcf_insulator']
    df_uniq.columns = ['gene', 'dhs_id', 'score_ctcf_insulator']
    df_uniq = df_uniq.drop_duplicates()
    df_genome_ctcf = pd.merge(df_pre, df_uniq, on=['gene', 'dhs_id'],
                              how='left')
    df_genome_ctcf = df_genome_ctcf.fillna(0)

    os.remove(file_tmp)
    os.remove(file_cre_ctcf)
    os.remove(file_ctcf_tmp)

    file_out = os.path.join(path_out, f"{term}_input_file.txt")
    df_genome_ctcf.to_csv(file_out, sep='\t', index=None, na_rep='NA')

    file_corr = file_out
    file_label = os.path.join(path_label, f"{term}/{term}.txt")
    file_fea_label = os.path.join(path_out, f"{term}_feature_label.txt")
    file_res = os.path.join(path_out, f"{term}_result.txt")
    label = term

    df_corr = pd.read_csv(file_corr, sep='\t')
    df_label = pd.read_csv(file_label, sep='\t')
    # only select distal enhancer
    df_corr = \
        df_corr.loc[df_corr['type_cre'] != 'Protein-Promoter(Enhancer)', ]
    # if df_label.shape[0] == 0:
    #     return
    df_label['label'] = np.full(df_label.shape[0], 1)

    df_combine = pd.merge(
        df_corr, df_label, how='left',
        on=['gene', 'dhs_id', 'type_cre', 'ref_dhs_id'])
    df_combine = df_combine.fillna(0)

    # first step
    array_pred = np.full(df_combine.shape[0], 0)
    array_pred[df_combine['score_ctcf_insulator'] <= 0] = 1
    df_combine['pred'] = array_pred
    df_combine.to_csv(file_fea_label, sep='\t', index=None)
    precision = precision_score(df_combine['label'], array_pred)
    recall = recall_score(df_combine['label'], array_pred)

    cols = df_combine.columns[5:-2]
    list_res = [{'feature': 'CTCF_pred', 'correlation': '',
                 'diff_median': precision, 'pval': recall,
                 'label': label}]
    df_combine_filter = df_combine.loc[df_combine['pred'] == 1, :]
    for col in cols:
        df_sub = df_combine_filter.loc[:, [col, 'label']]
        array_pos = df_sub.loc[df_sub['label'] == 1, col]
        array_neg = df_sub.loc[df_sub['label'] == 0, col]
        try:
            _, pval = mannwhitneyu(array_pos, array_neg, alternative='greater')
        except ValueError:
            pval = np.nan
        diff_median = np.median(array_pos) - np.median(array_neg)
        # feature, corr = col.split('|')
        list_res.append({'feature': col, 'correlation': 'Spearman',
                         'diff_median': diff_median, 'pval': pval,
                         'label': label})
    df_res = pd.DataFrame(list_res)
    df_res.to_csv(file_res, sep='\t', index=False, na_rep='NA')

    return


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

    # first step
    array_pred = np.full(df_combine.shape[0], 0)
    array_pred[df_combine['score_ctcf_insulator'] <= 0] = 1
    df_combine['pred'] = array_pred
    df_combine.to_csv(file_fea_label, sep='\t', index=None)
    precision = precision_score(df_combine['label'], array_pred)
    recall = recall_score(df_combine['label'], array_pred)

    cols = df_combine.columns[5:-2]
    list_res = [{'feature': 'CTCF_pred', 'correlation': '',
                 'diff_median': precision, 'pval': recall,
                 'label': label}]
    df_combine_filter = df_combine.loc[df_combine['pred'] == 1, :]
    for col in cols:
        df_sub = df_combine_filter.loc[:, [col, 'label']]
        array_pos = df_sub.loc[df_sub['label'] == 1, col]
        array_neg = df_sub.loc[df_sub['label'] == 0, col]
        try:
            _, pval = mannwhitneyu(array_pos, array_neg, alternative='greater')
        except ValueError:
            pval = np.nan
        diff_median = np.median(array_pos) - np.median(array_neg)
        # feature, corr = col.split('|')
        list_res.append({'feature': col, 'correlation': 'Spearman',
                         'diff_median': diff_median, 'pval': pval,
                         'label': label})
    df_res = pd.DataFrame(list_res)
    df_res.to_csv(file_res, sep='\t', index=False, na_rep='NA')

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
    file_assay_out = os.path.join(path_label, 'Assay_2step_comparation_1_0.txt')
    df_assay.to_csv(file_assay_out, sep='\t', index=None, na_rep='NA')

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
    file_term_out = os.path.join(path_label, 'Cell_2step_comparation_1_0.txt')
    df_term.to_csv(file_term_out, sep='\t', index=None, na_rep='NA')

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

    time_end = time()
    print(time_end - time_start)
