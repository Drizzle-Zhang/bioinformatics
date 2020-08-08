#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: find_KO.py
# @time: 7/25/20 7:49 PM

from time import time
import os
import pandas as pd
import numpy as np


def prepare_mat():
    file_ko_otu = '/lustre/tianlab/zhangyu/microbiome/picrust/picrust/' \
                  'data/ko_13_5_precalculated.tab'
    df_ko_otu = pd.read_csv(file_ko_otu, sep='\t', index_col=0)

    path_picrust = \
        '/lustre/tianlab/zhangyu/microbiome/result/9.PICRUSt/origin_data'
    df_ko_otu_num = df_ko_otu.iloc[:-2, :-1]
    file_ko_otu_num = os.path.join(path_picrust, 'mat_OTU_KO.txt')
    df_ko_otu_num.to_csv(file_ko_otu_num, sep='\t')

    df_nsti = df_ko_otu.iloc[:-2, -1]
    file_nsti = os.path.join(path_picrust, 'metadata_NSTI.txt')
    df_nsti.to_csv(file_nsti, sep='\t')

    df_kegg = df_ko_otu.iloc[-2:, :-1].T
    file_kegg = os.path.join(path_picrust, 'metadata_KEGG.txt')
    df_kegg.to_csv(file_kegg, sep='\t')

    return


def find_ko():
    path_picrust = \
        '/lustre/tianlab/zhangyu/microbiome/result/9.PICRUSt/origin_data'
    path_out = '/lustre/tianlab/zhangyu/microbiome/result/9.PICRUSt/' \
               'ko_predict/select_KOs'
    file_otu_taxonomy = '/lustre/tianlab/zhangyu/microbiome/result/2.OTUs/' \
                        'OTUs_stat/taxonomy.txt'
    file_ko_otu_num = os.path.join(path_picrust, 'mat_OTU_KO.txt')
    sub_family = 'f__Lactobacillaceae'

    df_otu_tax = pd.read_csv(file_otu_taxonomy, sep='\t')
    sub_otu_tax = df_otu_tax.loc[df_otu_tax['family'] == sub_family, 'OTU_ID']

    # filter low quality OTUs
    file_nsti = os.path.join(path_picrust, 'metadata_NSTI.txt')
    df_nsti = pd.read_csv(file_nsti, sep='\t', index_col=0, header=None)
    sub_nsti = df_nsti.loc[sub_otu_tax.tolist()].dropna()
    sub_nsti = sub_nsti.loc[sub_nsti < 0.1]
    list_otu = list(sub_nsti.index)

    df_ko_otu_num = pd.read_csv(file_ko_otu_num, sep='\t',
                                index_col=0, dtype=np.int32)
    sub_ko_otu = df_ko_otu_num.loc[list_otu, :]
    sub_ko_otu = sub_ko_otu.dropna()

    file_sub_mat = os.path.join(path_out, f'mat_OTU_KO_{sub_family}.txt')
    sub_ko_otu.to_csv(file_sub_mat, sep='\t')

    sub_ko_otu = sub_ko_otu.loc[:, np.sum(sub_ko_otu) != 0]
    list_ko = list(sub_ko_otu.columns)

    file_kegg = os.path.join(path_picrust, 'metadata_KEGG.txt')
    df_kegg = pd.read_csv(file_kegg, sep='\t', index_col=0)
    sub_kegg = df_kegg[list_ko, :]

    file_sub = os.path.join(path_out, f'KO_{sub_family}.txt')
    sub_kegg.to_csv(file_sub, sep='\t')

    return


def prepare_diff_ko():
    path_picrust = \
        '/lustre/tianlab/zhangyu/microbiome/result/9.PICRUSt/origin_data'
    path_out = '/lustre/tianlab/zhangyu/microbiome/result/9.PICRUSt/' \
               'ko_predict/select_KOs'

    file_ko_otu_num = os.path.join(path_picrust, 'mat_OTU_KO.txt')
    df_ko_otu_num = pd.read_csv(file_ko_otu_num, sep='\t',
                                index_col=0, dtype=np.int32)
    file_otu_taxonomy = '/lustre/tianlab/zhangyu/microbiome/result/2.OTUs/' \
                        'OTUs_stat/taxonomy.txt'
    df_otu_tax = pd.read_csv(file_otu_taxonomy, sep='\t')
    file_nsti = os.path.join(path_picrust, 'metadata_NSTI.txt')
    df_nsti = pd.read_csv(file_nsti, sep='\t', index_col=0, header=None)
    file_family_top30 = '/lustre/tianlab/zhangyu/microbiome/result/' \
                        '3.Community_Structure/barplot/Top30_family.csv'
    df_family_top30 = pd.read_csv(file_family_top30, sep=',')
    families_top10 = \
        [f"f__{family}"
         for family in (df_family_top30.loc[0:9, 'Taxonomy']).tolist()]

    # family level
    # families = \
    #     list(set(families_top20).intersection(
    #         set(df_otu_tax['family'].unique().tolist())
    #     ))
    list_exp = []
    for sub_family in families_top10:
        sub_otu_tax = df_otu_tax.loc[
            df_otu_tax['family'] == sub_family, 'OTU_ID']
        try:
            sub_nsti = df_nsti.loc[sub_otu_tax.tolist()].dropna()
        except KeyError:
            continue
        # sub_nsti = sub_nsti.loc[sub_nsti[1] < 0.06]
        list_otu = list(sub_nsti.index)
        sub_ko_otu = df_ko_otu_num.loc[list_otu, :]
        sub_ko_otu = sub_ko_otu.dropna()
        series_family = np.sum(sub_ko_otu)
        series_family.name = sub_family
        list_exp.append(series_family)

    df_exp = pd.concat(list_exp, axis=1)
    file_family_mat = os.path.join(path_out, f'mat_family_KO.txt')
    df_exp.to_csv(file_family_mat, sep='\t')

    return


if __name__ == '__main__':
    time_start = time()

    time_end = time()
    print(time_end - time_start)
