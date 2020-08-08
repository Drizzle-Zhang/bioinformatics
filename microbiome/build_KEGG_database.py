#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: build_KEGG_database.py
# @time: 7/26/20 1:04 AM

from time import time
import pandas as pd
import numpy as np
from multiprocessing import Pool
from functools import partial


def sub_kegg(df_kegg, level, term):
    array_bool = np.full(df_kegg.shape[0], 0)
    for i, sub_dict in enumerate(df_kegg.to_dict('records')):
        pathways = sub_dict['metadata_KEGG_Pathways']
        for pathway in pathways.strip().split('|'):
            list_term = pathway.split(';')
            if len(list_term) == 3:
                if list_term[0] == 'Unclassified':
                    if level == 3:
                        continue
                    else:
                        term_level = list_term[level]
                else:
                    term_level = list_term[level - 1]
                if term_level == term:
                    array_bool[i] = 1
                    break

    series_out = pd.Series(array_bool, index=df_kegg.index, name=term)

    return series_out


def build_kegg_db():
    file_kegg = '/home/drizzle_zhang/microbiome/result/9.PICRUSt/' \
                'origin_data/metadata_KEGG.txt'
    df_kegg = pd.read_csv(file_kegg, sep='\t', index_col=0)

    # KEGG level3
    level = 3
    file_l3 = '/home/drizzle_zhang/microbiome/result/9.PICRUSt/ko_predict/' \
              'KEGG_L3.txt'
    df_l3 = pd.read_csv(file_l3, sep='\t', index_col=0)
    list_kegg_l3 = list(df_l3.index)

    # generate KO-KEGG file
    pool = Pool(6)
    func_kegg = partial(sub_kegg, df_kegg, level)
    list_out = pool.map(func_kegg, list_kegg_l3)
    pool.close()

    df_db = pd.concat(list_out, axis=1)
    df_db = df_db.loc[:, np.sum(df_db) != 0]
    file_da_l3 = '/home/drizzle_zhang/microbiome/result/9.PICRUSt/' \
                 'origin_data/KO_KEGG_L3.txt'
    df_db.to_csv(file_da_l3, sep='\t')

    # generate gene set file
    file_gmt_l3 = '/home/drizzle_zhang/microbiome/result/9.PICRUSt/' \
                  'origin_data/KEGG_KO_L3.gmt'
    index_ko = df_kegg.index
    with open(file_gmt_l3, 'w') as w_gmt:
        for term in list_kegg_l3:
            list_term_ko = []
            list_l2 = []
            list_l1 = []
            for i, sub_dict in enumerate(df_kegg.to_dict('records')):
                pathways = sub_dict['metadata_KEGG_Pathways']
                for pathway in pathways.strip().split('|'):
                    list_term = pathway.split(';')
                    if len(list_term) == 3:
                        if list_term[0] == 'Unclassified':
                            if level == 3:
                                continue
                            else:
                                term_level = list_term[level]
                        else:
                            term_level = list_term[level - 1]
                        if term_level == term:
                            list_term_ko.append(index_ko[i])
                            if list_term[0] == 'Unclassified':
                                list_l1.append(list_term[level - 1])
                            else:
                                list_l1.append(list_term[level - 3])
                                list_l2.append(list_term[level - 2])
                            break
            if len(list_term_ko) == 0:
                continue
            str_gene_set = '\t'.join(list_term_ko)
            term_l1 = list(set(list_l1))
            term_l2 = list(set(list_l2))
            line_w = f"{term}\t{';'.join(term_l1)}|{';'.join(term_l2)}\t" \
                     f"{str_gene_set}\n"
            w_gmt.write(line_w)

    # KEGG level2
    level = 2
    file_l2 = '/home/drizzle_zhang/microbiome/result/9.PICRUSt/ko_predict/' \
              'KEGG_L2.txt'
    df_l2 = pd.read_csv(file_l2, sep='\t', index_col=0)
    list_kegg_l2 = list(df_l2.index)

    # generate KO-KEGG file
    pool = Pool(6)
    func_kegg = partial(sub_kegg, df_kegg, level)
    list_out = pool.map(func_kegg, list_kegg_l2)
    pool.close()

    df_db = pd.concat(list_out, axis=1)
    df_db = df_db.loc[:, np.sum(df_db) != 0]
    file_da_l2 = '/home/drizzle_zhang/microbiome/result/9.PICRUSt/' \
                 'origin_data/KO_KEGG_L2.txt'
    df_db.to_csv(file_da_l2, sep='\t')

    # generate gene set file
    file_gmt_l2 = '/home/drizzle_zhang/microbiome/result/9.PICRUSt/' \
                  'origin_data/KEGG_KO_L2.gmt'
    index_ko = df_kegg.index
    with open(file_gmt_l2, 'w') as w_gmt:
        for term in list_kegg_l2:
            list_term_ko = []
            list_l1 = []
            for i, sub_dict in enumerate(df_kegg.to_dict('records')):
                pathways = sub_dict['metadata_KEGG_Pathways']
                for pathway in pathways.strip().split('|'):
                    list_term = pathway.split(';')
                    if len(list_term) == 3:
                        if list_term[0] == 'Unclassified':
                            if level == 3:
                                continue
                            else:
                                term_level = list_term[level]
                        else:
                            term_level = list_term[level - 1]
                        if term_level == term:
                            list_term_ko.append(index_ko[i])
                            if list_term[0] == 'Unclassified':
                                list_l1.append(list_term[level - 1])
                            else:
                                list_l1.append(list_term[level - 2])
                            break
            if len(list_term_ko) == 0:
                continue
            str_gene_set = '\t'.join(list_term_ko)
            term_l1 = list(set(list_l1))
            line_w = f"{term}\t{';'.join(term_l1)}\t{str_gene_set}\n"
            w_gmt.write(line_w)

    return


if __name__ == '__main__':
    time_start = time()

    time_end = time()
    print(time_end - time_start)