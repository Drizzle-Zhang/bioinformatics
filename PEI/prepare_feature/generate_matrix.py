#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: generate_matrix.py
# @time: 2020/5/3 22:55

from time import time
import pandas as pd
import os
import numpy as np
from multiprocessing import Pool
from functools import partial
from statsmodels.api import distributions


def generate_promoter_file():
    os.system(
        f"bedtools intersect -wa -wb -a {file_all_index} -b {file_promoter} | "
        f"cut -f 4,5,9 > {file_dhs_promoter}")
    df_promoter = pd.read_csv(file_promoter, sep='\t', header=None)
    df_promoter['idx'] = df_promoter.index
    df_promoter = df_promoter.loc[:, [3, 'idx']]
    df_promoter['gene'] = df_promoter[3].apply(
        lambda x: x.split('<-')[0])
    df_promoter = df_promoter.drop_duplicates(subset='gene')
    df_dhs_promoter = pd.read_csv(file_dhs_promoter, sep='\t', header=None)
    # print(df_dhs_promoter.shape[0])
    # df_dhs_promoter = df_dhs_promoter.drop_duplicates(subset=[0, 1, 'gene'])
    df_dhs_promoter = pd.merge(df_dhs_promoter, df_promoter,
                               left_on=2, right_on=3)
    # print(df_dhs_promoter.shape[0])
    df_dhs_promoter = df_dhs_promoter.rename(
        columns={0: 'dhs_id', 1: 'dhs_idx', 2: 'gene_pos'})
    df_dhs_promoter.to_csv(
        file_dhs_promoter, sep='\t', index=None,
        columns=['dhs_idx', 'dhs_id', 'gene_pos', 'gene', 'idx'])

    return


def sub_dhs_term(df_all, path_tmp, dict_in):
    file_index = dict_in['file_index']
    term = dict_in['term']
    str_term = dict_in['str_term']
    df_term = pd.read_csv(file_index, sep='\t', header=None,
                          index_col=1, usecols=[0, 2], names=[term, ''])
    # print(np.min(df_term))
    # print(np.max(df_term))
    df_merge = pd.merge(df_all, df_term, how='outer',
                        left_index=True, right_index=True)
    if na_mode == 'minus':
        df_merge = df_merge.fillna(np.min(df_term)[0] - 1)
    elif na_mode == 'constant':
        df_merge = df_merge.fillna(0)
    file_tmp = os.path.join(path_tmp, str_term + '.txt')
    if df_merge.shape[0] != df_all.shape[0]:
        print(term)
        print(df_merge.shape[0])
        print(df_all.shape[0])
    df_merge.to_csv(file_tmp, sep='\t', index=None, columns=[term])

    return file_tmp


def dhs_matrix():
    df_all = pd.read_csv(
        file_all_index, sep='\t', header=None, usecols=[3], names=['dhs_id'])
    path_tmp = os.path.join(path_matrix, 'tmp_DHS')
    if not os.path.exists(path_tmp):
        os.mkdir(path_tmp)
    file_dhs_id = os.path.join(path_tmp, 'dhs_id.txt')
    df_all.to_csv(file_dhs_id, sep='\t', index=None, columns=['dhs_id'])

    # cell line
    df_meta_cell = pd.read_csv(
        os.path.join(path_dhs_cell, 'meta.reference.tsv'), sep='\t')
    list_tmp = [file_dhs_id]
    list_input = []
    for term in (df_meta_cell['Biosample term name'].unique()).tolist():
        str_term = term.replace(' ', '_').replace('/', '+')
        path_term = os.path.join(path_dhs_cell, str_term)
        file_index = os.path.join(path_term, 'ref_index.txt')
        list_input.append({'file_index': file_index, 'term': term,
                           'str_term': str_term})

    # tissue
    df_meta_tissue = pd.read_csv(
        os.path.join(path_dhs_tissue_stan, 'meta.reference.tsv'), sep='\t')
    for i in range(df_meta_tissue.shape[0]):
        life_organ = df_meta_tissue.loc[i, 'Biosample life_organ']
        term = df_meta_tissue.loc[i, 'Biosample term name']
        str_life_organ = life_organ.replace(' ', '_')
        str_term = term.replace(' ', '_').replace('/', '+').replace("'", "--")
        path_term = os.path.join(
            path_dhs_tissue_stan, f"{str_life_organ}/{str_term}")
        file_index = os.path.join(path_term, 'ref_index.txt')
        list_input.append({'file_index': file_index,
                           'term': f"{life_organ} | {term}",
                           'str_term': f"{str_life_organ}+{str_term}"})

    pool = Pool(num_cpu)
    func_dhs = partial(sub_dhs_term, df_all, path_tmp)
    list_res = pool.map(func_dhs, list_input)
    pool.close()

    list_tmp.extend(list_res)
    file_matrix_dhs = os.path.join(path_matrix, 'DHS_matrix.txt')
    os.system(f"paste {' '.join(list_tmp)} > {file_matrix_dhs}")

    os.system(f"rm -rf {path_tmp}")

    def calculate_score(df_in):
        df_num = df_in.iloc[:, 5:]
        df_out = np.max(df_num, axis=0)
        df_out['gene'] = df_in.iloc[0, 3]

        return df_out

    df_dhs = pd.read_csv(file_matrix_dhs, sep='\t', index_col=0)
    df_dhs_promoter = pd.read_csv(file_dhs_promoter, sep='\t')
    df_dhs_promoter = pd.merge(
        df_dhs_promoter, df_dhs, left_on='dhs_id', right_index=True)
    df_gene_dhs = df_dhs_promoter.groupby('idx').apply(calculate_score)
    df_gene_dhs.index = df_gene_dhs['gene']
    df_gene_dhs = df_gene_dhs.rename(columns={'gene': 'gene_redu'})
    df_gene_dhs = df_gene_dhs.drop('gene_redu', axis=1)
    df_gene_dhs.to_csv(
        os.path.join(path_matrix, 'DHS_matrix.promoter.txt'), sep='\t')

    return


def sub_h3k4me3(df_dhs_promoter, path_tmp, df_all_pro, dict_in):
    file_index = dict_in['file_index']
    path_h3k4me3 = dict_in['path_h3k4me3']
    file_h3k4me3 = dict_in['file_h3k4me3']
    term = dict_in['term']
    str_term = dict_in['str_term']
    df_index = pd.read_csv(
        file_index, sep='\t', header=None, usecols=[0, 2, 3],
        names=['dhs_id_term', 'dhs_id', 'dhs_idx']
    )
    df_index_gene = pd.merge(
        df_dhs_promoter, df_index, on=['dhs_id', 'dhs_idx'], how='left')
    df_h3k4me3 = pd.read_csv(
        file_h3k4me3, sep='\t', header=None, usecols=[3, 6],
        names=['dhs_id_term', 'h3k4me3_score'], na_values=-10000)
    df_index_gene_h3k4me3 = pd.merge(
        df_index_gene, df_h3k4me3, on='dhs_id_term', how='left')
    df_index_gene_h3k4me3 = df_index_gene_h3k4me3.sort_values('idx')
    df_index_gene_h3k4me3.to_csv(
        os.path.join(path_h3k4me3, 'gene_h3k4me3.txt'), sep='\t', index=None,
        na_rep='NA'
    )
    df_gene_h3k4me3 = \
        df_index_gene_h3k4me3.loc[:, ['idx', 'h3k4me3_score']]

    def calculate_score(df_in):
        if df_in.shape[0] == 1:
            return df_in
        else:
            max_score = np.max(df_in.loc[:, 'h3k4me3_score'], )
            row_out = df_in.loc[df_in.loc[:, 'h3k4me3_score'] == max_score, :]
            return row_out

    df_gene_h3k4me3['key'] = df_gene_h3k4me3['idx']
    print(term)
    if na_mode == 'minus':
        print(np.nanmin(df_gene_h3k4me3['h3k4me3_score']))
        df_gene_h3k4me3 = df_gene_h3k4me3.fillna(
            np.nanmin(df_gene_h3k4me3['h3k4me3_score']) - 2)
    elif na_mode == 'constant':
        df_gene_h3k4me3 = df_gene_h3k4me3.fillna(0)
    df_gene_uniq = df_gene_h3k4me3.groupby('key').apply(calculate_score)
    df_gene_uniq = df_gene_uniq.drop_duplicates()
    df_gene_uniq = df_gene_uniq.sort_values('idx')
    file_tmp = os.path.join(path_tmp, str_term + '.txt')
    df_gene_uniq = df_gene_uniq.rename(columns={'h3k4me3_score': term})
    if df_gene_uniq.shape[0] != df_all_pro.shape[0]:
        print(term)
        print(df_gene_uniq.shape[0])
        print(df_all_pro.shape[0])
    else:
        df_gene_uniq.to_csv(file_tmp, sep='\t', index=None, columns=[term])

    # df_out = df_gene_uniq.copy()
    # df_out.index = df_out['idx']
    # df_out = df_out[term]

    return file_tmp


def h3k4me3_matrix():
    df_dhs_promoter = pd.read_csv(file_dhs_promoter, sep='\t')

    path_tmp = os.path.join(path_matrix, 'tmp_H3K4me3')
    if not os.path.exists(path_tmp):
        os.mkdir(path_tmp)

    df_all_pro = df_dhs_promoter.loc[:, ['idx', 'gene']].copy()
    df_all_pro = df_all_pro.drop_duplicates()
    df_all_pro = df_all_pro.sort_values('idx')
    file_tmp_pro = os.path.join(path_tmp, 'promoter_idx.txt')
    df_all_pro.to_csv(file_tmp_pro, sep='\t', index=None, columns=['gene'])

    # cell line
    df_meta_cell = pd.read_csv(
        os.path.join(path_h3k4me3_cell, 'meta.reference.tsv'), sep='\t')
    list_tmp = [file_tmp_pro]
    list_input = []
    for term in (df_meta_cell['Biosample term name'].unique()).tolist():
        str_term = term.replace(' ', '_').replace('/', '+')
        path_term_dhs = os.path.join(path_dhs_cell, str_term)
        file_index = os.path.join(path_term_dhs, 'index.txt')
        if not os.path.exists(file_index):
            continue
        path_term_h3k4me3 = os.path.join(path_h3k4me3_cell, str_term)
        file_h3k4me3 = os.path.join(
            path_term_h3k4me3, 'DHS_promoter_H3K4me3.txt')
        if not os.path.exists(file_h3k4me3):
            continue
        list_input.append({'file_index': file_index, 'term': term,
                           'str_term': str_term, 'file_h3k4me3': file_h3k4me3,
                           'path_h3k4me3': path_term_h3k4me3})

    # tissue
    df_meta_tissue = pd.read_csv(
        os.path.join(path_h3k4me3_tissue, 'meta.reference.tsv'), sep='\t')
    df_meta_tissue = df_meta_tissue.loc[df_meta_tissue['Level'] == 'term', :]
    for i in df_meta_tissue.index:
        life_organ = df_meta_tissue.loc[i, 'Biosample life_organ']
        suborgan = df_meta_tissue.loc[i, 'Biosample suborgan']
        term = df_meta_tissue.loc[i, 'Biosample term name']
        file_ref_dhs = df_meta_tissue.loc[i, 'file_ref_dhs']
        str_life_organ = life_organ.replace(' ', '_')
        str_term = term.replace(' ', '_').replace('/', '+').replace("'", "--")
        path_term, _ = os.path.split(file_ref_dhs)
        file_index = os.path.join(path_term, 'index.txt')
        if suborgan == 'single':
            path_term_h3k4me3 = os.path.join(
                path_h3k4me3_tissue, f"{str_life_organ}/{str_term}")
        else:
            str_suborgan = suborgan.replace(' ', '_')
            path_term_h3k4me3 = os.path.join(
                path_h3k4me3_tissue,
                f"{str_life_organ}/{str_suborgan}/{str_term}")
        file_h3k4me3 = os.path.join(
            path_term_h3k4me3, 'DHS_promoter_H3K4me3.txt')
        list_input.append({'file_index': file_index,
                           'term': f"{life_organ} | {term}",
                           'str_term': f"{str_life_organ}+{str_term}",
                           'file_h3k4me3': file_h3k4me3,
                           'path_h3k4me3': path_term_h3k4me3})

    pool = Pool(num_cpu)
    func_h3k4me3 = partial(sub_h3k4me3, df_dhs_promoter, path_tmp, df_all_pro)
    list_res = pool.map(func_h3k4me3, list_input)
    pool.close()

    list_tmp.extend(list_res)
    file_matrix_h3k4me3 = os.path.join(path_matrix, 'H3K4me3_matrix.txt')
    os.system(f"paste {' '.join(list_tmp)} > {file_matrix_h3k4me3}")

    os.system(f"rm -rf {path_tmp}")

    return


def sub_h3k27ac(df_all, path_tmp, dict_in):
    file_index = dict_in['file_index']
    path_h3k27ac = dict_in['path_h3k27ac']
    file_h3k27ac = dict_in['file_h3k27ac']
    term = dict_in['term']
    str_term = dict_in['str_term']
    df_index = pd.read_csv(
        file_index, sep='\t', header=None, usecols=[0, 3],
        names=['dhs_id_term', 'dhs_idx']
    )
    df_index_dhs = pd.merge(
        df_all, df_index, left_index=True, right_on='dhs_idx')
    df_h3k27ac = pd.read_csv(
        file_h3k27ac, sep='\t', header=None, usecols=[3, 8],
        names=['dhs_id_term', 'h3k27ac_score'], na_values=-10000)
    df_index_dhs_h3k27ac = pd.merge(
        df_index_dhs, df_h3k27ac, on='dhs_id_term')
    df_index_dhs_h3k27ac = df_index_dhs_h3k27ac.sort_index()
    df_index_dhs_h3k27ac.to_csv(
        os.path.join(path_h3k27ac, 'dhs_h3k27ac.txt'),
        sep='\t', index=None, na_rep='.')
    df_dhs_h3k27ac = \
        df_index_dhs_h3k27ac.loc[:, ['dhs_idx', 'h3k27ac_score']]

    def calculate_score(df_in):
        if df_in.shape[0] == 1:
            return df_in
        else:
            max_score = np.max(df_in.loc[:, 'h3k27ac_score'])
            row_out = df_in.loc[df_in.loc[:, 'h3k27ac_score'] == max_score, :]
            return row_out

    df_dhs_h3k27ac['key'] = df_dhs_h3k27ac['dhs_idx']
    print(term)
    if na_mode == 'minus':
        print(np.nanmin(df_dhs_h3k27ac['h3k27ac_score']))
        df_dhs_h3k27ac = df_dhs_h3k27ac.fillna(
            np.nanmin(df_dhs_h3k27ac['h3k27ac_score']) - 2)
    elif na_mode == 'constant':
        df_dhs_h3k27ac = df_dhs_h3k27ac.fillna(0)
    df_dhs_uniq = df_dhs_h3k27ac.groupby('key').apply(calculate_score)
    df_dhs_uniq = df_dhs_uniq.drop_duplicates()
    df_dhs_uniq.index = df_dhs_uniq['dhs_idx']
    df_out = pd.merge(
        df_all, df_dhs_uniq, left_index=True, right_index=True, how='left')
    if na_mode == 'minus':
        df_out = df_out.fillna(np.nanmin(df_dhs_h3k27ac['h3k27ac_score']) - 2)
    elif na_mode == 'constant':
        df_out = df_out.fillna(0)
    # df_out = df_out.sort_values('dhs_idx')
    file_tmp = os.path.join(path_tmp, str_term + '.txt')
    df_out = df_out.rename(columns={'h3k27ac_score': term})
    if df_out.shape[0] != df_all.shape[0]:
        print(term)
        print(df_out.shape[0])
        print(df_all.shape[0])
    else:
        df_out.to_csv(file_tmp, sep='\t', index=None, columns=[term])

    return file_tmp


def h3k27ac_matrix():
    df_all = pd.read_csv(
        file_all_index, sep='\t', header=None, usecols=[3], names=['dhs_id'])
    path_tmp = os.path.join(path_matrix, 'tmp_H3K27ac')
    if not os.path.exists(path_tmp):
        os.mkdir(path_tmp)
    file_dhs_id = os.path.join(path_tmp, 'dhs_id.txt')
    df_all.to_csv(file_dhs_id, sep='\t', index=None, columns=['dhs_id'])

    # cell line
    df_meta_cell = pd.read_csv(
        os.path.join(path_h3k27ac_cell, 'meta.reference.tsv'), sep='\t')
    list_tmp = [file_dhs_id]
    list_input = []
    for term in (df_meta_cell['Biosample term name'].unique()).tolist():
        str_term = term.replace(' ', '_').replace('/', '+')
        path_term_dhs = os.path.join(path_dhs_cell, str_term)
        file_index = os.path.join(path_term_dhs, 'index.txt')
        if not os.path.exists(file_index):
            continue
        path_term_h3k27ac = os.path.join(path_h3k27ac_cell, str_term)
        file_h3k27ac = os.path.join(
            path_term_h3k27ac, 'DHS_promoter_H3K4me3_H3K27ac.txt')
        if not os.path.exists(file_h3k27ac):
            continue
        list_input.append({'file_index': file_index, 'term': term,
                           'str_term': str_term, 'file_h3k27ac': file_h3k27ac,
                           'path_h3k27ac': path_term_h3k27ac})

    # tissue
    df_meta_tissue_h3k4me3 = pd.read_csv(
        os.path.join(path_h3k4me3_tissue, 'meta.reference.tsv'), sep='\t')
    df_meta_tissue_h3k27ac = pd.read_csv(
        os.path.join(path_h3k27ac_tissue, 'meta.reference.tsv'), sep='\t')
    df_meta_tissue = pd.merge(
        df_meta_tissue_h3k27ac, df_meta_tissue_h3k4me3,
        on=['Biosample life_organ', 'Biosample suborgan',
            'Biosample term name', 'Level'])
    df_meta_tissue = df_meta_tissue.loc[df_meta_tissue['Level'] == 'term', :]
    for i in df_meta_tissue.index:
        life_organ = df_meta_tissue.loc[i, 'Biosample life_organ']
        suborgan = df_meta_tissue.loc[i, 'Biosample suborgan']
        term = df_meta_tissue.loc[i, 'Biosample term name']
        file_ref_dhs = df_meta_tissue.loc[i, 'file_ref_dhs']
        str_life_organ = life_organ.replace(' ', '_')
        str_term = term.replace(' ', '_').replace('/', '+').replace("'", "--")
        path_term, _ = os.path.split(file_ref_dhs)
        file_index = os.path.join(path_term, 'index.txt')
        if suborgan == 'single':
            path_term_h3k27ac = os.path.join(
                path_h3k27ac_tissue, f"{str_life_organ}/{str_term}")
        else:
            str_suborgan = suborgan.replace(' ', '_')
            path_term_h3k27ac = os.path.join(
                path_h3k27ac_tissue,
                f"{str_life_organ}/{str_suborgan}/{str_term}")
        file_h3k27ac = os.path.join(
            path_term_h3k27ac, 'DHS_promoter_H3K4me3_H3K27ac.txt')
        list_input.append({'file_index': file_index,
                           'term': f"{life_organ} | {term}",
                           'str_term': f"{str_life_organ}+{str_term}",
                           'file_h3k27ac': file_h3k27ac,
                           'path_h3k27ac': path_term_h3k27ac})

    pool = Pool(num_cpu)
    func_h3k27ac = partial(sub_h3k27ac, df_all, path_tmp)
    list_res = pool.map(func_h3k27ac, list_input)
    pool.close()

    list_tmp.extend(list_res)
    file_matrix_h3k27ac = os.path.join(path_matrix, 'H3K27ac_matrix.txt')
    os.system(f"paste {' '.join(list_tmp)} > {file_matrix_h3k27ac}")

    os.system(f"rm -rf {path_tmp}")

    return


def gtex_expression_matrix():
    mat_exp = pd.read_csv(file_expression, sep='\t', index_col=0)
    df_promoter = pd.read_csv(file_promoter_uniq, sep='\t', header=None)
    df_dhs_promoter = pd.read_csv(file_dhs_promoter, sep='\t')

    df_all_pro = df_dhs_promoter.loc[:, ['idx', 'gene']].copy()
    df_all_pro = df_all_pro.drop_duplicates()
    df_all_pro = pd.merge(df_all_pro, df_promoter, left_on=['idx', 'gene'],
                          right_on=[6, 7])
    df_all_pro = df_all_pro.sort_values('idx')

    df_ensg_gene = df_all_pro.loc[:, [4, 'gene']]
    df_pro_exp = pd.merge(df_ensg_gene, mat_exp, left_on=4, right_index=True)

    # ensg_id = df_all_pro[4].tolist()
    # genes = df_all_pro['gene']
    # mat_exp = mat_exp.loc[ensg_id, :]
    # gtex_genes = mat_exp['Description']
    # set_overlap = set(genes.tolist()).intersection(gtex_genes.tolist())

    df_pro_exp.index = df_pro_exp['gene']
    df_pro_exp = df_pro_exp.rename(columns={'gene': 'del_gene'})
    df_pro_exp = df_pro_exp.drop([4, 'del_gene', 'Description'], axis=1)

    # normalization
    def normalize(col_in):
        col_0 = col_in[col_in == 0]
        col_num = col_in[col_in != 0]
        ecdf = distributions.ECDF(col_num)
        col_ecdf = pd.Series(ecdf(col_num),
                             index=col_num.index, name=col_num.name)
        col_out = pd.concat([col_ecdf, col_0])

        return col_out

    df_quantile = df_pro_exp.apply(normalize)
    file_matrix_gtex_expression = \
        os.path.join(path_matrix, 'GTEx_expression_matrix.txt')
    df_quantile.to_csv(file_matrix_gtex_expression, sep='\t')

    return


if __name__ == '__main__':
    time_start = time()
    num_cpu = 40
    # na_mode = 'minus'
    na_mode = 'constant'
    path_root = '/lustre/tianlab/zhangyu/PEI'
    path_origin = path_root + '/origin_data'
    path_mid = path_root + '/mid_data_correct'

    file_all_index = \
        path_mid + '/database_feature/DHS_index/all_index.txt'
    path_matrix = path_mid + '/database_feature/matrix'
    file_promoter = \
        path_origin + '/gene/promoters.up2k.protein.gencode.v19.bed'
    file_dhs_promoter = \
        path_mid + '/database_feature/DHS_index/promoter_index.txt'
    generate_promoter_file()

    file_promoter_uniq = \
        path_origin + '/gene/promoters.up2k.protein.gencode.v19.unique.bed'
    df_promoter_main = pd.read_csv(file_promoter, sep='\t', header=None)
    df_promoter_main['idx'] = df_promoter_main.index
    df_promoter_main['gene'] = df_promoter_main[3].apply(
        lambda x: x.split('<-')[0])
    df_promoter_main = df_promoter_main.drop_duplicates(subset='gene')
    df_promoter_main.to_csv(
        file_promoter_uniq, sep='\t', index=None, header=None)

    # DHS
    path_dhs_cell = path_mid + '/cell_line/DHS/GRCh38tohg19_standard'
    path_dhs_tissue_stan = path_mid + '/tissue/DHS/GRCh38tohg19_standard'
    dhs_matrix()

    # H3K4me3
    path_h3k4me3_cell = path_mid + '/cell_line/DHS/reference_map'
    path_h3k4me3_tissue = path_mid + '/tissue/DHS/reference_map'
    # h3k4me3_matrix()

    # H3K27ac
    path_h3k27ac_cell = path_mid + '/cell_line/DHS/cRE_annotation'
    path_h3k27ac_tissue = path_mid + '/tissue/DHS/cRE_annotation'
    # h3k27ac_matrix()

    # GETx expression
    file_expression = \
        path_origin + \
        '/GTEx/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct'
    gtex_expression_matrix()

    time_end = time()
    print(time_end - time_start)
