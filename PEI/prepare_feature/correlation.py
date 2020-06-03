#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: correlation.py
# @time: 2020/5/5 23:28

from time import time
import pandas as pd
import os
import numpy as np
from multiprocessing import Pool
from functools import partial
from scipy.stats import spearmanr, pearsonr, kendalltau


def generate_promoter_dict():
    df_pro_idx = pd.read_csv(file_dhs_promoter, sep='\t')
    set_gene_pos = set(df_pro_idx['gene_pos'].tolist())

    dict_gene = {}
    with open(file_promoter, 'r') as r_pro:
        for line in r_pro:
            tmp_line = line.strip().split('\t')
            gene_pos = tmp_line[3]
            if gene_pos in set_gene_pos:
                gene = gene_pos.split('<-')[0]
                dict_gene[gene] = [tmp_line[0], tmp_line[1], tmp_line[2]]

    return dict_gene


def get_dhs(path_out, gene):
    chrom_pro = dict_gene_pos[gene][0]
    start_pro = dict_gene_pos[gene][1]
    end_pro = dict_gene_pos[gene][2]
    file_gene_pro = os.path.join(path_out, f"{gene}_promoter.bed")
    with open(file_gene_pro, 'w') as w_pro:
        w_pro.write(f"{chrom_pro}\t{start_pro}\t{end_pro}")
    file_gene_flank = os.path.join(path_out, f"{gene}_flank.bed")
    flank = 2_000_000 - 2000
    with open(file_gene_flank, 'w') as w_flank:
        w_flank.write(
            f"{chrom_pro}\t{max(int(start_pro) - flank, 0)}\t"
            f"{int(end_pro) + flank}")
    file_intersect1 = os.path.join(path_out, f'{gene}.intersect1')
    file_intersect2 = os.path.join(path_out, f'{gene}.intersect2')
    os.system(f"bedtools intersect -wb -a {file_gene_flank} "
              f"-b {file_all_index} -sorted "
              f"| cut -f 4,5,6,7,8 > {file_intersect1}")
    os.system(f"bedtools intersect -wb -a {file_intersect1} "
              f"-b {file_gene_pro} -sorted -v > {file_intersect2}")
    df_gene_dhs = pd.read_csv(file_intersect2, sep='\t', header=None)
    list_dhs = df_gene_dhs[3].tolist()

    os.remove(file_gene_pro)
    os.remove(file_gene_flank)
    os.remove(file_intersect1)
    os.remove(file_intersect2)

    return {'gene': gene, 'list_dhs': list_dhs}


def calculate_corr(path_out, dict_in):
    gene = dict_in['gene']
    vec_gene = dict_in['vec_gene']
    mat_dhs = dict_in['mat_dhs']
    # pearson correlation
    dict_corr = {}
    for col in mat_dhs.columns:
        vec_dhs = mat_dhs[col]
        rho, _ = pearsonr(vec_gene, vec_dhs)
        dict_corr[col] = rho
    corr_pearson = pd.Series(dict_corr)
    corr_pearson = corr_pearson.fillna(0)
    # spearman correlation
    dict_corr = {}
    for col in mat_dhs.columns:
        vec_dhs = mat_dhs[col]
        rho, _ = spearmanr(vec_gene, vec_dhs)
        dict_corr[col] = rho
    corr_spearman = pd.Series(dict_corr)
    corr_spearman = corr_spearman.fillna(0)
    # kendall correlation
    dict_corr = {}
    for col in mat_dhs.columns:
        vec_dhs = mat_dhs[col]
        rho, _ = kendalltau(vec_gene, vec_dhs)
        dict_corr[col] = rho
    corr_kendall = pd.Series(dict_corr)
    corr_kendall = corr_kendall.fillna(0)

    df_corr = pd.concat(
        [corr_pearson, corr_spearman, corr_kendall], axis=1, sort=False)
    df_corr.columns = ['pearson', 'spearman', 'kendall']
    df_corr['dhs_id'] = df_corr.index
    df_corr['gene'] = [gene for _ in range(df_corr.shape[0])]
    df_corr = df_corr.loc[
              :, ['gene', 'dhs_id', 'pearson', 'spearman', 'kendall']]
    file_corr = os.path.join(path_out, f"{gene}_corr.txt")
    df_corr.to_csv(file_corr, sep='\t', index=None, header=None)

    return


def discrete(score):
    if score <= 0:
        score_out = 0
    elif (score > 0) & (score <= 0.1):
        score_out = 1
    elif (score > 0.1) & (score <= 0.2):
        score_out = 2
    elif (score > 0.2) & (score <= 0.3):
        score_out = 3
    elif (score > 0.3) & (score <= 0.4):
        score_out = 4
    elif (score > 0.4) & (score <= 0.5):
        score_out = 5
    elif (score > 0.5) & (score <= 0.6):
        score_out = 6
    elif (score > 0.6) & (score <= 0.7):
        score_out = 7
    elif (score > 0.7) & (score <= 0.8):
        score_out = 8
    elif (score > 0.8) & (score <= 0.9):
        score_out = 9
    elif (score > 0.9) & (score <= 1):
        score_out = 10
    else:
        score_out = np.nan

    return score_out


def correlation(name_gene_in, name_dhs_in, file_mat_promoter,
                file_mat_dhs, path_out):
    df_mat_pro = pd.read_csv(file_mat_promoter, sep='\t', index_col=0)
    if (name_gene_in == 'expression') & (name_dhs_in == 'DHS'):
        meta_transfer_col = \
            path_origin + '/meta_file/meta_GTEx_DHS.txt'
    elif (name_gene_in == 'expression') & (name_dhs_in == 'H3K27ac'):
        meta_transfer_col = \
            path_origin + '/meta_file/meta_GTEx_H3K27ac.txt'
    else:
        meta_transfer_col = None
    if meta_transfer_col:
        df_meta_trans = pd.read_csv(meta_transfer_col, sep='\t')
        df_meta_trans = df_meta_trans.loc[
            (df_meta_trans['Biosample life stage'].apply(
                lambda x: isinstance(x, str))) |
            (df_meta_trans['Biosample term name'].apply(
                lambda x: isinstance(x, str))), :]
        col_encode = []
        for idx in df_meta_trans.index:
            life = df_meta_trans.loc[idx, 'Biosample life stage']
            organ = df_meta_trans.loc[idx, 'Biosample organ']
            term = df_meta_trans.loc[idx, 'Biosample term name']
            suborgan = df_meta_trans.loc[idx, 'Biosample suborgan']
            if isinstance(life, float):
                col_encode.append(term)
            else:
                life_organ = f"{life}_{organ}"
                if isinstance(term, float):
                    col_encode.append(f"{life_organ} | {suborgan}")
                else:
                    col_encode.append(f"{life_organ} | {term}")
        df_meta_trans['ENCODE tissue name'] = col_encode
        cols_gtex = df_meta_trans['GTEx tissue name'].tolist()
        df_mat_pro = df_mat_pro.loc[:, cols_gtex]
        df_mat_pro.columns = col_encode

    df_mat_dhs = pd.read_csv(file_mat_dhs, sep='\t', index_col=0)
    col_pro = df_mat_pro.columns
    col_dhs = df_mat_dhs.columns
    col_overlap = set(col_pro).intersection(col_dhs)
    # record tissue and cell overlaps
    file_corr_variables = os.path.join(path_out, 'variables_corr.txt')
    with open(file_corr_variables, 'w') as w_var:
        for var in col_overlap:
            w_var.write(var + '\n')
    df_mat_pro = df_mat_pro.loc[:, col_overlap]
    df_mat_dhs = df_mat_dhs.loc[:, col_overlap]
    genes = df_mat_pro.index
    path_gene = os.path.join(path_out, 'gene_corr')
    if not os.path.exists(path_gene):
        os.mkdir(path_gene)

    pool = Pool(num_cpu)
    func_get_dhs = partial(get_dhs, path_gene)
    list_dict = pool.map(func_get_dhs, genes)
    pool.close()
    #
    # # save tmp result
    # df_tmp = pd.DataFrame(list_dict)
    # df_tmp.to_csv(os.path.join(path_out, 'tmp.txt'), sep='\t')
    # df_tmp = pd.read_csv(os.path.join(path_out, 'tmp.txt'), sep='\t')
    # list_dict = df_tmp.to_dict('records')

    list_input = []
    for idx, dict_in in enumerate(list_dict):
        sub_dict = dict_in.copy()
        gene = sub_dict['gene']
        list_dhs = sub_dict['list_dhs']
        # tmp
        # list_dhs = list_dhs[2:-2].split("', '")
        vec_gene = df_mat_pro.loc[gene, :]
        mat_dhs = df_mat_dhs.loc[list_dhs, :].T
        sub_dict['vec_gene'] = vec_gene
        sub_dict['mat_dhs'] = mat_dhs
        list_input.append(sub_dict)
        print(idx, gene)
        if len(list_input) % 200 == 0:
            pool = Pool(num_cpu)
            func_calc = partial(calculate_corr, path_gene)
            pool.map(func_calc, list_input)
            pool.close()
            list_input = []

    pool = Pool(num_cpu)
    func_calc = partial(calculate_corr, path_gene)
    pool.map(func_calc, list_input)
    pool.close()

    corr_files = \
        [os.path.join(path_gene, f"{gene}_corr.txt") for gene in genes]
    list_cat = []
    for idx in range(len(corr_files)//200):
        if idx < len(corr_files)//200:
            sub_cat_in = ' '.join(corr_files[200*idx:200*(idx+1)])
        else:
            sub_cat_in = ' '.join(corr_files[200*idx:len(corr_files)])
        sub_cat_out = os.path.join(path_out, f'sub_{idx}.txt')
        if os.path.isfile(sub_cat_out):
            os.remove(sub_cat_out)
        os.system(f"cat {sub_cat_in} > {sub_cat_out}")
        list_cat.append(sub_cat_out)
    cat_in = ' '.join(list_cat)
    cat_out = os.path.join(path_out, 'correlation.txt')
    if os.path.isfile(cat_out):
        os.remove(cat_out)
    os.system(f"cat {cat_in} > {cat_out}")
    os.system(f"rm {' '.join(list_cat)}")

    return


if __name__ == '__main__':
    time_start = time()
    num_cpu = 40
    # path_root = '/local/zy/PEI'
    path_root = '/lustre/tianlab/zhangyu/PEI'
    path_origin = path_root + '/origin_data'
    path_mid = path_root + '/mid_data_correct'

    file_all_index = \
        path_mid + '/database_feature/DHS_index/all_index.txt'
    file_dhs_promoter = \
        path_mid + '/database_feature/DHS_index/promoter_index.txt'
    file_promoter = \
        path_origin + '/gene/promoters.up2k.protein.gencode.v19.bed'
    dict_gene_pos = generate_promoter_dict()

    path_matrix = path_mid + '/database_feature/matrix'
    # matrix_gene = ['DHS', 'H3K4me3', 'expression']
    matrix_gene = ['DHS']
    # files_gene = ['DHS_matrix.promoter.txt', 'H3K4me3_matrix.txt',
    # 'GTEx_expression_matrix.txt']
    files_gene = ['DHS_matrix.promoter.txt']
    # matrix_dhs = ['DHS', 'H3K27ac']
    # files_dhs = ['DHS_matrix.txt', 'H3K27ac_matrix.txt']
    matrix_dhs = ['DHS']
    files_dhs = ['DHS_matrix.txt']

    path_correlation = path_mid + '/database_feature/correlation'

    for i, name_gene in enumerate(matrix_gene):
        for j, name_dhs in enumerate(matrix_dhs):
            file_gene = os.path.join(path_matrix, files_gene[i])
            file_dhs = os.path.join(path_matrix, files_dhs[j])
            sub_path_out = os.path.join(
                path_correlation, f"{name_gene}_{name_dhs}")
            if not os.path.exists(sub_path_out):
                os.mkdir(sub_path_out)
            # if f"{name_gene}_{name_dhs}" == 'DHS_DHS':
            #     continue
            correlation(name_gene, name_dhs, file_gene, file_dhs, sub_path_out)
        #     if j == 0:
        #         break
        # if i == 0:
        #     break

    time_end = time()
    print(time_end - time_start)
