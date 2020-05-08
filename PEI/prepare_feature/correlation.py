#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: correlation.py
# @time: 2020/5/5 23:28

from time import time
import pandas as pd
import os
import numpy as np
from multiprocessing import Pool, Process
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


def sub_corr(df_mat_pro, df_mat_dhs, path_out, gene):
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

    vec_gene = df_mat_pro.loc[gene, :]
    mat_dhs = df_mat_dhs.loc[list_dhs, :].T
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

    return {'gene': gene, 'file_corr': file_corr}


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


def correlation(mat_promoter, mat_dhs, path_out):
    df_mat_pro = pd.read_csv(mat_promoter, sep='\t', index_col=0)
    df_mat_dhs = pd.read_csv(mat_dhs, sep='\t', index_col=0)
    col_pro = df_mat_pro.columns
    col_dhs = df_mat_dhs.columns
    col_overlap = set(col_pro).intersection(col_dhs)
    df_mat_pro = df_mat_pro.loc[:, col_overlap]
    df_mat_dhs = df_mat_dhs.loc[:, col_overlap]
    genes = df_mat_pro.index
    path_gene = os.path.join(path_out, 'gene_corr')
    if not os.path.exists(path_gene):
        os.mkdir(path_gene)

    # pool = Pool(num_cpu)
    # func_get_dhs = partial(get_dhs, path_gene)
    # list_dict = pool.map(func_get_dhs, genes)
    # pool.close()
    #
    # # save tmp result
    # df_tmp = pd.DataFrame(list_dict)
    # df_tmp.to_csv(os.path.join(path_out, 'tmp.txt'), sep='\t')
    df_tmp = pd.read_csv(os.path.join(path_out, 'tmp.txt'), sep='\t')
    list_dict = df_tmp.to_dict('records')

    subprocesses = []
    for idx, dict_in in enumerate(list_dict):
        gene = dict_in['gene']
        list_dhs = dict_in['list_dhs']
        # tmp
        list_dhs = list_dhs[2:-2].split("', '")
        vec_gene = df_mat_pro.loc[gene, :]
        mat_dhs = df_mat_dhs.loc[list_dhs, :].T
        dict_in['vec_gene'] = vec_gene
        dict_in['mat_dhs'] = mat_dhs
        if len(subprocesses) % num_cpu == 0:
            for process in subprocesses:
                process.join()
            subprocesses = []
        process = Process(target=calculate_corr, args=(path_gene, dict_in))
        process.start()
        print(idx, gene)
        subprocesses.append(process)

    for process in subprocesses:
        process.join()

    corr_files = \
        [os.path.join(path_gene, f"{gene}_corr.txt") for gene in genes]
    list_cat = []
    for idx in range(len(corr_files)//200):
        if idx < len(corr_files)//200:
            sub_cat_in = ' '.join(corr_files[200*idx:200*(idx+1)])
        else:
            sub_cat_in = ' '.join(corr_files[200*idx:len(corr_files)])
        sub_cat_out = os.path.join(path_out, f'sub_{i}.txt')
        os.system(f"cat {sub_cat_in} > {sub_cat_out}")
        list_cat.append(sub_cat_out)
    cat_in = ' '.join(list_cat)
    cat_out = os.path.join(path_out, 'correlation.txt')
    os.system(f"cat {cat_in} > {cat_out}")
    os.system(f"rm {' '.join(list_cat)}")

    return


if __name__ == '__main__':
    time_start = time()
    num_cpu = 40
    file_all_index = \
        '/local/zy/PEI/mid_data/database_feature/DHS_index/all_index.txt'
    file_dhs_promoter = \
        '/local/zy/PEI/mid_data/database_feature/DHS_index/promoter_index.txt'
    file_promoter = '/local/zy/PEI/origin_data/gene/' \
                    'promoters.up2k.protein.gencode.v19.bed'
    dict_gene_pos = generate_promoter_dict()

    path_matrix = '/local/zy/PEI/mid_data/database_feature/matrix'
    matrix_gene = ['DHS', 'H3K4me3']
    files_gene = ['DHS_matrix.promoter.txt', 'H3K27ac_matrix.txt']
    matrix_dhs = ['DHS', 'H3K27ac']
    files_dhs = ['DHS_matrix.txt', 'H3K4me3_matrix.txt']

    path_correlation = '/local/zy/PEI/mid_data/database_feature/correlation'

    for i, name_gene in enumerate(matrix_gene):
        for j, name_dhs in enumerate(matrix_dhs):
            file_gene = os.path.join(path_matrix, files_gene[i])
            file_dhs = os.path.join(path_matrix, files_dhs[j])
            sub_path_out = os.path.join(
                path_correlation, f"{name_gene}_{name_dhs}")
            if not os.path.exists(sub_path_out):
                os.mkdir(sub_path_out)
            correlation(file_gene, file_dhs, sub_path_out)
            if j == 0:
                break
        if i == 0:
            break

    time_end = time()
    print(time_end - time_start)
