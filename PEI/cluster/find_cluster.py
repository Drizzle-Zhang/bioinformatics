#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: find_cluster.py
# @time: 7/20/20 11:09 AM

from time import time
import pandas as pd
import os
import numpy as np
import json
from sklearn.decomposition import NMF, FactorAnalysis, LatentDirichletAllocation
from sklearn.linear_model import LassoCV, LinearRegression


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


def get_dhs_per_gene(path_out, gene):
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

    dict_out = {'gene': gene, 'list_dhs': list_dhs}
    str_json = json.dumps(dict_out)
    file_json = os.path.join(path_out, f"{gene}.json")
    with open(file_json, 'w') as w_json:
        w_json.write(str_json)

    return file_json


def decomposition():
    # file_mat_dhs = \
    #     '/lustre/tianlab/zhangyu/PEI/mid_data_correct/database_feature/' \
    #     'matrix/DHS_matrix.txt'
    name_gene_in = 'expression'
    name_dhs_in = 'H3K27ac'
    file_mat_promoter = \
        '/lustre/tianlab/zhangyu/PEI/mid_data_correct/database_feature/' \
        'matrix/GTEx_expression_matrix.txt'
    file_mat_dhs = \
        '/lustre/tianlab/zhangyu/PEI/mid_data_correct/database_feature/' \
        'matrix/H3K27ac_matrix.txt'
    tmp_path = '/lustre/tianlab/zhangyu/PEI/mid_data_correct/database_feature/decomposition/DHS'

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
    file_corr_variables = os.path.join(tmp_path, 'variables_corr.txt')
    with open(file_corr_variables, 'w') as w_var:
        for var in col_overlap:
            w_var.write(var + '\n')
    df_mat_pro = df_mat_pro.loc[:, col_overlap]
    df_mat_dhs = df_mat_dhs.loc[:, col_overlap]

    # GM12878 enhancer
    file_cre = '/lustre/tianlab/zhangyu/PEI/mid_data_correct/cell_line/' \
               'model_input/GM12878/GM12878_feature_label.txt'
    df_feature_label = pd.read_csv(file_cre, sep='\t', usecols=[0, 1, 2, 3])

    # path_json = '/lustre/tianlab/zhangyu/PEI/mid_data_correct/database_feature/correlation/DHS_DHS/gene_corr'
    path_json = '/lustre/tianlab/zhangyu/PEI/mid_data_correct/database_feature/correlation/DHS_H3K27ac/gene_corr'
    gene = 'SSH1'
    file_json = os.path.join(path_json, f"{gene}.json")
    with open(file_json, 'r') as r_json:
        str_json = r_json.readline()
        sub_dict = json.loads(str_json)
    gene = sub_dict['gene']
    list_dhs = sub_dict['list_dhs']

    # select distal enhancer
    df_promoter_dhs = pd.read_csv(file_dhs_promoter, sep='\t')
    set_pro_dhs = set(df_promoter_dhs['dhs_id'].tolist())
    sel_dhs = list(set(list_dhs).difference(set_pro_dhs))
    # sub_dhs = df_feature_label.loc[
    #     df_feature_label['gene'] == gene, 'ref_dhs_id'].tolist()
    # sel_dhs = list(set(list_dhs).intersection(set(sub_dhs)))

    mat_dhs = df_mat_dhs.loc[sel_dhs, :].T

    mat_dhs = mat_dhs.loc[:, np.sum(mat_dhs, axis=0) != 0]
    samples = mat_dhs.index
    list_dhs = mat_dhs.columns

    model_nmf = NMF(n_components=10, init='random', random_state=720)
    w_mat = model_nmf.fit_transform(mat_dhs)
    h_mat = model_nmf.components_
    recon_error = model_nmf.reconstruction_err_

    # model_fa = FactorAnalysis(n_components=10, random_state=720)
    # w_mat = model_fa.fit_transform(mat_dhs)
    # h_mat = model_fa.components_
    # score_fa = model_fa.score(mat_dhs)
    #
    # model_lda = LatentDirichletAllocation(n_components=4, random_state=720)
    # w_mat = model_lda.fit_transform(mat_dhs)
    # h_mat = model_lda.components_

    df_w_mat = pd.DataFrame(w_mat, index=samples)
    df_h_mat = pd.DataFrame(h_mat.T, index=list_dhs)
    #
    # file_w = os.path.join(tmp_path, f"{gene}_w.txt")
    # file_h = os.path.join(tmp_path, f"{gene}_h.txt")
    # df_w_mat.to_csv(file_w, sep='\t')
    # df_h_mat.to_csv(file_h, sep='\t')

    exp_gene = df_mat_pro.loc[gene, :]
    model_lasso = LassoCV(cv=5, random_state=720)
    reg = model_lasso.fit(df_w_mat, exp_gene)
    r_square = reg.score(df_w_mat, exp_gene)
    coef = np.abs(model_lasso.coef_)

    coef[w_mat[-1, :] == 0] = 0
    score_dhs = np.dot(df_h_mat, coef)
    df_score = pd.Series(score_dhs, index=list_dhs)

    sub_dhs = df_feature_label.loc[
        df_feature_label['gene'] == gene, 'ref_dhs_id'].tolist()
    sel_dhs = list(set(list_dhs).intersection(set(sub_dhs)))
    df_out = df_score.loc[sel_dhs]

    # model_linear = LinearRegression()
    # reg = model_linear.fit(df_w_mat.iloc[:, :10], exp_gene)
    # r_square = reg.score(df_w_mat.iloc[:, :10], exp_gene)
    # coef = model_linear.coef_


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

    time_end = time()
    print(time_end - time_start)
