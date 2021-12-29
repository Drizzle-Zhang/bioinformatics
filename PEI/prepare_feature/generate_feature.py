#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: generate_feature.py
# @time: 2020/5/6 18:29

from time import time
import pandas as pd
import numpy as np
import os
from multiprocessing import Pool
from subprocess import check_output
import glob


def gtex2encode(file_mat_promoter, meta_transfer_col):
    df_mat_pro = pd.read_csv(file_mat_promoter, sep='\t', index_col=0)
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

    return df_mat_pro


def sub_get_corr(file_enhancer, df_index, sub_df_gene, sub_path_out, folders,
                 gene, header_gene, file_header, method='Spearman'):
    chrom_pro = dict_gene_pos[gene][0]
    start_pro = dict_gene_pos[gene][1]
    end_pro = dict_gene_pos[gene][2]
    file_gene_pro = os.path.join(sub_path_out, f"{gene}_promoter.bed")
    with open(file_gene_pro, 'w') as w_pro:
        w_pro.write(f"{chrom_pro}\t{start_pro}\t{end_pro}")
    file_gene_flank = os.path.join(sub_path_out, f"{gene}_flank.bed")
    flank = 2_000_000 - 2000
    with open(file_gene_flank, 'w') as w_flank:
        w_flank.write(
            f"{chrom_pro}\t{max(int(start_pro) - flank, 0)}\t"
            f"{int(end_pro) + flank}")
    file_intersect1 = os.path.join(sub_path_out, f'{gene}.intersect1')
    file_intersect2 = os.path.join(sub_path_out, f'{gene}.intersect2')
    os.system(f"bedtools intersect -wb -a {file_gene_flank} "
              f"-b {file_enhancer} "
              f"| cut -f 4,5,6,7,8,9,10,11,12,13,14 > {file_intersect1}")
    os.system(f"bedtools intersect -wb -a {file_intersect1} "
              f"-b {file_gene_pro} -v > {file_intersect2}")
    len_ref = int(str(check_output(f"wc -l {file_intersect2}",
                                   shell=True).strip()).split(' ')[0][2:])
    if len_ref == 0:
        os.remove(file_gene_pro)
        os.remove(file_gene_flank)
        os.remove(file_intersect1)
        os.remove(file_intersect2)
        return
    df_gene_enh = pd.read_csv(file_intersect2, sep='\t', header=None)
    distance = (df_gene_enh[1] + df_gene_enh[2]) / 2 - (int(start_pro) + 2000)
    df_gene_dhs = df_gene_enh.rename(columns={3: 'dhs_id', 4: 'type_cre'})
    df_gene_dhs = df_gene_dhs.loc[:, ['dhs_id', 'type_cre']]
    df_merge = pd.merge(df_gene_dhs, df_index, on='dhs_id', how='left')
    df_merge['gene'] = [gene for _ in range(df_merge.shape[0])]
    df_merge = df_merge.loc[:, ['gene', 'dhs_id', 'type_cre', 'ref_dhs_id']]
    df_merge['distance'] = distance

    df_gene_enh = df_gene_enh.loc[:, [3, 5, 9]]
    df_gene_enh.columns = ['dhs_id', 'DHS', 'H3K27ac']
    df_gene_enh['gene'] = [gene for _ in range(df_gene_enh.shape[0])]

    os.remove(file_gene_pro)
    os.remove(file_gene_flank)
    os.remove(file_intersect1)
    os.remove(file_intersect2)

    df_merge_corr = df_merge.copy()
    for folder in folders:
        factor1, factor2 = folder.split('_')
        if (factor1 == 'expression') & (factor2 == 'DHS'):
            if 'expression_DHS' in sub_df_gene.index:
                vec_factor1 = sub_df_gene.loc['expression_DHS']
            else:
                vec_factor1 = None
        elif (factor1 == 'expression') & (factor2 == 'H3K27ac'):
            if 'expression_H3K27ac' in sub_df_gene.index:
                vec_factor1 = sub_df_gene.loc['expression_H3K27ac']
            else:
                vec_factor1 = None
        else:
            vec_factor1 = sub_df_gene.loc[factor1]
        vec_factor2 = df_gene_enh.loc[:, factor2]
        if vec_factor1:
            signal = np.sqrt(vec_factor1 * vec_factor2)
            signal.name = 'signal_' + folder
            df_signal = pd.concat(
                [df_gene_enh.loc[:, ['gene', 'dhs_id']], signal], axis=1)
            df_merge_corr = pd.merge(df_merge_corr, df_signal,
                                     how='left', on=['gene', 'dhs_id'])

        file_corr_gene = os.path.join(
            path_correlation, f"{folder}/gene_corr/{gene}_corr.txt")
        # df_gene_corr = pd.read_csv(
        #     file_corr_gene, sep='\t',
        #     names=['gene', 'ref_dhs_id', f'{folder}_Pearson',
        #            f'{folder}_Spearman', f'{folder}_Kendall'])
        df_gene_corr = pd.read_csv(
            file_corr_gene, sep='\t',
            names=['gene', 'ref_dhs_id', 'Pearson', 'Spearman', 'Kendall'])
        df_gene_corr = df_gene_corr.rename(columns={method: folder})
        df_gene_corr = df_gene_corr.loc[:, ['gene', 'ref_dhs_id', folder]]
        df_merge_corr = pd.merge(df_merge_corr, df_gene_corr,
                                 how='left', on=['gene', 'ref_dhs_id'])

    file_corr = os.path.join(sub_path_out, f"{gene}_corr.txt")
    df_merge_corr.to_csv(file_corr, sep='\t', index=None, header=None)

    if gene == header_gene:
        with open(file_header, 'w') as w_header:
            list_header = df_merge_corr.columns
            w_header.write('\t'.join(list_header) + '\n')

    return {'gene': gene, 'file_corr': file_corr}


def sub_generate_pair(dict_in):
    file_index = dict_in['file_index']
    path_out = dict_in['path_out']
    file_cre = dict_in['file_cre']
    term = dict_in['term']
    # str_term = dict_in['str_term']

    file_cre_pro = os.path.join(path_out, 'cRE_promoter.txt')
    file_cre_enh = os.path.join(path_out, 'cRE_enhancer.txt')

    set_pro = {'Protein-Promoter', 'Protein-Promoter(Enhancer)'}
    set_enh = {'Other-Promoter(Enhancer)', 'Protein-Promoter(Enhancer)',
               'Enhancer'}
    with open(file_cre, 'r') as r_cre:
        with open(file_cre_pro, 'w') as w_pro:
            with open(file_cre_enh, 'w') as w_enh:
                for line in r_cre:
                    list_line = line.strip().split('\t')
                    type_cre = list_line[4]
                    if type_cre in set_enh:
                        w_enh.write(line)
                    if type_cre in set_pro:
                        genes = list_line[6].split(',')
                        for gene in genes:
                            list_out = list_line.copy()
                            list_out[6] = gene
                            w_pro.write('\t'.join(list_out) + '\n')

    df_pro = pd.read_csv(file_cre_pro, sep='\t', header=None)
    genes = df_pro[6].apply(lambda x: x.split('<-')[0]).unique().tolist()
    genes = set(genes).intersection(set_genes_gtex)

    df_pro['gene'] = df_pro[6].apply(lambda x: x.split('<-')[0])
    df_genes = df_pro.loc[:, ['gene', 5, 7]]

    def max_scores(df_in):
        if df_in.shape[0] == 1:
            df_out = df_in.loc[:, [5, 7]]
        else:
            array_max = np.array(np.max(df_in.loc[:, [5, 7]])).reshape((1, 2))
            df_out = pd.DataFrame(array_max, columns=[5, 7])

        return df_out

    df_genes = df_genes.groupby('gene').apply(max_scores)
    df_genes.index = [index[0] for index in df_genes.index]
    df_genes.columns = ['DHS', 'H3K4me3']

    # add expression
    if term in df_expression_dhs.columns:
        df_exp_dhs = df_expression_dhs.loc[:, term]
        df_exp_dhs.name = 'expression_DHS'
        df_genes = pd.concat([df_genes, df_exp_dhs], axis=1, sort=False)
        df_genes = df_genes.dropna()
    if term in df_expression_h3k27ac.columns:
        df_exp_h3k27ac = df_expression_h3k27ac.loc[:, term]
        df_exp_h3k27ac.name = 'expression_H3K27ac'
        df_genes = pd.concat([df_genes, df_exp_h3k27ac], axis=1, sort=False)
        df_genes = df_genes.dropna()

    # index associating dhs_id and ref_dhs_id
    df_index = pd.read_csv(file_index, sep='\t', header=None,
                           usecols=[0, 2], names=['dhs_id', 'ref_dhs_id'])
    sub_path_out = os.path.join(path_out, 'tmp_gene')
    if not os.path.exists(sub_path_out):
        os.mkdir(sub_path_out)

    folders = os.listdir(path_correlation)
    # folders = ['expression_DHS']
    list_out = []
    header_gene = list(genes)[20]
    file_header = os.path.join(path_out, 'header.txt')
    for gene in list(genes):
        # gene = list(genes)[20]
        sub_df_gene = df_genes.loc[gene, :]
        dict_out = sub_get_corr(
            file_cre_enh, df_index, sub_df_gene, sub_path_out, folders,
            gene, header_gene, file_header, 'Spearman')
        # break
        list_out.append(dict_out)

    corr_files = [sub_dict['file_corr'] for sub_dict in list_out
                  if sub_dict is not None]
    # corr_files = glob.glob(os.path.join(sub_path_out, '*_corr.txt'))

    # for sub_dict in list_out:
    #     try:
    #         a = sub_dict['file_corr']
    #     except TypeError:
    #         print(sub_dict)
    #         print(sub_dict is not None)
    # file_header = os.path.join(path_out, 'header.txt')
    # with open(file_header, 'w') as w_header:
    #     list_header = \
    #         ['gene', 'dhs_id', 'type_cre', 'ref_dhs_id', 'distance'] + folders
    #     # all_corr_folders = []
    #     # for folder in folders:
    #     #     for corr in ['Pearson', 'Spearman', 'Kendall']:
    #     #         all_corr_folders.append(f"{folder}|{corr}")
    #     # for folder in folders:
    #     #     for corr in ['Spearman']:
    #     #         all_corr_folders.append(f"{folder}|{corr}")
    #     # list_header = \
    #     #     ['gene', 'dhs_id', 'type_cre', 'ref_dhs_id', 'distance'] + all_corr_folders
    #     w_header.write('\t'.join(list_header) + '\n')
    if not os.path.exists(file_header):
        return
    list_cat = [file_header]
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


def generate_pairs():
    # generate P-E pairs for each cell lines or tissues

    # cell line
    df_meta_cell = pd.read_csv(
        os.path.join(path_cre_cell, 'meta.reference.tsv'), sep='\t')
    list_input = []
    for term in (df_meta_cell['Biosample term name'].unique()).tolist():
        str_term = term.replace(' ', '_').replace('/', '+')
        path_term_dhs = os.path.join(path_dhs_cell, str_term)
        file_index = os.path.join(path_term_dhs, 'index.txt')
        if not os.path.exists(file_index):
            continue
        path_term_cre = os.path.join(path_cre_cell, str_term)
        file_cre = os.path.join(path_term_cre, 'cRE.txt')
        if not os.path.exists(file_cre):
            continue
        path_out = os.path.join(path_feature_cell, str_term)
        if not os.path.exists(path_out):
            os.mkdir(path_out)
        list_input.append({'file_index': file_index, 'term': term,
                           'str_term': str_term, 'file_cre': file_cre,
                           'path_out': path_out})

    pool = Pool(num_cpu)
    pool.map(sub_generate_pair, list_input)
    pool.close()

    return


if __name__ == '__main__':
    time_start = time()
    num_cpu = 40
    # path_root = '/local/zy/PEI'
    # path_origin = path_root + '/origin_data'
    # path_mid = path_root + '/mid_data'
    path_root = '/lustre/tianlab/zhangyu/PEI'
    path_origin = path_root + '/origin_data'
    path_mid = path_root + '/mid_data_correct'

    file_promoter = \
        path_origin + '/gene/promoters.up2k.protein.gencode.v19.bed'
    file_promoter_uniq = \
        path_origin + '/gene/promoters.up2k.protein.gencode.v19.unique.bed'
    df_promoter = pd.read_csv(file_promoter, sep='\t', header=None)
    df_promoter['idx'] = df_promoter.index
    df_promoter['gene'] = df_promoter[3].apply(
        lambda x: x.split('<-')[0])
    df_promoter = df_promoter.drop_duplicates(subset='gene')
    df_promoter.to_csv(file_promoter_uniq, sep='\t', index=None, header=None)

    # use genes in GTEx
    file_mat_exp_gtex = \
        path_mid + '/database_feature/matrix/GTEx_expression_matrix.txt'
    set_genes_gtex = set(pd.read_csv(
        file_mat_exp_gtex, sep='\t', usecols=[0]).iloc[:, 0].tolist())
    # genes not in GTEx
    # file_mat_h3k4me3 = \
    #     path_mid + '/database_feature/matrix/H3K4me3_matrix.txt'
    # set_genes_encode = set(pd.read_csv(
    #     file_mat_h3k4me3, sep='\t', usecols=[0]).iloc[:, 0].tolist())
    # set_not_gtex = set_genes_encode.difference(set_genes_gtex)

    dict_gene_pos = {}
    for sub_dict_gene in df_promoter.to_dict('records'):
        dict_gene_pos[sub_dict_gene['gene']] = \
            [sub_dict_gene[0], sub_dict_gene[1], sub_dict_gene[2]]

    # expression matrix
    # DHS
    meta_dhs = path_origin + '/meta_file/meta_GTEx_DHS.txt'
    df_expression_dhs = gtex2encode(file_mat_exp_gtex, meta_dhs)

    # H3K27ac
    meta_h3k27ac = path_origin + '/meta_file/meta_GTEx_H3K27ac.txt'
    df_expression_h3k27ac = gtex2encode(file_mat_exp_gtex, meta_h3k27ac)

    # cell line
    path_dhs_cell = \
        path_mid + '/cell_line/DHS/GRCh38tohg19_standard'
    path_dhs_tissue_cluster = \
        path_mid + '/tissue/DHS/GRCh38tohg19_cluster'
    path_dhs_tissue_stan = \
        path_mid + '/tissue/DHS/GRCh38tohg19_standard'

    path_cre_cell = path_mid + '/cell_line/DHS/cRE_annotation'
    path_cre_tissue = path_mid + '/tissue/DHS/cRE_annotation'

    path_feature_cell = path_mid + '/cell_line/model_input'
    path_feature_tissue = path_mid + '/tissue/model_input'

    path_correlation = path_mid + '/database_feature/correlation'
    generate_pairs()

    time_end = time()
    print(time_end - time_start)
