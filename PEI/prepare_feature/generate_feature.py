#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: generate_feature.py
# @time: 2020/5/6 18:29

from time import time
import pandas as pd
import os
from multiprocessing import Pool
from subprocess import check_output


def sub_get_corr(file_enhancer, df_index, sub_path_out, folders, gene,
                 method='Spearman'):
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
              f"-b {file_enhancer} -sorted "
              f"| cut -f 4,5,6,7,8,9,10,11,12,13,14 > {file_intersect1}")
    os.system(f"bedtools intersect -wb -a {file_intersect1} "
              f"-b {file_gene_pro} -sorted -v > {file_intersect2}")
    len_ref = int(str(check_output(f"wc -l {file_intersect2}",
                                   shell=True).strip()).split(' ')[0][2:])
    if len_ref == 0:
        return
    df_gene_dhs = pd.read_csv(file_intersect2, sep='\t', header=None)
    distance = (df_gene_dhs[1] + df_gene_dhs[2]) / 2 - (int(start_pro) + 2000)
    df_gene_dhs = df_gene_dhs.rename(columns={3: 'dhs_id', 4: 'type_cre'})
    df_gene_dhs = df_gene_dhs.loc[:, ['dhs_id', 'type_cre']]
    df_merge = pd.merge(df_gene_dhs, df_index, on='dhs_id', how='left')
    df_merge['gene'] = [gene for _ in range(df_merge.shape[0])]
    df_merge = df_merge.loc[:, ['gene', 'dhs_id', 'type_cre', 'ref_dhs_id']]
    df_merge['distance'] = distance

    os.remove(file_gene_pro)
    os.remove(file_gene_flank)
    os.remove(file_intersect1)
    os.remove(file_intersect2)

    df_merge_corr = df_merge
    for folder in folders:
        file_corr_gene = os.path.join(
            path_correlation, f"{folder}/gene_corr/{gene}_corr.txt")
        df_gene_corr = pd.read_csv(
            file_corr_gene, sep='\t',
            names=['gene', 'ref_dhs_id', 'Pearson', 'Spearman', 'Kendall'])
        df_gene_corr = df_gene_corr.rename(columns={method: folder})
        df_gene_corr = df_gene_corr.loc[:, ['gene', 'ref_dhs_id', folder]]
        df_merge_corr = pd.merge(df_merge_corr, df_gene_corr,
                                 how='left', on=['gene', 'ref_dhs_id'])

    file_corr = os.path.join(sub_path_out, f"{gene}_corr.txt")
    df_merge_corr.to_csv(file_corr, sep='\t', index=None, header=None)

    return {'gene': gene, 'file_corr': file_corr}


def sub_generate_pair(dict_in):
    file_index = dict_in['file_index']
    path_out = dict_in['path_out']
    file_cre = dict_in['file_cre']
    # term = dict_in['term']
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
                    elif type_cre in set_pro:
                        genes = list_line[6].split(',')
                        for gene in genes:
                            list_out = list_line.copy()
                            list_out[6] = gene
                            w_pro.write('\t'.join(list_out) + '\n')

    df_pro = pd.read_csv(file_cre_pro, sep='\t', header=None)
    genes = df_pro[6].apply(lambda x: x.split('<-')[0]).unique().tolist()

    df_index = pd.read_csv(file_index, sep='\t', header=None,
                           usecols=[0, 2], names=['dhs_id', 'ref_dhs_id'])
    sub_path_out = os.path.join(path_out, 'tmp_gene')
    if not os.path.exists(sub_path_out):
        os.mkdir(sub_path_out)

    folders = os.listdir(path_correlation)
    list_out = []
    for gene in genes:
        dict_out = sub_get_corr(
            file_cre_enh, df_index, sub_path_out, folders, gene, 'Spearman')
        list_out.append(dict_out)

    corr_files = [sub_dict['file_corr'] for sub_dict in list_out
                  if sub_dict is not None]
    # for sub_dict in list_out:
    #     try:
    #         a = sub_dict['file_corr']
    #     except TypeError:
    #         print(sub_dict)
    #         print(sub_dict is not None)
    file_header = os.path.join(path_out, 'header.txt')
    with open(file_header, 'w') as w_header:
        list_header = \
            ['gene', 'dhs_id', 'type_cre', 'ref_dhs_id', 'distance'] + folders
        w_header.write('\t'.join(list_header) + '\n')
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
    file_promoter = '/local/zy/PEI/origin_data/gene/' \
                    'promoters.up2k.protein.gencode.v19.bed'
    file_promoter_uniq = '/local/zy/PEI/origin_data/gene/' \
                         'promoters.up2k.protein.gencode.v19.unique.bed'
    df_promoter = pd.read_csv(file_promoter, sep='\t', header=None)
    df_promoter['idx'] = df_promoter.index
    df_promoter['gene'] = df_promoter[3].apply(
        lambda x: x.split('<-')[0])
    df_promoter = df_promoter.drop_duplicates(subset='gene')
    df_promoter.to_csv(file_promoter_uniq, sep='\t', index=None, header=None)

    dict_gene_pos = {}
    for sub_dict_gene in df_promoter.to_dict('records'):
        dict_gene_pos[sub_dict_gene['gene']] = \
            [sub_dict_gene[0], sub_dict_gene[1], sub_dict_gene[2]]

    path_dhs_cell = \
        '/local/zy/PEI/mid_data/cell_line/DHS/GRCh38tohg19_standard'
    path_dhs_tissue_cluster = \
        '/local/zy/PEI/mid_data/tissue/DHS/GRCh38tohg19_cluster'
    path_dhs_tissue_stan = \
        '/local/zy/PEI/mid_data/tissue/DHS/GRCh38tohg19_standard'

    path_cre_cell = '/local/zy/PEI/mid_data/cell_line/DHS/cRE_annotation'
    path_cre_tissue = '/local/zy/PEI/mid_data/tissue/DHS/cRE_annotation'

    path_feature_cell = '/local/zy/PEI/mid_data/cell_line/model_input'
    path_feature_tissue = '/local/zy/PEI/mid_data/tissue/model_input'

    path_correlation = '/local/zy/PEI/mid_data/database_feature/correlation'
    generate_pairs()

    time_end = time()
    print(time_end - time_start)
