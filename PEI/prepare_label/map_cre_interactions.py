#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: map_cre_interactions.py
# @time: 2020/4/20 16:11

from time import time
import re
import os
import pandas as pd
import numpy as np
from multiprocessing import Pool


def get_pairs(file_in, file_out):
    pattern_gene = re.compile(r'.+<-')
    file_out_tmp = file_out + '.tmp'
    with open(file_out_tmp, 'w') as w_out:
        fmt = "{gene}\t{dhs_id}\t{cre_type}\t{loop_score}\n"
        with open(file_in, 'r') as r_f:
            for line in r_f:
                list_line = line.strip().split('\t')
                loop_id = list_line[-2]
                loop_score = list_line[-1]
                cres1 = list_line[4].split(',')
                cres2 = list_line[10].split(',')
                set_cres1 = set(list_line[4].split(','))
                set_cres2 = set(list_line[10].split(','))
                if ('.' in set_cres1) | ('.' in set_cres2):
                    continue
                dhs_ids1 = list_line[3].split(',')
                dhs_ids2 = list_line[9].split(',')
                promoters1 = set(list_line[5].split(','))
                promoters2 = set(list_line[11].split(','))
                if loop_id[:2] == 'pp':
                    if (('Protein-Promoter' in set_cres1) |
                        ('Protein-Promoter(Enhancer)' in set_cres1)) & \
                            (('Protein-Promoter' in set_cres2) |
                             ('Protein-Promoter(Enhancer)' in set_cres2)):
                        genes1 = [pattern_gene.search(val).group()[:-2]
                                  for val in promoters1 if val != '.']
                        genes2 = [pattern_gene.search(val).group()[:-2]
                                  for val in promoters2 if val != '.']
                        for i, dhs_id in enumerate(dhs_ids1):
                            for gene2 in genes2:
                                dict_pro1 = dict(gene=gene2, dhs_id=dhs_id,
                                                 cre_type=cres1[i],
                                                 loop_score=loop_score)
                                w_out.write(fmt.format(**dict_pro1))
                        for i, dhs_id in enumerate(dhs_ids2):
                            for gene1 in genes1:
                                dict_pro2 = dict(gene=gene1, dhs_id=dhs_id,
                                                 cre_type=cres2[i],
                                                 loop_score=loop_score)
                                w_out.write(fmt.format(**dict_pro2))
                elif loop_id[:2] == 'po':
                    if (('Protein-Promoter' in set_cres1) |
                        ('Protein-Promoter(Enhancer)' in set_cres1)) & \
                            (('Enhancer' in set_cres2) |
                             ('Other-Promoter(Enhancer)' in set_cres2) |
                             ('Protein-Promoter(Enhancer)' in set_cres2)):
                        genes1 = [pattern_gene.search(val).group()[:-2]
                                  for val in promoters1 if val != '.']
                        for i, dhs_id in enumerate(dhs_ids2):
                            for gene1 in genes1:
                                dict_pro2 = dict(gene=gene1, dhs_id=dhs_id,
                                                 cre_type=cres2[i],
                                                 loop_score=loop_score)
                                w_out.write(fmt.format(**dict_pro2))
                else:
                    if (('Protein-Promoter' in set_cres1) |
                        ('Protein-Promoter(Enhancer)' in set_cres1)) & \
                            (('Protein-Promoter' in set_cres2) |
                             ('Protein-Promoter(Enhancer)' in set_cres2)):
                        genes1 = [pattern_gene.search(val).group()[:-2]
                                  for val in promoters1 if val != '.']
                        genes2 = [pattern_gene.search(val).group()[:-2]
                                  for val in promoters2 if val != '.']
                        for i, dhs_id in enumerate(dhs_ids1):
                            for gene2 in genes2:
                                dict_pro1 = dict(gene=gene2, dhs_id=dhs_id,
                                                 cre_type=cres1[i],
                                                 loop_score=loop_score)
                                w_out.write(fmt.format(**dict_pro1))
                        for i, dhs_id in enumerate(dhs_ids2):
                            for gene1 in genes1:
                                dict_pro2 = dict(gene=gene1, dhs_id=dhs_id,
                                                 cre_type=cres2[i],
                                                 loop_score=loop_score)
                                w_out.write(fmt.format(**dict_pro2))
                    elif (('Protein-Promoter' in set_cres1) |
                          ('Protein-Promoter(Enhancer)' in set_cres1)) & \
                            (('Enhancer' in set_cres2) |
                             ('Other-Promoter(Enhancer)' in set_cres2) |
                             ('Protein-Promoter(Enhancer)' in set_cres2)):
                        genes1 = [pattern_gene.search(val).group()[:-2]
                                  for val in promoters1 if val != '.']
                        for i, dhs_id in enumerate(dhs_ids2):
                            for gene1 in genes1:
                                dict_pro2 = dict(gene=gene1, dhs_id=dhs_id,
                                                 cre_type=cres2[i],
                                                 loop_score=loop_score)
                                w_out.write(fmt.format(**dict_pro2))
                    elif (('Enhancer' in set_cres1) |
                          ('Other-Promoter(Enhancer)' in set_cres1) |
                          ('Protein-Promoter(Enhancer)' in set_cres1)) & \
                            (('Protein-Promoter' in set_cres2) |
                             ('Protein-Promoter(Enhancer)' in set_cres2)):
                        genes2 = [pattern_gene.search(val).group()[:-2]
                                  for val in promoters2 if val != '.']
                        for i, dhs_id in enumerate(dhs_ids1):
                            for gene2 in genes2:
                                dict_pro1 = dict(gene=gene2, dhs_id=dhs_id,
                                                 cre_type=cres1[i],
                                                 loop_score=loop_score)
                                w_out.write(fmt.format(**dict_pro1))

    def drop_dup(x):
        if x.shape[0] == 1:
            return x
        else:
            max_overlap = np.max(x.iloc[:, -1])
            row_out = x.loc[x.iloc[:, -1] == max_overlap, :]
            return row_out

    df_pairs = pd.read_csv(file_out_tmp, sep='\t', header=None)
    df_pairs = df_pairs.drop_duplicates()
    df_pairs = df_pairs.rename(columns={0: 'key1', 1: 'key2'})
    df_pairs_out = df_pairs.groupby(['key1', 'key2']).apply(drop_dup)
    df_pairs_out.to_csv(file_out, sep='\t', index=None, header=None)

    # os.system(f"sort {file_out_tmp} | uniq > {file_out}")
    os.remove(file_out_tmp)

    return


def annotate_hic(file_cre, file_uniform, file_out, file_pair):
    file_bin1 = file_out + '.bin1'
    file_bin2 = file_out + '.bin2'
    os.system(f"cut -f 1,2,3,7,8 {file_uniform} > {file_bin1}")
    os.system(f"cut -f 4,5,6,7,8 {file_uniform} > {file_bin2}")
    file_intersect_bin1 = file_out + '.intersect.bin1'
    file_intersect_bin2 = file_out + '.intersect.bin2'
    os.system(f"bedtools intersect -a {file_bin1} -b {file_cre} -loj | "
              f"cut -f 1,2,3,4,5,9,10,12 | "
              f"grep -w 'Promoter\\|Enhancer' | "
              f"bedtools sort -i > {file_intersect_bin1}")
    os.system(f"bedtools intersect -a {file_bin2} -b {file_cre} -loj | "
              f"cut -f 1,2,3,4,5,9,10,12 | "
              f"grep -w 'Promoter\\|Enhancer' | "
              f"bedtools sort -i > {file_intersect_bin2}")

    def _merge_bin(df_in):
        dhs_ids = ','.join(df_in[5].tolist())
        cres = ','.join(df_in[6].tolist())
        promoters = ','.join(df_in[7].tolist())
        dict_out = {0: df_in.iloc[0, 0], 1: df_in.iloc[0, 1],
                    2: df_in.iloc[0, 2], 3: df_in.iloc[0, 3],
                    4: df_in.iloc[0, 4], 5: dhs_ids, 6: cres, 7: promoters}
        df_dict_out = pd.Series(dict_out)

        return df_dict_out

    df_intersect_bin1 = pd.read_csv(file_intersect_bin1, sep='\t', header=None)
    df_intersect_bin1 = df_intersect_bin1.rename(columns={3: 'key'})
    df_merge_bin1 = df_intersect_bin1.groupby('key').apply(_merge_bin)
    df_intersect_bin2 = pd.read_csv(file_intersect_bin2, sep='\t', header=None)
    df_intersect_bin2 = df_intersect_bin2.rename(columns={3: 'key'})
    df_merge_bin2 = df_intersect_bin2.groupby('key').apply(_merge_bin)
    df_bin1 = pd.read_csv(file_bin1, sep='\t', header=None)
    df_bin2 = pd.read_csv(file_bin2, sep='\t', header=None)
    df_res_bin1 = pd.merge(df_bin1, df_merge_bin1,
                           on=[0, 1, 2, 3, 4], how='outer')
    df_res_bin2 = pd.merge(df_bin2, df_merge_bin2,
                           on=[0, 1, 2, 3, 4], how='outer')
    df_out = pd.merge(df_res_bin1, df_res_bin2, on=[3, 4])
    df_out[8] = df_out[3]
    df_out = df_out.drop(3, axis=1)
    df_out[9] = df_out[4]
    df_out = df_out.drop(4, axis=1)
    df_out.to_csv(file_out, sep='\t', header=None, index=None, na_rep='.')

    os.remove(file_bin1)
    os.remove(file_bin2)
    os.remove(file_intersect_bin1)
    os.remove(file_intersect_bin2)

    get_pairs(file_out, file_pair)

    return


def prepare_interaction_data(dict_in):
    term = dict_in['Biosample term name']
    method = dict_in['Method']
    source = dict_in['Source']
    filename = dict_in['Filename']

    path_term = os.path.join(path_label, term)
    file_uniform = os.path.join(
        path_term, f"{method}__{source}__{filename[:-4]}.uniform")

    file_cre = os.path.join(path_ref_cellline, f"{term}/cRE.txt")
    if not os.path.isfile(file_cre):
        return

    file_out = os.path.join(
        path_term, f"{method}__{source}__{filename[:-4]}.interactions.cRE.txt")
    file_pair = os.path.join(
        path_term, f"{method}__{source}__{filename[:-4]}.pairs.gene.cRE.txt")
    annotate_hic(file_cre, file_uniform, file_out, file_pair)

    return


def generate_heatmap_data(file_heatmap):
    df_meta_sort = df_meta.sort_values(
        by=['Name', 'Method', 'Source'], axis=0)

    def calculate_similarity(file1, file2):
        df_1 = pd.read_csv(file1, sep='\t', header=None, usecols=[0, 1])
        df_2 = pd.read_csv(file2, sep='\t', header=None, usecols=[0, 1])
        df_merge = pd.merge(df_1, df_2, how='inner', on=[0, 1])
        score = df_merge.shape[0] / min(df_1.shape[0], df_2.shape[0])

        return score

    dicts_meta = df_meta_sort.to_dict('records')
    list_out = []
    for i in range(len(dicts_meta)):
        dict_in1 = dicts_meta[i]
        term1 = dict_in1['Biosample term name']
        method1 = dict_in1['Method']
        source1 = dict_in1['Source']
        filename1 = dict_in1['Filename']
        name1 = dict_in1['Name']
        path_term1 = os.path.join(path_label, term1)
        file_pair1 = os.path.join(
            path_term1,
            f"{method1}__{source1}__{filename1[:-4]}.pairs.gene.cRE.txt")
        if not os.path.isfile(file_pair1):
            continue
        label1 = f"{name1}_{method1}_{source1}"
        for j in range(len(dicts_meta)):
            dict_in2 = dicts_meta[j]
            term2 = dict_in2['Biosample term name']
            method2 = dict_in2['Method']
            source2 = dict_in2['Source']
            filename2 = dict_in2['Filename']
            name2 = dict_in2['Name']
            path_term2 = os.path.join(path_label, term2)
            file_pair2 = os.path.join(
                path_term2,
                f"{method2}__{source2}__{filename2[:-4]}.pairs.gene.cRE.txt")
            if not os.path.isfile(file_pair2):
                continue
            label2 = f"{name2}_{method2}_{source2}"
            list_out.append(
                {'label1': label1, 'label2': label2,
                 'score': calculate_similarity(file_pair1, file_pair2)})

    df_heatmap = pd.DataFrame(list_out)
    df_heatmap.to_csv(file_heatmap, sep='\t', index=None)

    return


def merge_files(file_num_pair):
    cell_lines = (df_meta['Biosample term name'].unique()).tolist()
    list_out = []
    for cell in cell_lines:
        sub_meta = df_meta.loc[df_meta['Biosample term name'] == cell, :]
        path_term = os.path.join(path_label, cell)
        list_df = []
        for dict_in in sub_meta.to_dict('records'):
            method = dict_in['Method']
            source = dict_in['Source']
            filename = dict_in['Filename']
            file_pair = os.path.join(
                path_term,
                f"{method}__{source}__{filename[:-4]}.pairs.gene.cRE.txt")
            sub_df = pd.read_csv(file_pair, sep='\t', header=None,
                                 usecols=[0, 1, 2])
            list_df.append(sub_df)
        df_cell = pd.concat(list_df)
        df_cell = df_cell.drop_duplicates()
        df_cell.to_csv(os.path.join(path_term, cell + '.txt'),
                       sep='\t', index=None, header=None)
        list_out.append({'Cell line': cell, 'Num of pairs': df_cell.shape[0]})

    df_heatmap = pd.DataFrame(list_out)
    df_heatmap.to_csv(file_num_pair, sep='\t', index=None)

    return


if __name__ == '__main__':
    time_start = time()
    path_ref_cellline = '/local/zy/PEI/mid_data/cell_line/DHS/cRE_annotation'
    path_label = \
        '/local/zy/PEI/mid_data/training_label/label_interactions_cutoff3'

    flie_meta = os.path.join(path_label, 'meta_label.txt')
    df_meta = pd.read_csv(flie_meta, sep='\t')

    pool = Pool(processes=40)
    pool.map(prepare_interaction_data, df_meta.to_dict('records'))
    pool.close()

    generate_heatmap_data(os.path.join(path_label, 'heatmap.txt'))

    merge_files(os.path.join(path_label, 'Num_pairs.txt'))

    time_end = time()
    print(time_end - time_start)