#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: node10_map_to_DHS.py
# @time: 12/6/19 3:50 PM

from time import time
import pandas as pd
import numpy as np
import os
from collections import defaultdict
from multiprocessing import Pool
from functools import partial
from subprocess import check_output


def sub_annotate_promoter(dict_in):
    file_dhs = dict_in['file_dhs']
    path_out = dict_in['path_out']
    path_h3k4me3 = dict_in['path_h3k4me3']
    sub_h3k4me3 = dict_in['sub_h3k4me3']
    loc_promoter = dict_in['loc_promoter']
    accessions = sub_h3k4me3['File accession'].tolist()

    def drop_dup(x):
        if x.shape[0] == 1:
            return x
        else:
            max_overlap = np.max(x.iloc[:, -1])
            row_out = x.loc[x.iloc[:, -1] == max_overlap, :]
            return row_out

    # promoter
    file_promoter = os.path.join(path_out, 'ref_promoter.txt')
    os.system(f"bedtools intersect -a {file_dhs} -b {loc_promoter} -wao "
              f"| cut -f 1,2,3,4,12,15 > {file_promoter}")
    # drop duplicates
    len_ref = int(str(check_output(f"wc -l {file_dhs}",
                                   shell=True).strip()).split(' ')[0][2:])
    len_pro = int(str(check_output(f"wc -l {file_promoter}",
                                   shell=True).strip()).split(' ')[0][2:])
    if len_ref == len_pro:
        file_ref = file_promoter
    else:
        file_promoter_uniq = os.path.join(path_out, 'ref_promoter.uniq.txt')
        file_promoter_sort = os.path.join(path_out, 'ref_promoter.sort.txt')
        df_plus = pd.read_csv(file_promoter, sep='\t', header=None)
        df_0 = df_plus.loc[df_plus[5] == 0, :]
        df_pn = df_plus.loc[df_plus[5] > 0, :]
        df_pn_uniq = df_pn.groupby(3).apply(drop_dup)
        df_pn_uniq = df_pn_uniq.drop_duplicates(subset=3)
        df_uniq = pd.concat([df_0, df_pn_uniq])
        df_uniq.to_csv(file_promoter_uniq, sep='\t', header=None, index=None)
        os.system(f"bedtools sort -i {file_promoter_uniq} > "
                  f"{file_promoter_sort}")
        os.remove(file_promoter)
        os.remove(file_promoter_uniq)
        os.system(f"mv {file_promoter_sort} {file_promoter}")
        file_ref = file_promoter

    for accession in accessions:
        file_accession = os.path.join(path_h3k4me3, accession + '.bed')
        file_plus = os.path.join(path_out, accession + '.plus')
        file_uniq = os.path.join(path_out, accession + '.uniq')
        file_sort = os.path.join(path_out, accession + '.sort')

        # map H3K4me3 to DHS
        col_num = int(
            check_output("head -n 1 " + file_ref + " | awk '{print NF}'",
                         shell=True).strip())
        use_col_list = list(range(1, col_num + 1))
        use_col_list.extend([col_num + 4, col_num + 5, col_num + 8])
        use_col = ','.join([str(num) for num in use_col_list])
        os.system(
            f"bedtools intersect -a {file_ref} -b {file_accession} -wao "
            f"| cut -f {use_col} > {file_plus}")
        # drop duplicates
        len_ref = int(str(check_output(f"wc -l {file_ref}",
                                       shell=True).strip()).split(' ')[0][2:])
        len_pro = int(str(check_output(f"wc -l {file_plus}",
                                       shell=True).strip()).split(' ')[0][2:])
        if len_ref == len_pro:
            file_ref = file_plus
        else:
            df_plus = pd.read_csv(file_plus, sep='\t', header=None)
            df_0 = df_plus.loc[df_plus[df_plus.shape[1] - 1] == 0, :]
            df_pn = df_plus.loc[df_plus[df_plus.shape[1] - 1] > 0, :]
            df_pn_uniq = df_pn.groupby(3).apply(drop_dup)
            df_pn_uniq = df_pn_uniq.drop_duplicates(subset=3)
            df_uniq = pd.concat([df_0, df_pn_uniq])
            df_uniq.to_csv(file_uniq, sep='\t', header=None, index=None)
            os.system(f"sort -k 1,1 -k2,2n {file_uniq} > {file_sort}")

            os.remove(file_ref)
            os.remove(file_plus)
            os.remove(file_uniq)
            file_ref = file_sort

    file_origin = os.path.join(path_out, 'DHS_promoter_H3K4me3.origin')
    os.system(f"mv {file_ref} {file_origin}")

    # adjust p value
    infer_num = np.sum(sub_h3k4me3['Inferred peak number'])
    file_num = sub_h3k4me3.shape[0]
    file_promoter_out = os.path.join(path_out, 'DHS_promoter_H3K4me3.txt')
    os.system(f"Rscript adjust_p_value_H3K4me3.R "
              f"{file_origin} {file_promoter_out} {infer_num} {file_num}")

    return


def annotate_promoter_to_dhs(path_dhs, path_cluster, path_h3k4me3,
                             loc_promoter, ref_dhs, ref_histone, path_out,
                             num_process):
    if os.path.exists(path_out):
        os.system(f"rm -rf {path_out}")
    os.mkdir(path_out)

    df_ref_histone = pd.read_csv(ref_histone, sep='\t')
    df_ref_histone = df_ref_histone.dropna()
    df_ref_dhs = pd.read_csv(ref_dhs, sep='\t')
    df_meta_h3k4me3 = pd.read_csv(
        os.path.join(path_h3k4me3, 'metadata.simple.tsv'), sep='\t'
    )

    organs = list(
        set([organ for organ in df_ref_histone['Biosample organ'].tolist()])
    )

    list_input = []
    for organ in organs:
        str_organ = organ.replace(' ', '_')
        sub_ref_histone = df_ref_histone.loc[
            df_ref_histone['Biosample organ'] == organ, :
        ]
        sub_ref_dhs = df_ref_dhs.loc[df_ref_dhs['Biosample organ'] == organ, :]
        path_organ = os.path.join(path_out, organ.replace(' ', '_'))
        if not os.path.exists(path_organ):
            os.mkdir(path_organ)
        suborgans = list(
            set([suborgan for suborgan in
                 sub_ref_histone['Biosample suborgan'].tolist()])
        )
        for suborgan in suborgans:
            str_suborgan = suborgan.replace(' ', '_')
            suborgan_histone = \
                sub_ref_histone.loc[
                    sub_ref_histone['Biosample suborgan'] == suborgan, :]
            path_suborgan = \
                os.path.join(path_organ, suborgan.replace(' ', '_'))
            if not os.path.exists(path_suborgan):
                os.mkdir(path_suborgan)
            for sub_dict in suborgan_histone.to_dict("records"):
                term = sub_dict['Biosample term name']
                life = sub_dict['Biosample life stage']
                str_term = sub_dict['Biosample term name'].replace(
                    ' ', '_').replace('/', '+').replace("'", '--')
                path_term = os.path.join(path_suborgan, f"{life}_{str_term}")
                if not os.path.exists(path_term):
                    os.mkdir(path_term)
                term_dhs = sub_ref_dhs.loc[
                    (sub_ref_dhs['Biosample life stage'] == life) &
                    (sub_ref_dhs['Biosample term name'] == term), :
                ]
                if term_dhs.shape[0] == 1:
                    file_dhs = os.path.join(
                        path_dhs,
                        f"{str_organ}/{life}/{str_term}/{str_term}.bed"
                    )
                elif term_dhs.shape[0] == 0:
                    file_dhs = os.path.join(
                        path_cluster,
                        f"{str_organ}/{str_suborgan}/{str_suborgan}.bed"
                    )
                sub_h3k4me3 = df_meta_h3k4me3.loc[
                    (df_meta_h3k4me3['Biosample life stage'] == life) &
                    (df_meta_h3k4me3['Biosample term name'] == term), :]
                list_input.append(dict(
                    file_dhs=file_dhs, path_out=path_term,
                    path_h3k4me3=path_h3k4me3, sub_h3k4me3=sub_h3k4me3,
                    loc_promoter=loc_promoter)
                )

    pool = Pool(processes=num_process)
    pool.map(sub_annotate_promoter, list_input)
    pool.close()

    return


def annotate


if __name__ == '__main__':
    time_start = time()
    # multiple processes
    num_cpu = 40

    # annotate DHS by term
    path_dhs_stan = '/local/zy/PEI/data/DHS/GRCh38tohg19_standard'
    path_dhs_cluster = '/local/zy/PEI/data/DHS/GRCh38tohg19_cluster'
    path_h3k4me3_stan = \
        '/local/zy/PEI/data/ENCODE/histone_ChIP-seq/' \
        'GRCh38tohg19/H3K4me3_standard'
    path_h3k27ac_stan = \
        '/local/zy/PEI/data/ENCODE/histone_ChIP-seq/' \
        'GRCh38tohg19/H3K27ac_standard'
    # select data having H3K4me3 and H3K27ac
    file_meta = '/local/zy/PEI/data/ENCODE/histone_ChIP-seq/' \
                'meta.reference.histone.tsv'
    df_h3k4me3 = pd.read_csv(
        os.path.join(path_h3k4me3_stan, 'meta.reference.tsv'), sep='\t'
    )
    df_h3k27ac = pd.read_csv(
        os.path.join(path_h3k27ac_stan, 'meta.reference.tsv'), sep='\t'
    )
    df_intersect = pd.merge(
        df_h3k4me3, df_h3k27ac,
        on=['Biosample organ', 'Biosample life stage', 'Biosample term name']
    )
    # df_intersect.to_csv(file_meta, sep='\t', index=None)

    # promoter reference
    promoter_file_hg19 = \
        '/local/zy/PEI/data/gene/promoters.up2k.protein.gencode.v19.bed'
    meta_suborgan_dhs = '/local/zy/PEI/data/DHS/meta.reference.tsv'
    path_ref_promoter = '/local/zy/PEI/data/DHS/ref_promoter'
    annotate_promoter_to_dhs(
        path_dhs_stan, path_dhs_cluster, path_h3k4me3_stan,
        promoter_file_hg19, meta_suborgan_dhs, file_meta, path_ref_promoter,
        num_cpu
    )

    time_end = time()
    print(time_end - time_start)
