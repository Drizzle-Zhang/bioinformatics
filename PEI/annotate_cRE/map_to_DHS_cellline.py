#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: map_to_DHS.py
# @time: 12/6/19 3:50 PM

from time import time
import pandas as pd
import numpy as np
import os
from multiprocessing import Pool
from functools import partial
from subprocess import check_output
import sys
# sys.path.append('/local/zy/my_git/bioinformatics/PEI/annotate_cRE')
sys.path.append(
    '/lustre/tianlab/zhangyu/my_git/bioinformatics/PEI/annotate_cRE')
from map_to_DHS import \
    drop_dup, sub_annotate_promoter, integrate_ctcf, sub_annotate_cre, \
    sub_annotate_cre_ctcf
# root_path = '/local/zy/my_git/bioinformatics/PEI/annotate_cRE'
root_path = '/lustre/tianlab/zhangyu/my_git/bioinformatics/PEI/annotate_cRE'


def annotate_promoter_to_dhs(path_dhs, path_h3k4me3,
                             loc_promoter, path_out, num_process):
    if os.path.exists(path_out):
        os.system(f"rm -rf {path_out}")
    os.mkdir(path_out)
    os.system(f"cp {os.path.join(path_h3k4me3, 'metadata.simple.tsv')} "
              f"{os.path.join(path_out, 'metadata.simple.tsv')}")
    os.system(f"cp {os.path.join(path_h3k4me3, 'meta.reference.tsv')} "
              f"{os.path.join(path_out, 'meta.reference.tsv')}")

    df_meta_h3k4me3 = pd.read_csv(
        os.path.join(path_h3k4me3, 'metadata.simple.tsv'), sep='\t'
    )

    list_input = []
    for term in list_terms:
        str_term = term.replace(' ', '_').replace('/', '+')
        path_term = os.path.join(path_out, str_term)
        if not os.path.exists(path_term):
            os.mkdir(path_term)
        file_dhs = os.path.join(path_dhs, f"{str_term}/{str_term}.bed")
        suborgan_h3k4me3 = df_meta_h3k4me3.loc[
            df_meta_h3k4me3['Biosample term name'] == term, :]
        list_input.append(dict(
            file_dhs=file_dhs, path_out=path_term,
            path_h3k4me3=path_h3k4me3, sub_h3k4me3=suborgan_h3k4me3,
            loc_promoter=loc_promoter)
        )

    pool = Pool(processes=num_process)
    pool.map(sub_annotate_promoter, list_input)
    pool.close()

    return


def map_h3k27ac(path_ref, path_h3k27ac, path_out, dict_in):
    accession_ids = dict_in['File accession']
    file_in = os.path.join(path_h3k27ac, accession_ids + '.bed')
    str_term = dict_in['Biosample term name'].replace(
                    ' ', '_').replace('/', '+')
    file_ref = os.path.join(
        path_ref, f"{str_term}/DHS_promoter_H3K4me3.txt"
    )

    # map H3K4me3 to DHS
    file_plus = os.path.join(path_out, accession_ids + '.plus')
    file_uniq = os.path.join(path_out, accession_ids + '.uniq')
    file_sort = os.path.join(path_out, accession_ids + '.sort')

    os.system(
        f"bedtools intersect -a {file_ref} -b {file_in} -wao "
        f"| cut -f 1,2,3,4,5,6,7,8,12,13,14,16 > {file_plus}")
    # drop duplicates
    df_plus = pd.read_csv(file_plus, sep='\t', header=None)
    df_0 = df_plus.loc[df_plus[df_plus.shape[1] - 1] == 0, :]
    df_pn = df_plus.loc[df_plus[df_plus.shape[1] - 1] > 0, :]
    df_pn_uniq = df_pn.groupby(3).apply(drop_dup)
    df_pn_uniq = df_pn_uniq.drop_duplicates(subset=3)
    df_uniq = pd.concat([df_0, df_pn_uniq])
    df_uniq.to_csv(file_uniq, sep='\t', header=None, index=None)
    os.system(f"sort -k 1,1 -k2,2n {file_uniq} > {file_sort}")

    os.remove(file_plus)
    os.remove(file_uniq)

    file_out = os.path.join(path_out, accession_ids + '.bed')
    os.system(f"mv {file_sort} {file_out}")

    return


def integrate_h3k27ac(path_h3k27ac_map, dict_in):
    file_ref = dict_in['file_ref']
    path_out = dict_in['path_out']
    sub_h3k27ac = dict_in['sub_h3k27ac']
    accessions = sub_h3k27ac['File accession'].tolist()
    ref_col_num = int(
        check_output("head -n 1 " + file_ref + " | awk '{print NF}'",
                     shell=True).strip())

    df_ref = pd.read_csv(file_ref, sep='\t', header=None)
    for accession in accessions:
        file_accession = os.path.join(path_h3k27ac_map, accession + '.bed')
        df_accession = pd.read_csv(
            file_accession, sep='\t', header=None,
            usecols=[0, 1, 2, 3, 8, 9, 10, 11])
        df_plus = pd.merge(df_ref, df_accession, on=[0, 1, 2, 3], how='left')
        df_ref = df_plus

        # check error
        try:
            assert ref_col_num == df_plus.shape[0]
        except AssertionError:
            print(dict_in)
            return

    file_origin = os.path.join(path_out, 'DHS_promoter_H3K4me3_H3K27ac.origin')
    df_ref.to_csv(file_origin, sep='\t', header=None, index=None)

    # adjust p value
    infer_num = np.sum(sub_h3k27ac['Inferred peak number'])
    file_num = sub_h3k27ac.shape[0]
    file_out = os.path.join(
        path_out, 'DHS_promoter_H3K4me3_H3K27ac.txt')
    os.system(f"Rscript {os.path.join(root_path, 'adjust_p_value_histone.R')} "
              f"{file_origin} {file_out} {infer_num} {file_num} {ref_col_num}")

    return


def annotate_cre(path_ref, path_h3k27ac, path_ctcf,
                 path_out_h3k27ac, path_cre, num_process):
    if os.path.exists(path_out_h3k27ac):
        os.system(f"rm -rf {path_out_h3k27ac}")
    os.mkdir(path_out_h3k27ac)
    os.system(f"cp {os.path.join(path_h3k27ac, 'metadata.simple.tsv')} "
              f"{os.path.join(path_out_h3k27ac, 'metadata.simple.tsv')}")
    os.system(f"cp {os.path.join(path_h3k27ac, 'meta.reference.tsv')} "
              f"{os.path.join(path_out_h3k27ac, 'meta.reference.tsv')}")

    df_meta_h3k27ac = pd.read_csv(
        os.path.join(path_h3k27ac, 'metadata.simple.tsv'), sep='\t'
    )
    df_merge = df_meta_h3k27ac.loc[
        df_meta_h3k27ac['Biosample term name'].apply(
            lambda x: x in set(list_terms)), :]

    # map H3K27ac to sample
    pool = Pool(processes=num_process)
    func_map = partial(map_h3k27ac, path_ref, path_h3k27ac, path_out_h3k27ac)
    pool.map(func_map, df_merge.to_dict('records'))
    pool.close()

    # integrate H3K27ac by term and suborgan
    if os.path.exists(path_cre):
        os.system(f"rm -rf {path_cre}")
    os.mkdir(path_cre)
    os.system(f"cp {os.path.join(path_h3k27ac, 'metadata.simple.tsv')} "
              f"{os.path.join(path_cre, 'metadata.simple.tsv')}")
    os.system(f"cp {os.path.join(path_h3k27ac, 'meta.reference.tsv')} "
              f"{os.path.join(path_cre, 'meta.reference.tsv')}")

    list_input = []
    for term in list_terms:
        str_term = term.replace(' ', '_').replace('/', '+')
        path_term = os.path.join(path_cre, str_term)
        if not os.path.exists(path_term):
            os.mkdir(path_term)
        file_ref = os.path.join(
            path_ref, f"{str_term}/DHS_promoter_H3K4me3.txt"
        )
        sub_h3k27ac = df_meta_h3k27ac.loc[
                      df_meta_h3k27ac['Biosample term name'] == term, :]
        list_input.append(dict(
            file_ref=file_ref, path_out=path_term,
            sub_h3k27ac=sub_h3k27ac)
        )

    pool = Pool(processes=num_process)
    func_integrate = partial(integrate_h3k27ac, path_out_h3k27ac)
    pool.map(func_integrate, list_input)
    pool.close()
    print('Annotation of H3K27ac is completed!')

    # integrate CTCF
    os.system(f"cp {os.path.join(path_ctcf, 'metadata.simple.tsv')} "
              f"{os.path.join(path_cre, 'metadata.simple.CTCF.tsv')}")
    os.system(f"cp {os.path.join(path_ctcf, 'meta.reference.tsv')} "
              f"{os.path.join(path_cre, 'meta.reference.CTCF.tsv')}")
    df_meta_ctcf = pd.read_csv(
        os.path.join(path_ctcf, 'metadata.simple.tsv'), sep='\t')
    list_terms_ctcf = \
        set(df_meta_ctcf['Biosample term name'].tolist()).intersection(
            list_terms)

    list_input_ctcf = []
    for term in list_terms_ctcf:
        str_term = term.replace(' ', '_').replace('/', '+')
        path_term = os.path.join(path_cre, str_term)
        if not os.path.exists(path_term):
            os.mkdir(path_term)
        file_ref = os.path.join(
            path_cre, f"{str_term}/DHS_promoter_H3K4me3_H3K27ac.txt"
        )
        sub_ctcf = df_meta_ctcf.loc[
                   df_meta_ctcf['Biosample term name'] == term, :]
        list_input_ctcf.append(dict(
            file_ref=file_ref, path_out=path_term,
            sub_ctcf=sub_ctcf)
        )

    pool = Pool(processes=num_process)
    func_integrate = partial(integrate_ctcf, path_ctcf)
    pool.map(func_integrate, list_input_ctcf)
    pool.close()
    print('Annotation of CTCF is completed!')

    pool = Pool(processes=num_process)
    pool.map(sub_annotate_cre, list_input)
    pool.close()

    pool = Pool(processes=num_process)
    pool.map(sub_annotate_cre_ctcf, list_input_ctcf)
    pool.close()

    return


if __name__ == '__main__':
    time_start = time()
    # multiple processes
    num_cpu = 40
    # path_root = '/local/zy/PEI'
    path_root = '/lustre/tianlab/zhangyu/PEI'
    path_origin = path_root + '/origin_data'
    path_mid = path_root + '/mid_data_correct'

    # annotate DHS by term
    path_dhs_stan = \
        path_mid + '/cell_line/DHS/GRCh38tohg19_standard'
    path_h3k4me3_stan = \
        path_mid + '/cell_line/ENCODE/histone_ChIP-seq/H3K4me3_standard'
    path_h3k27ac_stan = \
        path_mid + '/cell_line/ENCODE/histone_ChIP-seq/H3K27ac_standard'
    # select data having H3K4me3 and H3K27ac
    df_dhs = pd.read_csv(
        os.path.join(path_dhs_stan, 'metadata.simple.tsv'), sep='\t'
    )
    df_h3k4me3 = pd.read_csv(
        os.path.join(path_h3k4me3_stan, 'metadata.simple.tsv'), sep='\t'
    )
    df_h3k27ac = pd.read_csv(
        os.path.join(path_h3k27ac_stan, 'metadata.simple.tsv'), sep='\t'
    )
    list_terms = list(
        (set(df_dhs['Biosample term name'].tolist()).intersection(
            set(df_h3k4me3['Biosample term name'].tolist())
        )).intersection(
            set(df_h3k27ac['Biosample term name'].tolist())
        ))

    # promoter reference
    promoter_file_hg19 = \
        path_origin + '/gene/promoters.up2k.protein.gencode.v19.bed'
    meta_suborgan_dhs = \
        path_mid + '/cell_line/DHS/meta.reference.tsv'
    path_ref_promoter = path_mid + '/cell_line/DHS/reference_map'
    annotate_promoter_to_dhs(
        path_dhs_stan, path_h3k4me3_stan,
        promoter_file_hg19, path_ref_promoter, num_cpu
    )
    print('Annotation of promoters and H3K4me3 is completed!')

    # map H3K27ac to reference
    path_ctcf_stan = \
        path_mid + '/cell_line/ENCODE/TF_ChIP-seq/CTCF_standard'
    protein_exon = path_origin + '/gene/exon.protein.gencode.v19.bed'
    path_map_h3k27ac = path_mid + '/cell_line/DHS/map_H3K27ac'
    path_combine_h3k27ac = \
        path_mid + '/cell_line/DHS/cRE_annotation'
    annotate_cre(path_ref_promoter, path_h3k27ac_stan, path_ctcf_stan,
                 path_map_h3k27ac, path_combine_h3k27ac, num_cpu)
    print('Annotation of H3K27ac and CTCF is completed!')

    time_end = time()
    print(time_end - time_start)
