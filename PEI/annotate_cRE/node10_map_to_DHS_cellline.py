#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: node10_map_to_DHS.py
# @time: 12/6/19 3:50 PM

from time import time
import pandas as pd
import numpy as np
import os
from multiprocessing import Pool
from functools import partial
from subprocess import check_output
root_path = '/local/zy/my_git/bioinformatics/PEI/annotate_cRE'


def drop_dup(x):
    if x.shape[0] == 1:
        return x
    else:
        max_overlap = np.max(x.iloc[:, -1])
        row_out = x.loc[x.iloc[:, -1] == max_overlap, :]
        return row_out


def sub_annotate_promoter(dict_in):
    file_dhs = dict_in['file_dhs']
    path_out = dict_in['path_out']
    path_h3k4me3 = dict_in['path_h3k4me3']
    sub_h3k4me3 = dict_in['sub_h3k4me3']
    loc_promoter = dict_in['loc_promoter']
    accessions = sub_h3k4me3['File accession'].tolist()

    # promoter
    file_promoter = os.path.join(path_out, 'ref_promoter.txt')
    os.system(f"bedtools intersect -a {file_dhs} -b {loc_promoter} -wao "
              f"| cut -f 1,2,3,4,5,12,15 > {file_promoter}")
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
        df_0 = df_plus.loc[df_plus[df_plus.shape[1] - 1] == 0, :]
        df_pn = df_plus.loc[df_plus[df_plus.shape[1] - 1] > 0, :]
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
        use_col_list.extend([col_num + 4, col_num + 5, col_num + 6,
                             col_num + 8])
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
    os.system(f"Rscript {os.path.join(root_path, 'adjust_p_value_H3K4me3.R')} "
              f"{file_origin} {file_promoter_out} {infer_num} {file_num}")

    return


def annotate_promoter_to_dhs(path_dhs, path_h3k4me3,
                             loc_promoter, path_out, num_process):
    if os.path.exists(path_out):
        os.system(f"rm -rf {path_out}")
    os.mkdir(path_out)

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
    len_ref = int(str(check_output(f"wc -l {file_ref}",
                                   shell=True).strip()).split(' ')[0][2:])

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
            assert len_ref == df_plus.shape[0]
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
    os.system(f"Rscript {os.path.join(root_path, 'adjust_p_value_H3K27ac.R')} "
              f"{file_origin} {file_out} {infer_num} {file_num}")

    return


def integrate_ctcf(path_ctcf, dict_in):
    file_ref = dict_in['file_ref']
    path_out = dict_in['path_out']
    sub_ctcf = dict_in['sub_ctcf']
    accessions = sub_ctcf['File accession'].tolist()
    len_ref = int(str(check_output(f"wc -l {file_ref}",
                                   shell=True).strip()).split(' ')[0][2:])

    file_ref_ori = file_ref
    for accession in accessions:
        file_accession = os.path.join(path_ctcf, accession + '.bed')
        file_plus = os.path.join(path_out, accession + '.plus')
        file_uniq = os.path.join(path_out, accession + '.uniq')
        file_sort = os.path.join(path_out, accession + '.sort')

        # map H3K4me3 to DHS
        col_num = int(
            check_output("head -n 1 " + file_ref + " | awk '{print NF}'",
                         shell=True).strip())
        use_col_list = list(range(1, col_num + 1))
        use_col_list.extend([col_num + 4, col_num + 5, col_num + 6,
                             col_num + 8])
        use_col = ','.join([str(num) for num in use_col_list])
        os.system(
            f"bedtools intersect -a {file_ref} -b {file_accession} -wao "
            f"| cut -f {use_col} > {file_plus}")
        # drop duplicates
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

            if file_ref != file_ref_ori:
                os.remove(file_ref)
            os.remove(file_plus)
            os.remove(file_uniq)
            file_ref = file_sort

            # check error
            len_sort = int(str(check_output(
                f"wc -l {file_sort}", shell=True).strip()).split(' ')[0][2:])
            try:
                assert len_ref == len_sort
            except AssertionError:
                print(dict_in)
                return

    file_origin = os.path.join(
        path_out, 'DHS_promoter_H3K4me3_H3K27ac_CTCF.origin')
    os.system(f"mv {file_ref} {file_origin}")

    df_origin = pd.read_csv(file_origin, sep='\t', header=None)
    file_num = sub_ctcf.shape[0]
    cols = [11 + 4*i for i in range(file_num)]
    df_out = df_origin.iloc[:, :10]
    df_scores = df_origin.loc[:, cols]
    df_scores = df_scores.applymap(lambda x: float(x) if x != '.' else -10000)
    df_out[10] = np.max(df_scores, axis=1)
    file_out = os.path.join(
        path_out, 'DHS_promoter_H3K4me3_H3K27ac_CTCF.txt')
    df_out.to_csv(file_out, sep='\t', header=None, index=None)

    return


def sub_annotate_cre(dict_in):
    path = dict_in['path_out']
    file_in = os.path.join(path, 'DHS_promoter_H3K4me3_H3K27ac_CTCF.txt')
    if not os.path.exists(file_in):
        return
    file_out = os.path.join(path, 'cRE.txt')
    with open(file_in, 'r') as r_f:
        with open(file_out, 'w') as w_f:
            fmt_dhs = "{chrom}\t{start}\t{end}\t{dhs_id}\t{cre}\t" \
                      "{dhs_score}\t{promoter_id}\t{score_h3k4me3}\t" \
                      "{p_h3k4me3}\t{score_h3k27ac}\t{p_h3k27ac}\t" \
                      "{score_ctcf}\n"
            for line in r_f:
                list_line = line.strip().split('\t')
                chrom = list_line[0]
                start = list_line[1]
                end = list_line[2]
                dhs_id = list_line[3]
                dhs_score = float(list_line[4])
                promoter_id = list_line[5]
                score_h3k4me3 = float(list_line[6])
                p_h3k4me3 = float(list_line[7])
                score_h3k27ac = float(list_line[8])
                p_h3k27ac = float(list_line[9])
                score_ctcf = float(list_line[10])
                if (promoter_id != '.') & (p_h3k4me3 != 0) & \
                        (p_h3k27ac != 0):
                    cre = 'Protein-Promoter(Enhancer)'
                elif (promoter_id == '.') & (p_h3k4me3 != 0) & \
                        (p_h3k27ac != 0):
                    cre = 'Other-Promoter(Enhancer)'
                elif (promoter_id != '.') & (p_h3k4me3 != 0) & \
                        (p_h3k27ac == 0):
                    cre = 'Protein-Promoter'
                elif (p_h3k4me3 == 0) & (p_h3k27ac != 0):
                    cre = 'Enhancer'
                elif (p_h3k27ac == 0) & (score_ctcf != -10000):
                    cre = 'Insulator'
                else:
                    cre = '.'
                w_f.write(fmt_dhs.format(**locals()))

    # file_no_exon = os.path.join(path, 'cRE_no_exon.txt')
    # file_no_exon_1 = file_no_exon + '.tmp1'
    # file_no_exon_2 = file_no_exon + '.tmp2'
    # file_no_exon_3 = file_no_exon + '.tmp3'
    # os.system(f"bedtools intersect -a {file_out} "
    #           f"-b {protein_exon} -v > {file_no_exon_1}")
    # # promoter_file_hg19
    # os.system(f"bedtools intersect -a {file_out} "
    #           f"-b {promoter_file_hg19} > {file_no_exon_2}")
    # if os.path.exists(file_no_exon_3):
    #     os.remove(file_no_exon_3)
    # else:
    #     os.system(f"cat {file_no_exon_1} {file_no_exon_2} > "
    #               f"{file_no_exon_3}")
    # os.system(f"bedtools sort -i {file_no_exon_3} | uniq > {file_no_exon}")
    # os.remove(file_no_exon_1)
    # os.remove(file_no_exon_2)
    # os.remove(file_no_exon_3)

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

    pool = Pool(processes=num_process)
    pool.map(sub_annotate_cre, list_input)
    pool.close()

    return


if __name__ == '__main__':
    time_start = time()
    # multiple processes
    num_cpu = 40

    # annotate DHS by term
    path_dhs_stan = \
        '/local/zy/PEI/mid_data/cell_line/DHS/GRCh38tohg19_standard'
    path_h3k4me3_stan = \
        '/local/zy/PEI/mid_data/cell_line/ENCODE/histone_ChIP-seq/' \
        'H3K4me3_standard'
    path_h3k27ac_stan = \
        '/local/zy/PEI/mid_data/cell_line/ENCODE/histone_ChIP-seq/' \
        'H3K27ac_standard'
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
        '/local/zy/PEI/origin_data/gene/' \
        'promoters.up2k.protein.gencode.v19.merge.bed'
    meta_suborgan_dhs = \
        '/local/zy/PEI/mid_data/cell_line/DHS/meta.reference.tsv'
    path_ref_promoter = '/local/zy/PEI/mid_data/cell_line/DHS/reference_map'
    annotate_promoter_to_dhs(
        path_dhs_stan, path_h3k4me3_stan,
        promoter_file_hg19, path_ref_promoter, num_cpu
    )

    # map H3K27ac to reference
    path_ctcf_stan = \
        '/local/zy/PEI/mid_data/cell_line/ENCODE/TF_ChIP-seq/CTCF_standard'
    protein_exon = '/local/zy/PEI/origin_data/gene/' \
                   'exon.protein.gencode.v19.bed'
    path_map_h3k27ac = '/local/zy/PEI/mid_data/cell_line/DHS/map_H3K27ac'
    path_combine_h3k27ac = \
        '/local/zy/PEI/mid_data/cell_line/DHS/cRE_annotation'
    annotate_cre(path_ref_promoter, path_h3k27ac_stan, path_ctcf_stan,
                 path_map_h3k27ac, path_combine_h3k27ac, num_cpu)

    time_end = time()
    print(time_end - time_start)
