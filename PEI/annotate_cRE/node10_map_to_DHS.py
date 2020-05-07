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

    def drop_dup_promoter(x):
        if x.shape[0] == 1:
            return x
        else:
            row_out = x.iloc[0, 0:5]
            row_out[5] = ','.join([x.iloc[i, 5] for i in range(x.shape[0])])
            row_out[6] = \
                ','.join([str(x.iloc[i, 6]) for i in range(x.shape[0])])
            dict_out = row_out.to_dict()
            row_out = pd.DataFrame(dict_out, index=[row_out[3]])
            return row_out

    if len_ref == len_pro:
        file_ref = file_promoter
    else:
        file_promoter_uniq = os.path.join(path_out, 'ref_promoter.uniq.txt')
        file_promoter_sort = os.path.join(path_out, 'ref_promoter.sort.txt')
        df_plus = pd.read_csv(file_promoter, sep='\t', header=None,
                              dtype={6: 'str'})
        df_0 = df_plus.loc[df_plus[df_plus.shape[1] - 1] == '0', :]
        df_pn = (df_plus.loc[df_plus[df_plus.shape[1] - 1] != '0', :]).copy()
        df_pn['key'] = df_pn[3]
        df_pn_uniq = df_pn.groupby('key').apply(drop_dup_promoter)
        df_pn_uniq = df_pn_uniq.drop('key', axis=1)
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
            df_plus = pd.read_csv(file_plus, sep='\t', header=None,
                                  dtype={6: 'str'})
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


def annotate_promoter_to_dhs(path_cluster, path_dhs, path_h3k4me3,
                             loc_promoter, ref_histone, ref_dhs, path_out,
                             num_process):
    if os.path.exists(path_out):
        os.system(f"rm -rf {path_out}")
    os.mkdir(path_out)
    os.system(f"cp {os.path.join(path_h3k4me3, 'metadata.simple.tsv')} "
              f"{os.path.join(path_out, 'metadata.simple.tsv')}")

    df_ref_histone = pd.read_csv(ref_histone, sep='\t')
    df_ref_histone = df_ref_histone.fillna('single')
    df_ref_histone['Biosample life_organ'] = df_ref_histone.apply(
        lambda x: x['Biosample life stage'] + '_' + x['Biosample organ'],
        axis=1
    )
    df_ref_histone.to_csv(
        os.path.join(path_out, 'meta.reference.histone.tsv'), sep='\t')
    df_ref_dhs = pd.read_csv(ref_dhs, sep='\t')
    df_ref_dhs['Biosample life_organ'] = df_ref_dhs.apply(
        lambda x: x['Biosample life stage'] + '_' + x['Biosample organ'],
        axis=1
    )
    df_meta_h3k4me3 = pd.read_csv(
        os.path.join(path_h3k4me3, 'metadata.simple.tsv'), sep='\t'
    )

    life_organs = list(set(df_ref_histone['Biosample life_organ'].tolist()))

    list_input = []
    list_ref = []
    for life_organ in life_organs:
        str_life_organ = life_organ.replace(' ', '_')
        sub_ref_histone = df_ref_histone.loc[
            df_ref_histone['Biosample life_organ'] == life_organ, :
        ]
        sub_ref_dhs = df_ref_dhs.loc[
            df_ref_dhs['Biosample life_organ'] == life_organ, :]
        if sub_ref_dhs.shape[0] == 0:
            continue
        path_life_organ = os.path.join(path_out, str_life_organ)
        if not os.path.exists(path_life_organ):
            os.mkdir(path_life_organ)
        suborgans = list(set(sub_ref_histone['Biosample suborgan'].tolist()))
        for suborgan in suborgans:
            suborgan_histone = sub_ref_histone.loc[
                sub_ref_histone['Biosample suborgan'] == suborgan, :]
            suborgan_dhs = sub_ref_dhs.loc[
                sub_ref_dhs['Biosample suborgan'] == suborgan, :]
            str_suborgan = suborgan.replace(' ', '_')
            terms = list(set(suborgan_histone['Biosample term name'].tolist()))
            if suborgan != 'single':
                path_suborgan = os.path.join(path_life_organ, str_suborgan)
                if suborgan_dhs.shape[0] == 0:
                    file_dhs_suborgan = os.path.join(
                        path_cluster,
                        f"{str_life_organ}/{str_life_organ}.bed"
                    )
                else:
                    file_dhs_suborgan = os.path.join(
                        path_cluster,
                        f"{str_life_organ}/{str_suborgan}/{str_suborgan}.bed"
                    )
                suborgan_h3k4me3 = pd.merge(
                    suborgan_histone, df_meta_h3k4me3,
                    on=['Biosample life stage', 'Biosample term name']
                )
                suborgan_h3k4me3 = \
                    suborgan_h3k4me3.drop_duplicates('File accession')
                if suborgan_h3k4me3.shape[0] == 0:
                    continue
                list_input.append(dict(
                    file_dhs=file_dhs_suborgan, path_out=path_suborgan,
                    path_h3k4me3=path_h3k4me3, sub_h3k4me3=suborgan_h3k4me3,
                    loc_promoter=loc_promoter)
                )
                list_ref.append(
                    {'Biosample life_organ': life_organ,
                     'Biosample suborgan': suborgan,
                     'Biosample term name': 'empty',
                     'file_ref_dhs': file_dhs_suborgan, 'Level': 'suborgan'})
            else:
                path_suborgan = path_life_organ
            if not os.path.exists(path_suborgan):
                os.mkdir(path_suborgan)
            for term in terms:
                term_histone = suborgan_histone.loc[
                    suborgan_histone['Biosample term name'] == term, :]
                term_dhs = sub_ref_dhs.loc[
                    sub_ref_dhs['Biosample term name'] == term, :]
                str_term = term.replace(
                    ' ', '_').replace('/', '+').replace("'", '--')
                path_term = os.path.join(path_suborgan, str_term)
                if not os.path.exists(path_term):
                    os.mkdir(path_term)
                if term_dhs.shape[0] > 0:
                    file_dhs_term = os.path.join(
                        path_dhs,
                        f"{str_life_organ}/{str_term}/{str_term}.bed"
                    )
                else:
                    if suborgan != 'single':
                        file_dhs_term = file_dhs_suborgan
                    else:
                        file_dhs_term = os.path.join(
                            path_cluster,
                            f"{str_life_organ}/{str_life_organ}.bed"
                        )
                term_h3k4me3 = pd.merge(
                    term_histone, df_meta_h3k4me3,
                    on=['Biosample life stage', 'Biosample term name']
                )
                term_h3k4me3 = \
                    term_h3k4me3.drop_duplicates('File accession')
                if term_h3k4me3.shape[0] == 0:
                    continue
                list_input.append(dict(
                    file_dhs=file_dhs_term, path_out=path_term,
                    path_h3k4me3=path_h3k4me3, sub_h3k4me3=term_h3k4me3,
                    loc_promoter=loc_promoter)
                )
                list_ref.append(
                    {'Biosample life_organ': life_organ,
                     'Biosample suborgan': suborgan,
                     'Biosample term name': term,
                     'file_ref_dhs': file_dhs_term, 'Level': 'term'})

    df_ref_h3k4me3 = pd.DataFrame(list_ref)
    df_ref_h3k4me3.to_csv(
        os.path.join(path_out, 'meta.reference.tsv'), sep='\t', index=None)

    pool = Pool(processes=num_process)
    pool.map(sub_annotate_promoter, list_input)
    pool.close()

    return


def map_h3k27ac(path_ref, path_h3k27ac, path_out, dict_in):
    accession_ids = dict_in['File accession']
    file_in = os.path.join(path_h3k27ac, accession_ids + '.bed')
    str_life_organ = dict_in['Biosample life_organ'].replace(' ', '_')
    str_suborgan = dict_in['Biosample suborgan'].replace(' ', '_')
    str_term = dict_in['Biosample term name'].replace(
                    ' ', '_').replace('/', '+').replace("'", '--')
    if str_suborgan == 'single':
        file_ref = os.path.join(
            path_ref,
            f"{str_life_organ}/{str_term}/DHS_promoter_H3K4me3.txt"
        )
    else:
        file_ref = os.path.join(
            path_ref, f"{str_life_organ}/{str_suborgan}/"
                      f"{str_term}/DHS_promoter_H3K4me3.txt"
        )
    if not os.path.exists(file_ref):
        return

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


def integrate_h3k27ac(path_h3k27ac, dict_in):
    file_ref = dict_in['file_ref']
    path_out = dict_in['path_out']
    sub_h3k27ac = dict_in['sub_h3k27ac']
    accessions = sub_h3k27ac['File accession'].tolist()
    len_ref = int(str(check_output(f"wc -l {file_ref}",
                                   shell=True).strip()).split(' ')[0][2:])

    file_ref_ori = file_ref
    for accession in accessions:
        file_accession = os.path.join(path_h3k27ac, accession + '.bed')
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

    file_origin = os.path.join(path_out, 'DHS_promoter_H3K4me3_H3K27ac.origin')
    os.system(f"mv {file_ref} {file_origin}")

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
    file_in = os.path.join(path, 'DHS_promoter_H3K4me3_H3K27ac.txt')
    if not os.path.exists(file_in):
        return
    file_out = os.path.join(path, 'cRE.txt')
    with open(file_in, 'r') as r_f:
        with open(file_out, 'w') as w_f:
            fmt_dhs = "{chrom}\t{start}\t{end}\t{dhs_id}\t{cre}\t" \
                      "{dhs_score}\t{promoter_id}\t{score_h3k4me3}\t" \
                      "{p_h3k4me3}\t{score_h3k27ac}\t{p_h3k27ac}\n"
            for line in r_f:
                list_line = line.strip().split('\t')
                chrom = list_line[0]
                start = list_line[1]
                end = list_line[2]
                dhs_id = list_line[3]
                dhs_score = list_line[4]
                promoter_id = list_line[5]
                score_h3k4me3 = list_line[6]
                p_h3k4me3 = list_line[7]
                score_h3k27ac = list_line[8]
                p_h3k27ac = list_line[9]
                if (promoter_id != '.') & (p_h3k4me3 != '0') & \
                        (p_h3k27ac != '0'):
                    cre = 'Protein-Promoter(Enhancer)'
                elif (promoter_id == '.') & (p_h3k4me3 != '0') & \
                        (p_h3k27ac != '0'):
                    cre = 'Other-Promoter(Enhancer)'
                elif (promoter_id != '.') & (p_h3k4me3 != '0') & \
                        (p_h3k27ac == '0'):
                    cre = 'Protein-Promoter'
                elif (p_h3k4me3 == '0') & (p_h3k27ac != '0'):
                    cre = 'Enhancer'
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


def sub_annotate_cre_ctcf(dict_in):
    path = dict_in['path_out']
    file_in = os.path.join(path, 'DHS_promoter_H3K4me3_H3K27ac_CTCF.txt')
    if not os.path.exists(file_in):
        return
    file_out = os.path.join(path, 'cRE.CTCF.txt')
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
                if (promoter_id != '.') & (p_h3k4me3 != '0') & \
                        (p_h3k27ac != '0'):
                    cre = 'Protein-Promoter(Enhancer)'
                elif (promoter_id == '.') & (p_h3k4me3 != '0') & \
                        (p_h3k27ac != '0'):
                    cre = 'Other-Promoter(Enhancer)'
                elif (promoter_id != '.') & (p_h3k4me3 != '0') & \
                        (p_h3k27ac == '0'):
                    cre = 'Protein-Promoter'
                elif (p_h3k4me3 == '0') & (p_h3k27ac != '0'):
                    cre = 'Enhancer'
                elif (p_h3k27ac == '0') & (score_ctcf != '-10000'):
                    cre = 'Insulator'
                else:
                    cre = '.'
                w_f.write(fmt_dhs.format(**locals()))

    return


def annotate_cre(path_ref, path_h3k27ac, path_cre, num_process):
    # if os.path.exists(path_out_h3k27ac):
    #     os.system(f"rm -rf {path_out_h3k27ac}")
    # os.mkdir(path_out_h3k27ac)
    # os.system(f"cp {os.path.join(path_h3k27ac, 'metadata.simple.tsv')} "
    #           f"{os.path.join(path_out_h3k27ac, 'metadata.simple.tsv')}")
    # os.system(f"cp {os.path.join(path_h3k27ac, 'meta.reference.tsv')} "
    #           f"{os.path.join(path_out_h3k27ac, 'meta.reference.tsv')}")
    # os.system(
    #     f"cp {os.path.join(path_ref, 'meta.reference.histone.tsv')} "
    #     f"{os.path.join(path_out_h3k27ac, 'meta.reference.histone.tsv')}")

    df_ref_histone = pd.read_csv(
        os.path.join(path_ref, 'meta.reference.histone.tsv'), sep='\t')
    df_meta_h3k27ac = pd.read_csv(
        os.path.join(path_h3k27ac, 'metadata.simple.tsv'), sep='\t'
    )
    df_merge = pd.merge(
        df_ref_histone, df_meta_h3k27ac,
        on=['Biosample life stage', 'Biosample term name'])
    df_merge = df_merge.drop_duplicates('File accession')

    # map H3K27ac to sample
    # pool = Pool(processes=num_process)
    # func_map = partial(map_h3k27ac, path_ref, path_h3k27ac, path_out_h3k27ac)
    # pool.map(func_map, df_merge.to_dict('records'))
    # pool.close()

    # integrate H3K27ac by term and suborgan
    if os.path.exists(path_cre):
        os.system(f"rm -rf {path_cre}")
    os.mkdir(path_cre)
    os.system(f"cp {os.path.join(path_h3k27ac, 'metadata.simple.tsv')} "
              f"{os.path.join(path_cre, 'metadata.simple.tsv')}")
    os.system(f"cp {os.path.join(path_h3k27ac, 'meta.reference.tsv')} "
              f"{os.path.join(path_cre, 'meta.reference.tsv')}")
    os.system(f"cp {os.path.join(path_ref, 'meta.reference.histone.tsv')} "
              f"{os.path.join(path_cre, 'meta.reference.histone.tsv')}")
    life_organs = list(set(df_merge['Biosample life_organ'].tolist()))

    list_input = []
    list_ref = []
    for life_organ in life_organs:
        str_life_organ = life_organ.replace(' ', '_')
        sub_ref_histone = \
            df_ref_histone.loc[
             df_ref_histone['Biosample life_organ'] == life_organ, :
            ]
        path_life_organ = os.path.join(path_cre, str_life_organ)
        if not os.path.exists(path_life_organ):
            os.mkdir(path_life_organ)
        suborgans = list(set(sub_ref_histone['Biosample suborgan'].tolist()))
        for suborgan in suborgans:
            str_suborgan = suborgan.replace(' ', '_')
            suborgan_histone = sub_ref_histone.loc[
                    sub_ref_histone['Biosample suborgan'] == suborgan, :]
            terms = list(set(suborgan_histone['Biosample term name'].tolist()))
            if suborgan != 'single':
                path_suborgan = os.path.join(path_life_organ, str_suborgan)
                suborgan_h3k27ac = pd.merge(
                    suborgan_histone, df_meta_h3k27ac,
                    on=['Biosample life stage', 'Biosample term name']
                )
                suborgan_h3k27ac = \
                    suborgan_h3k27ac.drop_duplicates('File accession')
                if suborgan_h3k27ac.shape[0] == 0:
                    continue
                file_ref_suborgan = os.path.join(
                    path_ref,
                    f"{str_life_organ}/{str_suborgan}/DHS_promoter_H3K4me3.txt"
                )
                if os.path.exists(file_ref_suborgan):
                    list_input.append(dict(
                        file_ref=file_ref_suborgan,
                        path_out=path_suborgan,
                        sub_h3k27ac=suborgan_h3k27ac)
                    )
                    list_ref.append(
                        {'Biosample life_organ': life_organ,
                         'Biosample suborgan': suborgan,
                         'Biosample term name': 'empty',
                         'file_ref_h3k4me3': file_ref_suborgan,
                         'Level': 'suborgan'})
                else:
                    print(life_organ, suborgan)
            else:
                path_suborgan = path_life_organ
            if not os.path.exists(path_suborgan):
                os.mkdir(path_suborgan)
            for term in terms:
                term_histone = suborgan_histone.loc[
                    suborgan_histone['Biosample term name'] == term, :]
                str_term = term.replace(
                    ' ', '_').replace('/', '+').replace("'", '--')
                if suborgan != 'single':
                    file_ref = os.path.join(
                        path_ref,
                        f"{str_life_organ}/{str_suborgan}/{str_term}/"
                        f"DHS_promoter_H3K4me3.txt"
                    )
                else:
                    file_ref = os.path.join(
                        path_ref,
                        f"{str_life_organ}/{str_term}/DHS_promoter_H3K4me3.txt"
                    )
                if not os.path.exists(file_ref):
                    print(life_organ, term)
                    continue
                path_term = os.path.join(path_suborgan, str_term)
                if not os.path.exists(path_term):
                    os.mkdir(path_term)
                term_h3k27ac = pd.merge(
                    term_histone, df_meta_h3k27ac,
                    on=['Biosample life stage', 'Biosample term name']
                )
                term_h3k27ac = term_h3k27ac.drop_duplicates('File accession')
                if term_h3k27ac.shape[0] == 0:
                    continue
                list_input.append(dict(
                    file_ref=file_ref, path_out=path_term,
                    sub_h3k27ac=term_h3k27ac)
                )
                list_ref.append(
                    {'Biosample life_organ': life_organ,
                     'Biosample suborgan': suborgan,
                     'Biosample term name': term,
                     'file_ref_h3k4me3': file_ref, 'Level': 'term'})

    df_ref_h3k27ac = pd.DataFrame(list_ref)
    df_ref_h3k27ac.to_csv(
        os.path.join(path_cre, 'meta.reference.tsv'), sep='\t', index=None)

    pool = Pool(processes=num_process)
    func_integrate = partial(integrate_h3k27ac, path_h3k27ac)
    pool.map(func_integrate, list_input)
    pool.close()
    print('Annotation of H3K27ac is completed!')

    # integrate CTCF
    os.system(f"cp {os.path.join(path_ctcf_stan, 'metadata.simple.tsv')} "
              f"{os.path.join(path_cre, 'metadata.simple.CTCF.tsv')}")
    os.system(f"cp {os.path.join(path_ctcf_stan, 'meta.reference.tsv')} "
              f"{os.path.join(path_cre, 'meta.reference.CTCF.tsv')}")
    df_meta_ctcf = pd.read_csv(
        os.path.join(path_ctcf_stan, 'metadata.simple.tsv'), sep='\t')
    df_merge_copy = df_merge.loc[:,
                    ['Biosample life stage', 'Biosample term name',
                     'Biosample life_organ']].copy()
    df_merge_ctcf = pd.merge(
        df_merge_copy, df_meta_ctcf,
        on=['Biosample life stage', 'Biosample term name'])
    df_merge_ctcf = df_merge_ctcf.drop_duplicates('File accession')
    life_organs = list(set(df_merge_ctcf['Biosample life_organ'].tolist()))

    list_input_ctcf = []
    for life_organ in life_organs:
        str_life_organ = life_organ.replace(' ', '_')
        sub_ref_histone = \
            df_ref_histone.loc[
             df_ref_histone['Biosample life_organ'] == life_organ, :
            ]
        path_life_organ = os.path.join(path_cre, str_life_organ)
        if not os.path.exists(path_life_organ):
            os.mkdir(path_life_organ)
        suborgans = list(set(sub_ref_histone['Biosample suborgan'].tolist()))
        for suborgan in suborgans:
            str_suborgan = suborgan.replace(' ', '_')
            suborgan_histone = \
                sub_ref_histone.loc[
                 sub_ref_histone['Biosample suborgan'] == suborgan, :]
            terms = list(set(suborgan_histone['Biosample term name'].tolist()))
            if suborgan != 'single':
                path_suborgan = os.path.join(path_life_organ, str_suborgan)
                suborgan_ctcf = pd.merge(
                    suborgan_histone, df_merge_ctcf,
                    on=['Biosample life stage', 'Biosample term name']
                )
                suborgan_ctcf = \
                    suborgan_ctcf.drop_duplicates('File accession')
                if suborgan_ctcf.shape[0] == 0:
                    continue
                file_ref_suborgan = os.path.join(
                    path_cre,
                    f"{str_life_organ}/{str_suborgan}/"
                    f"DHS_promoter_H3K4me3_H3K27ac.txt"
                )
                if os.path.exists(file_ref_suborgan):
                    list_input_ctcf.append(dict(
                        file_ref=file_ref_suborgan,
                        path_out=path_suborgan,
                        sub_ctcf=suborgan_ctcf)
                    )
                else:
                    print(life_organ, suborgan)
            else:
                path_suborgan = path_life_organ
            if not os.path.exists(path_suborgan):
                os.mkdir(path_suborgan)
            for term in terms:
                term_histone = \
                    suborgan_histone.loc[
                        suborgan_histone['Biosample term name'] == term, :]
                str_term = term.replace(
                    ' ', '_').replace('/', '+').replace("'", '--')
                if suborgan != 'single':
                    file_ref = os.path.join(
                        path_cre,
                        f"{str_life_organ}/{str_suborgan}/{str_term}/"
                        f"DHS_promoter_H3K4me3_H3K27ac.txt"
                    )
                else:
                    file_ref = os.path.join(
                        path_cre,
                        f"{str_life_organ}/{str_term}/"
                        f"DHS_promoter_H3K4me3_H3K27ac.txt"
                    )
                if not os.path.exists(file_ref):
                    print(life_organ, term)
                    continue
                path_term = os.path.join(path_suborgan, str_term)
                if not os.path.exists(path_term):
                    os.mkdir(path_term)
                term_ctcf = pd.merge(
                    term_histone, df_meta_ctcf,
                    on=['Biosample life stage', 'Biosample term name']
                )
                term_ctcf = term_ctcf.drop_duplicates('File accession')
                if term_ctcf.shape[0] == 0:
                    continue
                list_input_ctcf.append(dict(
                    file_ref=file_ref, path_out=path_term,
                    sub_ctcf=term_ctcf)
                )

    pool = Pool(processes=num_process)
    func_integrate = partial(integrate_ctcf, path_ctcf_stan)
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

    # annotate DHS by term
    path_dhs_stan = '/local/zy/PEI/mid_data/tissue/DHS/GRCh38tohg19_standard'
    path_dhs_cluster = '/local/zy/PEI/mid_data/tissue/DHS/GRCh38tohg19_cluster'
    path_h3k4me3_stan = \
        '/local/zy/PEI/mid_data/tissue/ENCODE/histone_ChIP-seq/' \
        'H3K4me3_standard'
    path_h3k27ac_stan = \
        '/local/zy/PEI/mid_data/tissue/ENCODE/histone_ChIP-seq/' \
        'H3K27ac_standard'
    path_ctcf_stan = \
        '/local/zy/PEI/mid_data/tissue/ENCODE/TF_ChIP-seq/CTCF_standard'

    # promoter reference
    promoter_file_hg19 = \
        '/local/zy/PEI/origin_data/gene/promoters.up2k.protein.gencode.v19.bed'
    meta_suborgan_dhs = \
        '/local/zy/PEI/origin_data/meta_file/meta.reference.tsv'
    meta_suborgan_histone = \
        '/local/zy/PEI/origin_data/meta_file/meta.reference.histone.tsv'
    path_ref_promoter = '/local/zy/PEI/mid_data/tissue/DHS/reference_map'
    annotate_promoter_to_dhs(
        path_dhs_cluster, path_dhs_stan, path_h3k4me3_stan,
        promoter_file_hg19, meta_suborgan_histone, meta_suborgan_dhs,
        path_ref_promoter, num_cpu
    )
    print('Annotation of promoters and H3K4me3 is completed!')

    # map H3K27ac to reference
    protein_exon = \
        '/local/zy/PEI/origin_data/gene/exon.protein.gencode.v19.bed'
    # path_map_h3k27ac = '/local/zy/PEI/mid_data/tissue/DHS/map_H3K27ac'
    path_combine_h3k27ac = '/local/zy/PEI/mid_data/tissue/DHS/cRE_annotation'
    annotate_cre(path_ref_promoter, path_h3k27ac_stan, path_combine_h3k27ac,
                 num_cpu)
    print('Annotation of genome is completed!')

    time_end = time()
    print(time_end - time_start)
