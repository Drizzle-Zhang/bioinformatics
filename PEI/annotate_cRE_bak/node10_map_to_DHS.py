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


def annotate_promoter_to_dhs(path_cluster, path_h3k4me3,
                             loc_promoter, ref_histone, path_out, num_process):
    if os.path.exists(path_out):
        os.system(f"rm -rf {path_out}")
    os.mkdir(path_out)

    df_ref_histone = pd.read_csv(ref_histone, sep='\t')
    df_ref_histone = df_ref_histone.dropna()
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
            file_dhs = os.path.join(
                path_cluster,
                f"{str_organ}/{str_suborgan}/{str_suborgan}.bed"
            )
            suborgan_h3k4me3 = pd.merge(
                suborgan_histone, df_meta_h3k4me3,
                on=['Biosample life stage', 'Biosample term name']
            )
            list_input.append(dict(
                file_dhs=file_dhs, path_out=path_suborgan,
                path_h3k4me3=path_h3k4me3, sub_h3k4me3=suborgan_h3k4me3,
                loc_promoter=loc_promoter)
            )
            for sub_dict in suborgan_histone.to_dict("records"):
                term = sub_dict['Biosample term name']
                life = sub_dict['Biosample life stage']
                str_term = sub_dict['Biosample term name'].replace(
                    ' ', '_').replace('/', '+').replace("'", '--')
                path_term = os.path.join(path_suborgan, f"{life}_{str_term}")
                if not os.path.exists(path_term):
                    os.mkdir(path_term)
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


def map_h3k27ac(path_ref, path_h3k27ac, path_out, dict_in):
    accession_ids = dict_in['File accession']
    file_in = os.path.join(path_h3k27ac, accession_ids + '.bed')
    str_organ = dict_in['Biosample organ'].replace(' ', '_')
    str_suborgan = dict_in['Biosample suborgan'].replace(' ', '_')
    str_term = dict_in['Biosample term name'].replace(
                    ' ', '_').replace('/', '+').replace("'", '--')
    file_ref = os.path.join(
        path_ref,
        f"{str_organ}/{str_suborgan}/{dict_in['Biosample life stage']}_"
        f"{str_term}/DHS_promoter_H3K4me3.txt"
    )

    # map H3K4me3 to DHS
    file_plus = os.path.join(path_out, accession_ids + '.plus')
    file_uniq = os.path.join(path_out, accession_ids + '.uniq')
    file_sort = os.path.join(path_out, accession_ids + '.sort')

    os.system(
        f"bedtools intersect -a {file_ref} -b {file_in} -wao "
        f"| cut -f 1,2,3,4,5,6,10,11,14 > {file_plus}")
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
            file_accession, sep='\t', header=None, usecols=[0, 1, 2, 3, 6, 7])
        df_plus = pd.merge(df_ref, df_accession, on=[0, 1, 2, 3])
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
    os.system(f"Rscript adjust_p_value_H3K27ac.R "
              f"{file_origin} {file_out} {infer_num} {file_num}")

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
                      "{promoter_id}\t{score_h3k4me3}\t{score_h3k27ac}\n"
            for line in r_f:
                list_line = line.strip().split('\t')
                chrom = list_line[0]
                start = list_line[1]
                end = list_line[2]
                dhs_id = list_line[3]
                promoter_id = list_line[4]
                score_h3k4me3 = list_line[5]
                score_h3k27ac = list_line[6]
                if (promoter_id != '.') & (score_h3k4me3 != '0') & \
                        (score_h3k27ac != '0'):
                    cre = 'Promoter'
                elif score_h3k27ac != '0':
                    cre = 'Enhancer'
                else:
                    cre = '.'
                w_f.write(fmt_dhs.format(**locals()))

    return


def annotate_cre(path_ref, ref_histone, path_h3k27ac,
                 path_out_h3k27ac, path_cre, num_process):
    if os.path.exists(path_out_h3k27ac):
        os.system(f"rm -rf {path_out_h3k27ac}")
    os.mkdir(path_out_h3k27ac)
    os.system(f"cp {os.path.join(path_h3k27ac, 'metadata.simple.tsv')} "
              f"{os.path.join(path_out_h3k27ac, 'metadata.simple.tsv')}")
    os.system(f"cp {os.path.join(path_h3k27ac, 'meta.reference.tsv')} "
              f"{os.path.join(path_out_h3k27ac, 'meta.reference.tsv')}")

    df_ref_histone = pd.read_csv(ref_histone, sep='\t')
    df_ref_histone = df_ref_histone.dropna()
    df_meta_h3k27ac = pd.read_csv(
        os.path.join(path_h3k27ac, 'metadata.simple.tsv'), sep='\t'
    )
    df_merge = pd.merge(
        df_ref_histone, df_meta_h3k27ac,
        on=['Biosample organ', 'Biosample life stage', 'Biosample term name'])

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
    organs = list(
        set([organ for organ in df_ref_histone['Biosample organ'].tolist()])
    )

    list_input = []
    for organ in organs:
        str_organ = organ.replace(' ', '_')
        sub_ref_histone = df_ref_histone.loc[
            df_ref_histone['Biosample organ'] == organ, :
        ]
        path_organ = os.path.join(path_cre, organ.replace(' ', '_'))
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
            file_ref_suborgan = os.path.join(
                path_ref,
                f"{str_organ}/{str_suborgan}/DHS_promoter_H3K4me3.txt"
            )
            suborgan_h3k27ac = pd.merge(
                suborgan_histone, df_meta_h3k27ac,
                on=['Biosample life stage', 'Biosample term name']
            )
            list_input.append(dict(
                file_ref=file_ref_suborgan, path_out=path_suborgan,
                sub_h3k27ac=suborgan_h3k27ac)
            )
            for sub_dict in suborgan_histone.to_dict("records"):
                term = sub_dict['Biosample term name']
                life = sub_dict['Biosample life stage']
                str_term = sub_dict['Biosample term name'].replace(
                    ' ', '_').replace('/', '+').replace("'", '--')
                path_term = os.path.join(path_suborgan, f"{life}_{str_term}")
                file_ref = os.path.join(
                    path_ref,
                    f"{str_organ}/{str_suborgan}/{life}_{str_term}/"
                    f"DHS_promoter_H3K4me3.txt"
                )
                if not os.path.exists(path_term):
                    os.mkdir(path_term)
                sub_h3k27ac = df_meta_h3k27ac.loc[
                    (df_meta_h3k27ac['Biosample life stage'] == life) &
                    (df_meta_h3k27ac['Biosample term name'] == term), :]
                list_input.append(dict(
                    file_ref=file_ref, path_out=path_term,
                    sub_h3k27ac=sub_h3k27ac)
                )

    pool = Pool(processes=num_process)
    func_integrate = partial(integrate_h3k27ac, path_out_h3k27ac)
    pool.map(func_integrate, list_input)
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
    set_h3k27ac = set([
        f"{sub_dict['Biosample life stage']}_{sub_dict['Biosample term name']}"
        for sub_dict in df_h3k27ac.to_dict('records')
    ])
    set_intersect = set([
        f"{sub_dict['Biosample life stage']}_{sub_dict['Biosample term name']}"
        for sub_dict in df_intersect.to_dict('records')
    ])
    set_diff = set_h3k27ac.difference(set_intersect)
    # df_intersect.to_csv(file_meta, sep='\t', index=None)

    # promoter reference
    promoter_file_hg19 = \
        '/local/zy/PEI/data/gene/' \
        'promoters.up2k.protein.gencode.v19.merge.bed'
    meta_suborgan_dhs = '/local/zy/PEI/data/DHS/meta.reference.tsv'
    path_ref_promoter = '/local/zy/PEI/data/DHS/reference_map'
    annotate_promoter_to_dhs(
        path_dhs_cluster, path_h3k4me3_stan,
        promoter_file_hg19, file_meta, path_ref_promoter, num_cpu
    )

    # map H3K27ac to reference
    path_map_h3k27ac = '/local/zy/PEI/data/DHS/map_H3K27ac'
    path_combine_h3k27ac = '/local/zy/PEI/data/DHS/cRE_annotation'
    annotate_cre(path_ref_promoter, file_meta, path_h3k27ac_stan,
                 path_map_h3k27ac, path_combine_h3k27ac, num_cpu)

    time_end = time()
    print(time_end - time_start)
