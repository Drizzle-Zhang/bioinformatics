#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: preparation.py
# @time: 12/3/19 10:52 AM

from time import time
import pandas as pd
import numpy as np
import os
from multiprocessing import Pool
from functools import partial
import sys
from statsmodels.api import distributions
# sys.path.append('/local/zy/my_git/bioinformatics/PEI/annotate_cRE')
sys.path.append(
    '/lustre/tianlab/zhangyu/my_git/bioinformatics/PEI/annotate_cRE')
from preparation import \
    filter_meta, build_dict_attr, add_attr, hg38tohg19, \
    merge_peak_bed, overlap_matrix, merge_experiment, merge_standard_bed
root_path = '/lustre/tianlab/zhangyu/my_git/bioinformatics/PEI/annotate_cRE'


def modify_meta(df_meta):
    # only select cell line data
    df_meta = df_meta.loc[df_meta['Biosample type'] == 'cell line', :]

    return df_meta


def unique_bed_files(path_in, path_out, flank_percent, num_process):
    df_meta = pd.read_csv(
        os.path.join(path_in, 'metadata.simple.tsv'), sep='\t')

    if os.path.exists(path_out):
        os.system(f"rm -rf {path_out}")
    os.mkdir(path_out)

    os.system(f"cp {os.path.join(path_in, 'metadata.simple.tsv')} "
              f"{os.path.join(path_out, 'metadata.simple.tsv')}")

    meta_out = df_meta.loc[
               :, ['Biosample term name', 'Biosample life stage',
                   'Biosample organ', 'Biosample cell']]
    meta_out = meta_out.drop_duplicates()
    meta_out.to_csv(os.path.join(path_out, 'meta.reference.tsv'),
                    sep='\t', index=None)

    list_input = []
    terms = set(df_meta['Biosample term name'].tolist())
    for term in terms:
        term_meta = \
            df_meta.loc[df_meta['Biosample term name'] == term, :]
        life_stage = pd.unique(
            meta_out.loc[
                meta_out['Biosample term name'] == term,
                'Biosample life stage'])[0]
        accession_ids = \
            list(set(term_meta['Experiment accession'].tolist()))
        path_term = \
            os.path.join(path_out, term.replace(' ', '_').replace('/', '+'))
        if not os.path.exists(path_term):
            os.makedirs(path_term)
        list_input.append(
            dict(path=path_term,
                 term_name=term.replace(' ', '_').replace(
                     '/', '+').replace("'", '--'),
                 life_stage=life_stage,
                 accession_ids=accession_ids,
                 flank_percent=flank_percent))
        # if term == "HTR-8/SVneo":
        #     break

    merge_peak_bed(path_in, list_input, num_process)
    
    pool = Pool(processes=num_process)
    func_overlap = partial(overlap_matrix, path_in)
    list_df = pool.map(func_overlap, list_input)
    pool.close()

    df_overlap = pd.concat(list_df, sort=False)
    df_overlap.to_csv(
        os.path.join(path_out, 'overlap.txt'), sep='\t', index=None
    )

    return


def sub_stan(type_bed, path_in, path_out, dict_in):
    if type_bed == 'DHS':
        term = dict_in['Biosample term name']
        term_name = term.replace(' ', '_').replace('/', '+').replace("'", "--")
        file_label = term
        file = f"{term_name}/{term_name}.bed"
        folder3 = os.path.join(path_out, term_name)
        file_in = os.path.join(path_in, file)
        file_out_unsort = os.path.join(path_out, file + '.unsort')
        file_out = os.path.join(path_out, file)
        # make folder
        if not os.path.exists(folder3):
            os.mkdir(folder3)
        file_score = os.path.join(path_in, f"{term_name}/score.txt")
        file_quantile = os.path.join(folder3, 'quantile.txt')
        os.system(
            f"Rscript {os.path.join(root_path, 'correct_DHS_score.R')} "
            f"{file_score} {file_quantile}")
        df_quantile = pd.read_csv(
            file_quantile, sep='\t', header=None, index_col=0,
            names=['dhs_id', 'quantile'])
        set_dhs_id = set(df_quantile.index)
    else:
        file_in = os.path.join(path_in, dict_in['File accession'] + '.bed')
        # normalization
        df_bed = pd.read_csv(file_in, sep='\t', header=None, usecols=[6])
        ecdf = distributions.ECDF(df_bed.iloc[:, 0])
        df_quantile = pd.Series(ecdf(df_bed.iloc[:, 0]), index=df_bed.index)
        file_out_unsort = \
            os.path.join(path_out, dict_in['File accession'] + '.unsort.bed')
        file_out = os.path.join(path_out, dict_in['File accession'] + '.bed')

    with open(file_in, 'r') as r_f:
        with open(file_out_unsort, 'w') as w_f:
            fmt_dhs = "{chrom}\t{start}\t{end}\t{label}\t{score}\t.\t" \
                      "{file_label}\t{accessions}\n"
            fmt_chip = "{chrom}\t{start}\t{end}\t{label}\t" \
                       "{score}\t{pvalue}\t{accessions}\n"
            for idx, line in enumerate(r_f):
                list_line = line.strip().split('\t')
                chrom = list_line[0]
                start = list_line[1]
                end = list_line[2]
                accessions = list_line[3]
                label = f"{type_bed}<-{chrom}:{start}-{end}"
                if type_bed == 'DHS':
                    dhs_id = f"{chrom}:{start}-{end}"
                    if dhs_id in set_dhs_id:
                        score = df_quantile.loc[dhs_id, 'quantile']
                    else:
                        continue
                    w_f.write(fmt_dhs.format(**locals()))
                elif type_bed in {'H3K4me3', 'H3K27ac'}:
                    score = df_quantile.loc[idx]
                    pvalue = list_line[7]
                    w_f.write(fmt_chip.format(**locals()))
                elif type_bed == 'CTCF':
                    score = df_quantile.loc[idx]
                    pvalue = list_line[8]
                    w_f.write(fmt_chip.format(**locals()))

    os.system(f"bedtools sort -i {file_out_unsort} > {file_out}")
    os.remove(file_out_unsort)

    return


def standardize_bed(path_in, path_out, type_bed, num_process):
    if os.path.exists(path_out):
        os.system(f"rm -rf {path_out}")
    os.mkdir(path_out)
    os.system(f"cp {os.path.join(path_in, 'metadata.simple.tsv')} "
              f"{os.path.join(path_out, 'metadata.simple.tsv')}")
    if os.path.exists(os.path.join(path_in, 'meta.reference.tsv')):
        os.system(f"cp {os.path.join(path_in, 'meta.reference.tsv')} "
                  f"{os.path.join(path_out, 'meta.reference.tsv')}")
    else:
        df_meta = pd.read_csv(
            os.path.join(path_in, 'metadata.simple.tsv'), sep='\t')
        meta_out = df_meta.loc[
                   :, ['Biosample term name', 'Biosample life stage',
                       'Biosample organ', 'Biosample cell']]
        meta_out = meta_out.drop_duplicates()
        meta_out.to_csv(os.path.join(path_out, 'meta.reference.tsv'),
                        sep='\t', index=None)

    if type_bed == 'DHS':
        df_meta = pd.read_csv(
            os.path.join(path_in, 'meta.reference.tsv'), sep='\t'
        )
        df_meta = df_meta['Biosample term name'].drop_duplicates()
        list_input = []
        list_term = df_meta.tolist()
        for term in list_term:
            list_input.append({'Biosample term name': term})
            folder1 = os.path.join(
                path_out, f"{term.replace(' ', '_').replace('/', '+')}")
            if not os.path.exists(folder1):
                os.mkdir(folder1)

    else:
        df_meta = pd.read_csv(
            os.path.join(path_in, 'metadata.simple.tsv'), sep='\t'
        )
        list_input = df_meta.to_dict('records')

    pool = Pool(processes=num_process)
    func_stan = partial(sub_stan, type_bed, path_in, path_out)
    pool.map(func_stan, list_input)
    pool.close()

    return


def merge_all_cells(path_stan, num_process):
    df_meta = pd.read_csv(
        os.path.join(path_stan, 'meta.reference.tsv'), sep='\t', usecols=[0])
    accession_ids = \
        [f"{term.replace(' ', '_').replace('/', '+')}/"
         f"{term.replace(' ', '_').replace('/', '+')}"
         for term in (df_meta['Biosample term name'].unique()).tolist()]

    dict_merge = [dict(
        path=path_stan,
        term_name='all_celllines',
        accession_ids=accession_ids,
        flank_percent=0.5)]
    merge_standard_bed(path_stan, dict_merge, num_process)

    return


if __name__ == '__main__':
    time_start = time()
    # parameters
    num_cpu = 40
    path_root = '/lustre/tianlab/zhangyu/PEI'
    path_origin = path_root + '/origin_data'
    path_mid = path_root + '/mid_data_correct'

    # build life stage dictionary
    path_lifestage = path_origin + '/ENCODE/metadata/life_stage'
    dict_lifestage = build_dict_attr(path_lifestage)

    # build organ dictionary
    path_meta_organ = path_origin + '/ENCODE/metadata/organ'
    dict_organ = build_dict_attr(path_meta_organ)

    # build organ dictionary
    path_meta_cell = path_origin + '/ENCODE/metadata/cell'
    dict_cell = build_dict_attr(path_meta_cell)

    # build organ dictionary
    path_meta_lab = path_origin + '/ENCODE/metadata/lab'
    dict_lab = build_dict_attr(path_meta_lab)
    print("Preparation of dictionary files and reference files is completed")

    # DHS reference
    # metafile
    path_dhs = path_origin + '/ENCODE/DNase-seq/all'
    ori_meta_dhs = os.path.join(path_dhs, 'metadata.tsv')
    df_meta_dhs = filter_meta(ori_meta_dhs)
    df_meta_dhs = add_attr(df_meta_dhs, dict_lifestage, 'Biosample life stage')
    df_meta_dhs = add_attr(df_meta_dhs, dict_organ, 'Biosample organ')
    df_meta_dhs = add_attr(df_meta_dhs, dict_cell, 'Biosample cell')
    df_meta_dhs = add_attr(df_meta_dhs, dict_lab, 'Lab')
    df_meta_dhs = modify_meta(df_meta_dhs)
    meta_dhs = path_mid + '/cell_line/metadata.simple.tsv'
    df_meta_dhs.to_csv(meta_dhs, sep='\t', index=None)
    print("DHS metadata ---- completed")

    # hg38 to hg19
    path_hg38tohg19 = path_mid + '/cell_line/ENCODE/DNase-seq/GRCh38tohg19'
    hg38tohg19(path_dhs, path_hg38tohg19, meta_dhs, num_cpu)
    print("Format conversion: hg38 -> hg19 ---- completed")

    # integrate files from same experiment
    path_exp_dhs = \
        path_mid + '/cell_line/ENCODE/DNase-seq/GRCh38tohg19_experiment'
    merge_experiment(path_hg38tohg19, path_exp_dhs, 0.5, num_cpu)
    print("Integration of files from same experiment ---- completed")

    # build DHS reference
    path_dhs_hg38tohg19 = path_mid + '/cell_line/DHS/GRCh38tohg19/'
    unique_bed_files(path_exp_dhs, path_dhs_hg38tohg19, 0.5, num_cpu)
    print("Integration of files from same term ---- completed")

    # standardization
    path_dhs_stan = path_mid + '/cell_line/DHS/GRCh38tohg19_standard'
    standardize_bed(path_dhs_hg38tohg19, path_dhs_stan, 'DHS', num_cpu)
    print('Standardization of DHS completed!')

    # merge dhs files from all cell lines
    merge_all_cells(path_dhs_stan, num_cpu)
    print('Merge of DHS completed!')

    # preparation of bed files of histone and TF
    # H3K4me3
    path_h3k4me3 = path_origin + '/ENCODE/histone_ChIP-seq/H3K4me3'
    ori_meta_h3k4me3 = os.path.join(path_h3k4me3, 'metadata.tsv')
    df_meta_h3k4me3 = pd.read_csv(ori_meta_h3k4me3, sep='\t')
    df_meta_h3k4me3 = df_meta_h3k4me3.loc[
        (df_meta_h3k4me3['File Status'] == 'released') &
        ((df_meta_h3k4me3['Output type'] == 'replicated peaks') |
         (df_meta_h3k4me3['Output type'] == 'stable peaks')) &
        (df_meta_h3k4me3['Biosample treatments'].apply(
            lambda x: np.isnan(x) if isinstance(x, float) else False)) &
        (df_meta_h3k4me3['Biosample genetic modifications methods'].apply(
            lambda x: np.isnan(x) if isinstance(x, float) else False)),
        ['File accession', 'Experiment accession', 'Biosample term id',
         'Biosample term name', 'Biosample type', 'Biosample treatments',
         'Biosample genetic modifications methods', 'Output type', 'Assembly',
         'Biological replicate(s)', 'File Status']]
    df_meta_h3k4me3.index = df_meta_h3k4me3['File accession']
    df_meta_h3k4me3 = \
        add_attr(df_meta_h3k4me3, dict_lifestage, 'Biosample life stage')
    df_meta_h3k4me3 = add_attr(df_meta_h3k4me3, dict_organ, 'Biosample organ')
    df_meta_h3k4me3 = add_attr(df_meta_h3k4me3, dict_cell, 'Biosample cell')
    df_meta_h3k4me3 = add_attr(df_meta_h3k4me3, dict_lab, 'Lab')
    df_meta_h3k4me3 = modify_meta(df_meta_h3k4me3)
    meta_h3k4me3 = path_mid + '/cell_line/ENCODE/' \
                              'histone_ChIP-seq/metadata.simple.H3K4me3.tsv'
    df_meta_h3k4me3.to_csv(meta_h3k4me3, sep='\t', index=None)
    print("H3K4me3 metadata ---- completed")

    # hg38 to hg19
    path_hg38tohg19 = path_mid + '/cell_line/ENCODE/histone_ChIP-seq/H3K4me3'
    hg38tohg19(path_h3k4me3, path_hg38tohg19, meta_h3k4me3, num_cpu)
    print("Format conversion: hg38 -> hg19 ---- completed")

    # standardization
    path_h3k4me3_stan = \
        path_mid + '/cell_line/ENCODE/histone_ChIP-seq/H3K4me3_standard'
    standardize_bed(path_hg38tohg19, path_h3k4me3_stan, 'H3K4me3', num_cpu)
    print('Standardization of H3K4me3 completed!')

    # H3K27ac
    path_h3k27ac = path_origin + '/ENCODE/histone_ChIP-seq/H3K27ac'
    ori_meta_h3k27ac = os.path.join(path_h3k27ac, 'metadata.tsv')
    df_meta_h3k27ac = pd.read_csv(ori_meta_h3k27ac, sep='\t')
    df_meta_h3k27ac = df_meta_h3k27ac.loc[
        (df_meta_h3k27ac['File Status'] == 'released') &
        ((df_meta_h3k27ac['Output type'] == 'replicated peaks') |
         (df_meta_h3k27ac['Output type'] == 'stable peaks')) &
        (df_meta_h3k27ac['Biosample treatments'].apply(
            lambda x: np.isnan(x) if isinstance(x, float) else False)) &
        (df_meta_h3k27ac['Biosample genetic modifications methods'].apply(
            lambda x: np.isnan(x) if isinstance(x, float) else False)),
        ['File accession', 'Experiment accession', 'Biosample term id',
         'Biosample term name', 'Biosample type', 'Biosample treatments',
         'Biosample genetic modifications methods', 'Output type', 'Assembly',
         'Biological replicate(s)', 'File Status']]
    df_meta_h3k27ac.index = df_meta_h3k27ac['File accession']
    df_meta_h3k27ac = \
        add_attr(df_meta_h3k27ac, dict_lifestage, 'Biosample life stage')
    df_meta_h3k27ac = add_attr(df_meta_h3k27ac, dict_organ, 'Biosample organ')
    df_meta_h3k27ac = add_attr(df_meta_h3k27ac, dict_cell, 'Biosample cell')
    df_meta_h3k27ac = add_attr(df_meta_h3k27ac, dict_lab, 'Lab')
    df_meta_h3k27ac = modify_meta(df_meta_h3k27ac)
    meta_h3k27ac = path_mid + '/cell_line/ENCODE/' \
                              'histone_ChIP-seq/metadata.simple.H3K27ac.tsv'
    df_meta_h3k27ac.to_csv(meta_h3k27ac, sep='\t', index=None)
    print("H3K27ac metadata ---- completed")

    # hg38 to hg19
    path_hg38tohg19 = path_mid + '/cell_line/ENCODE/histone_ChIP-seq/H3K27ac'
    hg38tohg19(path_h3k27ac, path_hg38tohg19, meta_h3k27ac, num_cpu)
    print("Format conversion: hg38 -> hg19 ---- completed")

    # standardization
    path_h3k27ac_stan = \
        path_mid + '/cell_line/ENCODE/histone_ChIP-seq/H3K27ac_standard'
    standardize_bed(path_hg38tohg19, path_h3k27ac_stan, 'H3K27ac', num_cpu)
    print('Standardization of H3K27ac completed!')

    # CTCF
    path_ctcf = path_origin + '/ENCODE/TF_ChIP-seq/CTCF'
    ori_meta_ctcf = os.path.join(path_ctcf, 'metadata.tsv')
    df_meta_ctcf = pd.read_csv(ori_meta_ctcf, sep='\t')
    df_meta_ctcf = df_meta_ctcf.loc[
        (df_meta_ctcf['File Status'] == 'released') &
        ((df_meta_ctcf['Output type'] == 'optimal IDR thresholded peaks') |
         (df_meta_ctcf['Output type'] ==
          'pseudoreplicated IDR thresholded peaks')) &
        (df_meta_ctcf['Biosample treatments'].apply(
            lambda x: np.isnan(x) if isinstance(x, float) else False)) &
        (df_meta_ctcf['Biosample genetic modifications methods'].apply(
            lambda x: np.isnan(x) if isinstance(x, float) else False)),
        ['File accession', 'Experiment accession', 'Biosample term id',
         'Biosample term name', 'Biosample type', 'Biosample treatments',
         'Biosample genetic modifications methods', 'Output type', 'Assembly',
         'Biological replicate(s)', 'File Status']]
    df_meta_ctcf.index = df_meta_ctcf['File accession']
    df_meta_ctcf = \
        add_attr(df_meta_ctcf, dict_lifestage, 'Biosample life stage')
    df_meta_ctcf = add_attr(df_meta_ctcf, dict_organ, 'Biosample organ')
    df_meta_ctcf = add_attr(df_meta_ctcf, dict_cell, 'Biosample cell')
    df_meta_ctcf = add_attr(df_meta_ctcf, dict_lab, 'Lab')
    df_meta_ctcf = modify_meta(df_meta_ctcf)
    meta_ctcf = \
        path_mid + '/cell_line/ENCODE/TF_ChIP-seq/metadata.simple.CTCF.tsv'
    df_meta_ctcf.to_csv(meta_ctcf, sep='\t', index=None)
    print("CTCF metadata ---- completed")

    # hg38 to hg19
    path_hg38tohg19 = path_mid + '/cell_line/ENCODE/TF_ChIP-seq/CTCF'
    hg38tohg19(path_ctcf, path_hg38tohg19, meta_ctcf, num_cpu)
    print("Format conversion: hg38 -> hg19 ---- completed")

    # standardization
    path_ctcf_stan = path_mid + '/cell_line/ENCODE/TF_ChIP-seq/CTCF_standard'
    standardize_bed(path_hg38tohg19, path_ctcf_stan, 'CTCF', num_cpu)
    print('Standardization of CTCF completed!')

    time_end = time()
    print(time_end - time_start)
