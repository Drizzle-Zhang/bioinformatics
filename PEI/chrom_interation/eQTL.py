#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: eQTL.py
# @time: 2020/3/4 13:25

from time import time
import re
import os
import pandas as pd
from multiprocessing import Pool
import numpy as np
from subprocess import check_output
from functools import partial
import glob


def sub_uniform_eqtl(dict_in):
    file_in = dict_in['sub_file_in']
    file_out = dict_in['sub_file_out']
    tss_upstream = dict_in['tss_upstream']
    flank = dict_in['flank']
    pattern_start = re.compile(r'_[0-9]+_')
    pattern_chrom = re.compile(r'.+?_')
    with open(file_out, 'w') as w_out:
        fmt = "{chrom1}\t{start1}\t{end1}\t{variant_id}\t" \
              "{gene_symbol}\t{ensg}\t{tss}\t{tss_distance}\t{eqtl_params}\n"
        with open(file_in, 'r') as r_in:
            for line in r_in:
                list_line = line.strip().split('\t')
                variant_id = list_line[0]
                if variant_id == 'variant_id':
                    continue
                ensg = list_line[1]
                start = int(pattern_start.search(variant_id).group()[1:-1])
                df_ensg = df_protein_genes.loc[
                          df_protein_genes['gene_id'] == ensg, :]
                if df_ensg.shape[0] != 1:
                    continue
                chrom1 = f"chr{pattern_chrom.match(variant_id).group()[:-1]}"
                start1 = str(start - flank)
                end1 = str(start + flank)
                eqtl_params = ';'.join(list_line[3:])
                gene_symbol = df_ensg['gene_name'].iloc[0]
                if df_ensg['strand'].iloc[0] == '+':
                    tss = df_ensg['start'].iloc[0]
                else:
                    tss = df_ensg['end'].iloc[0]
                tss_distance = start - int(tss)
                if abs(tss_distance) < tss_upstream:
                    continue
                w_out.write(fmt.format(**locals()))

    return


def transform_eqtl(file_in, file_cre, file_out, file_pair_eqtl, file_pair_cre,
                   tss_upstream=5000, flank=1000):
    path_out, _ = os.path.split(file_out)
    path_subfile_in = os.path.join(path_out, 'subfiles_in')
    if not os.path.exists(path_subfile_in):
        os.mkdir(path_subfile_in)
    path_subfile_out = os.path.join(path_out, 'subfiles_out')
    if not os.path.exists(path_subfile_out):
        os.mkdir(path_subfile_out)
    os.system(f"split -a 4 -d -l 10000 {file_in} {path_subfile_in}/subfile")

    list_dict_in = []
    for i, subfile in enumerate(glob.glob(path_subfile_in + '/*')):
        list_dict_in.append(
            dict(sub_file_in=subfile, tss_upstream=tss_upstream, flank=flank,
                 sub_file_out=f"{path_subfile_out}/subfile{i}")
        )

    pool = Pool(40)
    pool.map(sub_uniform_eqtl, list_dict_in)
    pool.close()

    os.system(
        f"cat {' '.join(glob.glob(path_subfile_out + '/*'))} > {file_out}")
    os.system(f"rm -rf {path_subfile_in} {path_subfile_out}")

    os.system(f"bedtools intersect -a {file_out} -b {file_cre} -loj | "
              f"cut -f 4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20 | "
              f"grep -w 'Promoter\\|Enhancer' > {file_pair_eqtl}")
    df_eqtl = pd.read_csv(file_pair_eqtl, sep='\t', header=None)
    df_eqtl = df_eqtl.drop_duplicates()
    dhs_center = np.mean(df_eqtl[[6, 7]], axis=1)
    dhs_center = np.round(dhs_center)
    df_eqtl[16] = np.abs(df_eqtl[3] - dhs_center)
    df_eqtl.to_csv(file_pair_eqtl, sep='\t', header=None, index=None)
    df_eqtl = df_eqtl.iloc[:, [1] + list(range(8, 17))]
    df_eqtl = df_eqtl.drop_duplicates()
    df_eqtl.to_csv(file_pair_cre, sep='\t', header=None, index=None)

    return


def sub_hg18tohg19(path_hg18, path_hg19, file_name):
    file_hg18 = os.path.join(path_hg18, file_name)
    file_hg19 = os.path.join(path_hg19, file_name)
    file_hg18_1 = os.path.join(path_hg19, file_name + '.hg18.1')
    file_hg18_2 = os.path.join(path_hg19, file_name + '.hg18.2')
    file_hg19_1 = os.path.join(path_hg19, file_name + '.hg19.1')
    file_hg19_2 = os.path.join(path_hg19, file_name + '.hg19.2')
    pattern_chrom = re.compile(r'ld_.+?_')
    chrom = pattern_chrom.match(file_name).group()[3:-1]

    w_hg18_1 = open(file_hg18_1, 'w')
    w_hg18_2 = open(file_hg18_2, 'w')
    fmt_1 = "{chrom}\t{start1}\t{end1}\t{ld_id}\t{score}\t.\n"
    fmt_2 = "{chrom}\t{start2}\t{end2}\t{ld_id}\t{score}\t.\n"
    with open(file_hg18, 'r') as r_hg18:
        for line in r_hg18:
            list_line = line.strip().split(' ')
            start1 = int(list_line[0])
            end1 = start1 + 1
            start2 = int(list_line[1])
            end2 = start2 + 1
            ld_id = f"{list_line[3]}-{list_line[4]}"
            score = list_line[6]
            w_hg18_1.write(fmt_1.format(**locals()))
            w_hg18_2.write(fmt_2.format(**locals()))

    w_hg18_1.close()
    w_hg18_2.close()

    file_chain = \
        '/local/zy/tools/files_liftOver/hg18ToHg19.over.chain.gz'
    file_ummap_1 = os.path.join(path_hg19, file_name + '.unmap.1')
    file_ummap_2 = os.path.join(path_hg19, file_name + '.unmap.2')
    os.system(f"/local/zy/tools/liftOver {file_hg18_1} {file_chain} "
              f"{file_hg19_1} {file_ummap_1}")
    os.system(f"/local/zy/tools/liftOver {file_hg18_2} {file_chain} "
              f"{file_hg19_2} {file_ummap_2}")

    file_hg19_sort_1 = os.path.join(path_hg19, file_name + '.hg19.sort.1')
    file_hg19_sort_2 = os.path.join(path_hg19, file_name + '.hg19.sort.2')
    os.system(f"sort -k 4 -t '\t' {file_hg19_1} > {file_hg19_sort_1}")
    os.system(f"sort -k 4 -t '\t' {file_hg19_2} > {file_hg19_sort_2}")

    file_hg19_join = file_hg19 + '.hg19.join'
    df_hg19_1 = pd.read_csv(file_hg19_sort_1, sep='\t', header=None,
                            dtype={1: np.str})
    df_hg19_1['site1'] = df_hg19_1.apply(lambda x: x[0] + '_' + x[1], axis=1)
    df_hg19_2 = pd.read_csv(file_hg19_sort_2, sep='\t', header=None,
                            dtype={1: np.str})
    df_hg19_2['site2'] = df_hg19_2.apply(lambda x: x[0] + '_' + x[1], axis=1)
    df_merge = pd.merge(df_hg19_1, df_hg19_2, on=[3, 4])
    df_merge = df_merge.loc[:, ['site1', 'site2', 3, 4]]
    df_merge.to_csv(file_hg19, sep='\t', header=None, index=None)

    os.remove(file_hg18_1)
    os.remove(file_hg18_2)
    os.remove(file_hg19_1)
    os.remove(file_hg19_2)
    os.remove(file_ummap_1)
    os.remove(file_ummap_2)
    os.remove(file_hg19_sort_1)
    os.remove(file_hg19_sort_2)

    return


def ld_hg18tohg19(path_ld, race, path_out):
    files = os.listdir(path_ld)
    files_input = []
    files_output = []
    for file in files:
        if file[:2] != 'ld':
            continue
        file_race = file.split('_')[2][:-4]
        if file_race == race:
            files_input.append(file)
            files_output.append(os.path.join(path_out, file))

    pool = Pool(23)
    func_hg18tohg19 = partial(sub_hg18tohg19, path_ld, path_out)
    pool.map(func_hg18tohg19, files_input)
    pool.close()

    file_ld = os.path.join(path_out, f"ld_{race}.txt")
    os.system(f"cat {' '.join(files_output)} > {file_ld}")

    return


def negative_eqtl_sample(file_sig_eqtl, file_all_eqtl, file_sample_eqtl,
                         file_ld, cutoff_rs=0.2):
    # calculate nominal p-value cutoff
    df_sig_eqtl = pd.read_csv(file_sig_eqtl, sep='\t')
    len_sig_eqtl = df_sig_eqtl.shape[0]

    # ld data
    df_ld = pd.read_csv(file_ld, sep='\t', header=None)
    df_ld = df_ld.loc[df_ld[3] > cutoff_rs, :]

    # delete some snp (ld with eQTL)
    sig_snps = set((df_sig_eqtl['variant_id'].apply(
        lambda x: 'chr' + '_'.join(x.split('_')[:2]))).tolist())

    snps_1 = (df_ld.loc[df_ld[0].apply(lambda x: x in sig_snps), 1]).tolist()
    snps_2 = (df_ld.loc[df_ld[1].apply(lambda x: x in sig_snps), 0]).tolist()
    set_snps = set(snps_1 + snps_2)

    df_all_distance = pd.read_csv(file_all_eqtl, sep='\t')
    all_cols = ['variant_id', 'gene_id', 'tss_distance', 'ma_samples',
                'ma_count', 'maf', 'pval_nominal', 'slope', 'slope_se']
    df_all_distance = df_all_distance[all_cols]
    df_all_distance_filter = df_all_distance.loc[
        (df_all_distance['slope_se'] < 0.1) &
        (np.abs(df_all_distance['slope']) < 0.05), :]
    df_all_distance_filter = \
        df_all_distance_filter.loc[
         df_all_distance_filter['gene_id'].apply(
             lambda x: x in ensgs_positive), :]

    # sample
    # num_sample = 25*len_sig_eqtl
    # df_all_distance_sample = \
    #     df_all_distance_filter.sample(n=num_sample, random_state=1234)
    df_all_distance_sample = df_all_distance_filter.loc[
                             df_all_distance_filter['variant_id'].apply(
                                 lambda x: ('chr' + '_'.join(
                                     x.split('_')[:2])) not in set_snps
                             ), :]
    df_all_distance_sample = df_all_distance_sample.sort_index()
    df_all_distance_sample.to_csv(
        file_sample_eqtl, sep='\t', header=None, index=None)

    return


def negative_eqtl_select(file_pair_cre, file_pair_cre_neg, file_sample_neg,
                         times=1.0):
    df_positive = pd.read_csv(file_pair_cre, sep='\t', header=None,
                              usecols=[0, 1, 9])
    num_pos = df_positive.shape[0]
    num_select = times * num_pos
    df_negative = pd.read_csv(file_pair_cre_neg, sep='\t', header=None)
    # df_negative = \
    #     df_negative.loc[
    #      df_negative[0].apply(lambda x: x in genes_positive), :]
    df_negative = df_negative.sample(num_select, random_state=123)

    # list_df = []
    # for i in range(0, 100):
    #     down_limit = i * len_bin
    #     up_limit = (i + 1) * len_bin
    #     sub_pos = df_positive.loc[
    #               (df_positive[9] > down_limit) & (df_positive[9] <= up_limit),
    #               :]
    #     num_sample = int(num_select * (sub_pos.shape[0]/num_pos))
    #     sub_neg = df_negative.loc[
    #               (df_negative[9] > down_limit) & (df_negative[9] <= up_limit),
    #               :]
    #     sub_neg_sample = sub_neg.sample(
    #         min(num_sample, sub_neg.shape[0]), random_state=123)
    #     list_df.append(sub_neg_sample)
    #
    # df_negative = pd.concat(list_df, sort=False)
    df_negative.to_csv(file_sample_neg, sep='\t', header=None, index=None)

    return


if __name__ == '__main__':
    time_start = time()
    # hyper-params
    len_bin = 10000
    # protein coding genes
    protein_gene = \
        '/local/zy/PEI/origin_data/gene/genes.protein.gencode.v19.bed'
    df_protein_genes = pd.read_csv(protein_gene, sep='\t', header=None)
    df_protein_genes.columns = ['chrom', 'start', 'end',
                                'gene_name', 'gene_id', 'strand']

    # GM12878
    file_sig_eqtl_gm = \
        '/local/zy/PEI/origin_data/GTEx/GTEx_Analysis_v7_eQTL/' \
        'Cells_EBV-transformed_lymphocytes.v7.signif_variant_gene_pairs.txt'
    file_cre_gm = '/local/zy/PEI/mid_data/cell_line/DHS/cRE_annotation/' \
                  'GM12878/cRE.txt'
    file_uniform_gm_pos = \
        '/local/zy/PEI/mid_data/cell_line/eQTL/GM12878/positive/uniform.txt'
    file_eqtl_pair_gm_pos = \
        '/local/zy/PEI/mid_data/cell_line/eQTL/GM12878/positive/' \
        'eQTL_cRE_pairs.txt'
    file_cre_pair_gm_pos = \
        '/local/zy/PEI/mid_data/cell_line/eQTL/GM12878/positive/cRE_pairs.txt'
    transform_eqtl(file_sig_eqtl_gm, file_cre_gm, file_uniform_gm_pos,
                   file_eqtl_pair_gm_pos, file_cre_pair_gm_pos)

    # LD data
    path_ld_data = '/local/zy/PEI/origin_data/HapMap/LD_data'
    path_ld_treated = '/local/zy/PEI/mid_data/LD_data/hg19'
    # ld_hg18tohg19(path_ld_data, 'CEU', path_ld_treated)
    file_ld_ceu = '/local/zy/PEI/mid_data/LD_data/hg19/ld_CEU.txt'

    # positive genes
    df_pos_eqtl = pd.read_csv(file_eqtl_pair_gm_pos, sep='\t', header=None)
    ensgs_positive = set(df_pos_eqtl[2].tolist())
    genes_positive = set(df_pos_eqtl[1].tolist())

    # transform negative eQTL
    file_all_eqtl_gm = \
        '/local/zy/PEI/origin_data/GTEx/' \
        'GTEx_Analysis_v7_eQTL_all_associations/' \
        'Cells_EBV-transformed_lymphocytes.allpairs.txt'
    file_sample_eqtl_gm = \
        '/local/zy/PEI/mid_data/cell_line/eQTL/GM12878/negative/' \
        'allpairs.sample.txt'
    negative_eqtl_sample(
        file_sig_eqtl_gm, file_all_eqtl_gm, file_sample_eqtl_gm, file_ld_ceu)
    file_uniform_gm_neg = \
        '/local/zy/PEI/mid_data/cell_line/eQTL/GM12878/negative/uniform.txt'
    file_eqtl_pair_gm_neg = \
        '/local/zy/PEI/mid_data/cell_line/eQTL/GM12878/negative/' \
        'eQTL_cRE_pairs.txt'
    file_cre_pair_gm_neg = \
        '/local/zy/PEI/mid_data/cell_line/eQTL/GM12878/negative/cRE_pairs.txt'
    transform_eqtl(file_sample_eqtl_gm, file_cre_gm, file_uniform_gm_neg,
                   file_eqtl_pair_gm_neg, file_cre_pair_gm_neg)

    file_sample_gm_neg = \
        '/local/zy/PEI/mid_data/cell_line/eQTL/GM12878/negative/' \
        'cRE_pairs.sample.1.txt'
    negative_eqtl_select(file_cre_pair_gm_pos, file_cre_pair_gm_neg,
                         file_sample_gm_neg, times=1)
    file_sample_gm_neg_5 = \
        '/local/zy/PEI/mid_data/cell_line/eQTL/GM12878/negative/' \
        'cRE_pairs.sample.5.txt'
    negative_eqtl_select(file_cre_pair_gm_pos, file_cre_pair_gm_neg,
                         file_sample_gm_neg_5, times=5)

    time_end = time()
    print(time_end - time_start)

    # analysis
    # file_eqtl_pair_gm_neg = \
    #     '/local/zy/PEI/mid_data/cell_line/eQTL/GM12878/negative/' \
    #     'eQTL_cRE_pairs.txt'
    # file_eqtl_pair_gm_pos = \
    #     '/local/zy/PEI/mid_data/cell_line/eQTL/GM12878/positive/' \
    #     'eQTL_cRE_pairs.txt'
    # df_try = pd.read_csv(file_eqtl_pair_gm_pos, sep='\t', header=None)
    # file_pchic = \
    #     "/local/zy/PEI/mid_data/cell_line/gene_cre_pairs/GM12878/pcHi-C/" \
    #     "Jung_NG_2019/pairs.gene.cRE.txt"
    # df_pchic = pd.read_csv(file_pchic, sep='\t', header=None, usecols=[0, 1])
    # df_pchic[3] = np.ones(df_pchic.shape[0], dtype='int')
    # df_merge = pd.merge(df_try, df_pchic, left_on=[1, 8], right_on=[0, 1],
    #                     how='left')
    # df_merge = df_merge.fillna(0)
    # df_merge = df_merge.drop_duplicates()
    # df_n = df_merge.loc[df_merge[5] > 0.9, :]
    # df_1 = df_n.loc[df_n['3_y'] == 1, :]
    # df_merge_1 = df_merge.loc[df_merge['3_y'] == 1, [1, 8]]
    # df_merge_1_unique = df_merge_1.drop_duplicates()
