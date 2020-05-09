#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: compare_hic.py
# @time: 12/16/19 4:08 PM

from time import time
import re
import os
import pandas as pd
from multiprocessing import Pool
import numpy as np
from subprocess import check_output
from functools import partial


def transform_ng2019(file_pp, file_po, file_out, cutoff):
    file_tmp = file_out + '.tmp'
    with open(file_tmp, 'w') as w_out:
        fmt = "{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t" \
              "{loop_id}\t{score}\n"
        pattern_chrom = re.compile(r'.+:')
        pattern_start = re.compile(r':.+-')
        pattern_end = re.compile(r'-.+')
        with open(file_pp, 'r') as r_pp:
            for i, line in enumerate(r_pp):
                list_line = line.strip().split('\t')
                if list_line[0] == 'type':
                    continue
                score = list_line[11]
                if float(score) < cutoff:
                    continue
                loop_id = 'pp' + str(i)
                chrom1 = pattern_chrom.search(list_line[1]).group()[:-1]
                start1 = pattern_start.search(list_line[1]).group()[1:-1]
                end1 = pattern_end.search(list_line[1]).group()[1:]
                chrom2 = pattern_chrom.search(list_line[2]).group()[:-1]
                start2 = pattern_start.search(list_line[2]).group()[1:-1]
                end2 = pattern_end.search(list_line[2]).group()[1:]
                w_out.write(fmt.format(**locals()))
        with open(file_po, 'r') as r_po:
            for i, line in enumerate(r_po):
                list_line = line.strip().split('\t')
                if list_line[0] == 'type':
                    continue
                score = list_line[11]
                if float(score) < cutoff:
                    continue
                loop_id = 'po' + str(i)
                chrom1 = pattern_chrom.search(list_line[1]).group()[:-1]
                start1 = pattern_start.search(list_line[1]).group()[1:-1]
                end1 = pattern_end.search(list_line[1]).group()[1:]
                chrom2 = pattern_chrom.search(list_line[2]).group()[:-1]
                start2 = pattern_start.search(list_line[2]).group()[1:-1]
                end2 = pattern_end.search(list_line[2]).group()[1:]
                w_out.write(fmt.format(**locals()))

    os.system(f"sort -k 1,1 -k2,2n {file_tmp} > {file_out}")
    os.remove(file_tmp)

    return


def transform_3div(file_in, file_out, cutoff):
    with open(file_out, 'w') as w_out:
        fmt = "{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t" \
              "{loop_id}\t{score}\n"
        with open(file_in, 'r') as r_pp:
            for i, line in enumerate(r_pp):
                list_line = line.strip().split('\t')
                if list_line[0] == 'idx':
                    continue
                score = list_line[4]
                if float(score) < cutoff:
                    continue
                loop_id = str(i)
                chrom1 = list_line[15]
                start1 = list_line[16]
                end1 = int(start1) + 5000
                chrom2 = chrom1
                start2 = list_line[17]
                end2 = int(start2) + 5000
                w_out.write(fmt.format(**locals()))

    return


def transform_3dgb(file_in, file_out):
    with open(file_out, 'w') as w_out:
        fmt = "{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t" \
              "{loop_id}\t{score}\n"
        with open(file_in, 'r') as r_pp:
            for i, line in enumerate(r_pp):
                list_line = line.strip().split('\t')
                score = '.'
                loop_id = 'po' + str(i)
                chrom1 = list_line[0]
                start1 = list_line[1]
                end1 = list_line[2]
                chrom2 = list_line[3]
                start2 = list_line[4]
                end2 = list_line[5]
                w_out.write(fmt.format(**locals()))

    return


def transform_psych(file_in, file_out):
    with open(file_out, 'w') as w_out:
        fmt = "{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t" \
              "{loop_id}\t{score}\n"
        with open(file_in, 'r') as r_pp:
            for i, line in enumerate(r_pp):
                list_line = line.strip().split('\t')
                if len(list_line) != 5:
                    continue
                score = '.'
                loop_id = str(i)
                chrom1 = 'chr' + list_line[0]
                start1 = list_line[1]
                end1 = list_line[2]
                chrom2 = 'chr' + list_line[0]
                start2 = list_line[3]
                end2 = list_line[4]
                w_out.write(fmt.format(**locals()))

    return


def transform_plac(file_in, file_out):
    with open(file_out, 'w') as w_out:
        fmt = "{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t" \
              "{loop_id}\t{score}\n"
        with open(file_in, 'r') as r_pp:
            for i, line in enumerate(r_pp):
                list_line = line.strip().split('\t')
                if len(list_line) > 6:
                    score = list_line[8]
                else:
                    score = '.'
                loop_id = str(i)
                chrom1 = list_line[0]
                start1 = list_line[1]
                end1 = list_line[2]
                chrom2 = list_line[3]
                start2 = list_line[4]
                end2 = list_line[5]
                w_out.write(fmt.format(**locals()))

    return


def annotate_hic(file_cre, file_hic, file_out):
    file_bin1 = file_out + '.bin1'
    file_bin2 = file_out + '.bin2'
    os.system(f"cut -f 1,2,3,7 {file_hic} > {file_bin1}")
    os.system(f"cut -f 4,5,6,7 {file_hic} > {file_bin2}")
    file_intersect_bin1 = file_out + '.intersect.bin1'
    file_intersect_bin2 = file_out + '.intersect.bin2'
    os.system(f"bedtools intersect -a {file_bin1} -b {file_cre} -loj | "
              f"cut -f 1,2,3,4,8,9,10 | "
              f"grep -w 'Promoter\\|Enhancer' | "
              f"bedtools sort -i > {file_intersect_bin1}")
    os.system(f"bedtools intersect -a {file_bin2} -b {file_cre} -loj | "
              f"cut -f 1,2,3,4,8,9,10 | "
              f"grep -w 'Promoter\\|Enhancer' | "
              f"bedtools sort -i > {file_intersect_bin2}")

    def merge_bin(df_in):
        dhs_ids = ','.join(df_in[4].tolist())
        cres = ','.join(df_in[5].tolist())
        promoters = ','.join(df_in[6].tolist())
        dict_out = {0: df_in.iloc[0, 0], 1: df_in.iloc[0, 1],
                    2: df_in.iloc[0, 2], 3: df_in.iloc[0, 3],
                    4: dhs_ids, 5: cres, 6: promoters}
        df_dict_out = pd.Series(dict_out)

        return df_dict_out

    df_intersect_bin1 = pd.read_csv(file_intersect_bin1, sep='\t', header=None)
    df_intersect_bin1 = df_intersect_bin1.rename(columns={3: 'key'})
    df_merge_bin1 = df_intersect_bin1.groupby('key').apply(merge_bin)
    df_intersect_bin2 = pd.read_csv(file_intersect_bin2, sep='\t', header=None)
    df_intersect_bin2 = df_intersect_bin2.rename(columns={3: 'key'})
    df_merge_bin2 = df_intersect_bin2.groupby('key').apply(merge_bin)
    df_bin1 = pd.read_csv(file_bin1, sep='\t', header=None)
    df_bin2 = pd.read_csv(file_bin2, sep='\t', header=None)
    df_res_bin1 = pd.merge(df_bin1, df_merge_bin1,
                           on=[0, 1, 2, 3], how='outer')
    df_res_bin2 = pd.merge(df_bin2, df_merge_bin2,
                           on=[0, 1, 2, 3], how='outer')
    df_out = pd.merge(df_res_bin1, df_res_bin2, on=3)
    df_out[7] = df_out[3]
    df_out = df_out.drop(3, axis=1)
    df_out.to_csv(file_out, sep='\t', header=None, index=None, na_rep='.')

    os.remove(file_bin1)
    os.remove(file_bin2)
    os.remove(file_intersect_bin1)
    os.remove(file_intersect_bin2)

    return


def sub_stat(egenes, eqtl_file, dict_in):
    file_in = dict_in['file']
    file_out = dict_in['file_out']
    file_egene = dict_in['file_egene']
    file_overlap = dict_in['file_overlap']
    df_file = pd.read_csv(file_in, sep='\t', header=None)
    resolution = (np.mean(df_file[2] - df_file[1]) +
                  np.mean(df_file[8] - df_file[7]))/2
    len_file = df_file.shape[0]
    count_effect = 0
    pattern_gene = re.compile(r'.+<-')
    with open(file_out, 'w') as w_out:
        fmt = "{gene}\t{dhs_id}\t{cre_type}\n"
        with open(file_in, 'r') as r_f:
            for line in r_f:
                list_line = line.strip().split('\t')
                line_id = list_line[-1]
                cres1 = list_line[4].split(',')
                cres2 = list_line[10].split(',')
                set_cres1 = set(list_line[4].split(','))
                set_cres2 = set(list_line[10].split(','))
                count_add = 0
                if line_id[0:2] == 'pp':
                    if ('Promoter' in set_cres1) & \
                            (('Promoter' in set_cres2) |
                             ('Enhancer' in set_cres2)):
                        count_add = 1
                elif line_id[0:2] == 'po':
                    if ('Promoter' in set_cres1) & \
                            (('Promoter' in set_cres2) |
                             ('Enhancer' in set_cres2)):
                        count_add = 1
                else:
                    if ('.' in set_cres1) | ('.' in set_cres2):
                        continue
                    else:
                        count_add = 1
                count_effect = count_effect + count_add
                if count_add == 1:
                    dhs_ids1 = list_line[3].split(',')
                    dhs_ids2 = list_line[9].split(',')
                    promoters1 = set(list_line[5].split(','))
                    promoters2 = set(list_line[11].split(','))
                    if ('Promoter' in set_cres1) & ('Promoter' in set_cres2):
                        genes1 = ','.join(
                            [pattern_gene.search(val).group()[:-2]
                             for val in promoters1 if val != '.'])
                        genes2 = ','.join(
                            [pattern_gene.search(val).group()[:-2]
                             for val in promoters2 if val != '.'])
                        for i, dhs_id in enumerate(dhs_ids1):
                            dict_pro1 = dict(gene=genes1, dhs_id=dhs_id,
                                             cre_type=cres1[i])
                            w_out.write(fmt.format(**dict_pro1))
                        for i, dhs_id in enumerate(dhs_ids2):
                            dict_pro2 = dict(gene=genes2, dhs_id=dhs_id,
                                             cre_type=cres2[i])
                            w_out.write(fmt.format(**dict_pro2))
                    elif ('Promoter' in set_cres1) & ('Enhancer' in set_cres2):
                        genes1 = ','.join(
                            [pattern_gene.search(val).group()[:-2]
                             for val in promoters1 if val != '.'])
                        for i, dhs_id in enumerate(dhs_ids2):
                            dict_pro2 = dict(gene=genes1, dhs_id=dhs_id,
                                             cre_type=cres2[i])
                            w_out.write(fmt.format(**dict_pro2))
                    elif ('Enhancer' in set_cres1) & ('Promoter' in set_cres2):
                        genes2 = ','.join(
                            [pattern_gene.search(val).group()[:-2]
                             for val in promoters2 if val != '.'])
                        for i, dhs_id in enumerate(dhs_ids1):
                            dict_pro1 = dict(gene=genes2, dhs_id=dhs_id,
                                             cre_type=cres1[i])
                            w_out.write(fmt.format(**dict_pro1))

    # drop duplicates
    df_pair = pd.read_csv(file_out, sep='\t', header=None)
    df_pair = df_pair.drop_duplicates()
    df_pair.to_csv(file_out, sep='\t', header=None, index=None)

    len_cre_pairs = int(str(check_output(
        f"wc -l {file_out}", shell=True).strip()).split(' ')[0][2:])

    with open(file_egene, 'w') as w_egene:
        fmt = "{gene}\t{dhs_id}\n"
        with open(file_out, 'r') as r_out:
            for line in r_out:
                list_line = line.strip().split('\t')
                genes = list_line[0].split(',')
                for gene in genes:
                    if gene in egenes:
                        dict_out = dict(gene=gene, dhs_id=list_line[1])
                        w_egene.write(fmt.format(**dict_out))

    df_egene_pair = pd.read_csv(file_egene, sep='\t', header=None)
    df_egene_pair = df_egene_pair.drop_duplicates()
    df_egene_pair.to_csv(file_egene, sep='\t', header=None, index=None)
    num_egene = df_egene_pair[0].drop_duplicates().shape[0]
    df_eqtl = pd.read_csv(eqtl_file, sep='\t', header=None)
    df_overlap = pd.merge(df_egene_pair, df_eqtl, on=[0, 1])
    df_overlap.to_csv(file_overlap, sep='\t', header=None, index=None)
    num_overlap_egene = df_overlap[0].drop_duplicates().shape[0]
    # ratio_overlap =
    return {'label': dict_in['label'], 'count': count_effect,
            'total': len_file, 'ratio': count_effect/len_file,
            'resolution': resolution, 'num_pairs': len_cre_pairs,
            'num_egene_pairs': df_egene_pair.shape[0],
            'overlap_num': df_overlap.shape[0], 'num_egene': num_egene,
            'num_overlap_egene': num_overlap_egene}


def stat_effective_interaction(files, labels, out_files, egenes, eqtl_file,
                               egene_files, overlap_files, file_res):
    list_input = []
    for i, file in enumerate(files):
        list_input.append(dict(file=file, label=labels[i],
                               file_out=out_files[i],
                               file_egene=egene_files[i],
                               file_overlap=overlap_files[i]))

    pool = Pool(processes=10)
    func_stat = partial(sub_stat, egenes, eqtl_file)
    list_dict = pool.map(func_stat, list_input)
    pool.close()

    # merge all files
    list_egene_files = []
    for file in egene_files:
        df_sub = pd.read_csv(file, sep='\t', header=None)
        list_egene_files.append(df_sub)

    df_egene_pair = (pd.concat(list_egene_files)).drop_duplicates()
    num_egene = df_egene_pair[0].drop_duplicates().shape[0]
    df_eqtl = pd.read_csv(eqtl_file, sep='\t', header=None)
    df_overlap = pd.merge(df_egene_pair, df_eqtl, on=[0, 1])
    # df_overlap.to_csv(file_overlap, sep='\t', header=None, index=None)
    num_overlap_egene = df_overlap[0].drop_duplicates().shape[0]
    # ratio_overlap =
    list_dict.append(
        {'label': 'All', 'count': 0, 'total': 0, 'ratio': 0, 'resolution': 0,
         'num_pairs': 0, 'num_egene_pairs': df_egene_pair.shape[0],
         'overlap_num': df_overlap.shape[0], 'num_egene': num_egene,
         'num_overlap_egene': num_overlap_egene})

    df_out = pd.DataFrame(list_dict)
    df_out.to_csv(file_res, sep='\t', index=None)

    return


def transform_eqtl(file_in, file_out, file_pair_eqtl, file_pair_cre,
                   file_cre, df_egene, cutoff=5000):
    pattern_start = re.compile(r'_[0-9]+_')
    pattern_chrom = re.compile(r'.+?_')
    with open(file_out, 'w') as w_out:
        fmt = "{chrom1}\t{start1}\t{end1}\t{variant_id}\t" \
              "{gene_symbol}\t{ensg}\t{tss_distance}\t{pval_beta}\n"
        with open(file_in, 'r') as r_in:
            for line in r_in:
                list_line = line.strip().split('\t')
                if list_line[0] == 'variant_id':
                    continue
                if abs(int(list_line[2])) < cutoff:
                    continue
                ensg = list_line[1]
                start = int(pattern_start.search(list_line[0]).group()[1:-1])
                df_ensg = df_egene.loc[df_egene['gene_id'] == ensg, :]
                if df_ensg.shape[0] != 1:
                    continue
                gene_start = df_ensg['gene_start'].iloc[0]
                gene_end = df_ensg['gene_end'].iloc[0]
                if (start >= gene_start) & (start <= gene_end):
                    continue
                chrom1 = f"chr{pattern_chrom.match(list_line[0]).group()[:-1]}"
                variant_id = list_line[0]
                start1 = str(start - 500)
                end1 = str(start + 500)
                tss_distance = list_line[2]
                pval_beta = list_line[-1]
                gene_symbol = df_ensg['gene_name'].iloc[0]
                w_out.write(fmt.format(**locals()))

    os.system(f"bedtools intersect -a {file_out} -b {file_cre} -loj | "
              f"cut -f 4,5,6,12,13 | "
              f"grep -w 'Promoter\\|Enhancer' > {file_pair_eqtl}")
    df_eqtl = pd.read_csv(file_pair_eqtl, sep='\t', header=None)
    df_eqtl = df_eqtl.drop_duplicates()
    df_eqtl.to_csv(file_pair_eqtl, sep='\t', header=None, index=None)
    df_eqtl = df_eqtl.iloc[:, [1, 3]]
    df_eqtl = df_eqtl.drop_duplicates()
    df_eqtl.to_csv(file_pair_cre, sep='\t', header=None, index=None)

    return


def overlap_matrix(list_label, list_file, mat_out):
    list_vec = []
    for i in range(len(list_label)):
        dict_vec = {}
        df_row = pd.read_csv(list_file[i], sep='\t', header=None)
        for j, label in enumerate(list_label):
            df_col = pd.read_csv(list_file[j], sep='\t', header=None)
            df_overlap = pd.merge(df_row, df_col, on=[0, 1], how='inner')
            dict_vec[label] = \
                df_overlap.shape[0] / min(df_row.shape[0], df_col.shape[0])
        list_vec.append(dict_vec)

    df_mat = pd.DataFrame(list_vec, columns=list_label)
    df_mat.index = list_label
    df_mat.to_csv(mat_out, sep='\t')

    return


def sub_cutoff_ng(file_pp, file_po, file_cre, egenes, eqtl_file,
                  path_out, cutoff):
    uniform = os.path.join(path_out, 'uniform_' + str(cutoff) + '.txt')
    transform_ng2019(file_pp, file_po, uniform, cutoff)
    file_annotation = os.path.join(
        path_out, 'interactions_cRE' + str(cutoff) + '.txt')
    annotate_hic(file_cre, uniform, file_annotation)
    file_pair = os.path.join(path_out, 'Gene_cRE_' + str(cutoff) + '.txt')
    file_egene_pair = os.path.join(
        path_out, 'eGene_cRE_' + str(cutoff) + '.txt')
    file_overlap = os.path.join(
        path_out, 'eGene_cRE_overlap_' + str(cutoff) + '.txt')
    dict_in = dict(file=file_annotation, label=str(cutoff),
                   file_out=file_pair, file_egene=file_egene_pair,
                   file_overlap=file_overlap)
    dict_out = sub_stat(egenes, eqtl_file, dict_in)

    return dict_out


def sub_cutoff_3div(file_in, file_cre, egenes, eqtl_file,
                    path_out, cutoff):
    uniform = os.path.join(path_out, 'uniform_' + str(cutoff) + '.txt')
    transform_3div(file_in, uniform, cutoff)
    file_annotation = os.path.join(
        path_out, 'interactions_cRE' + str(cutoff) + '.txt')
    annotate_hic(file_cre, uniform, file_annotation)
    file_pair = os.path.join(path_out, 'Gene_cRE_' + str(cutoff) + '.txt')
    file_egene_pair = os.path.join(
        path_out, 'eGene_cRE_' + str(cutoff) + '.txt')
    file_overlap = os.path.join(
        path_out, 'eGene_cRE_overlap_' + str(cutoff) + '.txt')
    dict_in = dict(file=file_annotation, label=str(cutoff),
                   file_out=file_pair, file_egene=file_egene_pair,
                   file_overlap=file_overlap)
    dict_out = sub_stat(egenes, eqtl_file, dict_in)

    return dict_out


if __name__ == '__main__':
    time_start = time()
    # unify data format and cRE annotation
    file_cre_braincc = '/local/zy/PEI/data/DHS/cRE_annotation/' \
                       'brain/adult_cerebral_cortex/cRE.txt '
    # promoter capture Hi-C
    ng_pp = '/local/zy/PEI/compare_DLPFC/NGpcHiC/X5628FC.pp.txt'
    ng_po = '/local/zy/PEI/compare_DLPFC/NGpcHiC/X5628FC.po.txt'

    ng_uniform = '/local/zy/PEI/compare_DLPFC/NGpcHiC/interactions_0.01.txt'
    interaction_annotation_ng = \
        '/local/zy/PEI/compare_DLPFC/NGpcHiC/interactions_0.01.cRE.txt'
    transform_ng2019(ng_pp, ng_po, ng_uniform, 2)
    annotate_hic(file_cre_braincc, ng_uniform, interaction_annotation_ng)

    # 3DIV
    file_origin_3div = \
        '/local/zy/PEI/compare_DLPFC/3DIV/HiCaptureSeq_DL_cutoff_10.txt'

    uniform_3div = '/local/zy/PEI/compare_DLPFC/3DIV/interactions_14.txt'
    interaction_annotation_3div = \
        '/local/zy/PEI/compare_DLPFC/3DIV/interactions_14.cRE.txt'
    transform_3div(file_origin_3div, uniform_3div, 14)
    annotate_hic(file_cre_braincc, uniform_3div, interaction_annotation_3div)

    # 3DGB
    file_origin_3dgb = '/local/zy/PEI/compare_DLPFC/3DGB/' \
                       'Schmitt_2016.Cortex_DLPFC.hg19.peakachu-merged.loops'
    uniform_3dgb = '/local/zy/PEI/compare_DLPFC/3DGB/interactions.txt'
    interaction_annotation_3dgb = \
        '/local/zy/PEI/compare_DLPFC/3DGB/interactions.cRE.txt'
    transform_3dgb(file_origin_3dgb, uniform_3dgb)
    annotate_hic(file_cre_braincc, uniform_3dgb, interaction_annotation_3dgb)

    # psychENCODE
    file_origin_psych = '/local/zy/PEI/compare_DLPFC/psychENCODE/' \
                        'Promoter-anchored_chromatin_loops.bed'
    uniform_psych = '/local/zy/PEI/compare_DLPFC/psychENCODE/interactions.txt'
    interaction_annotation_psych = \
        '/local/zy/PEI/compare_DLPFC/psychENCODE/interactions.cRE.txt'
    transform_psych(file_origin_psych, uniform_psych)
    annotate_hic(file_cre_braincc, uniform_psych, interaction_annotation_psych)

    # PLAC-seq microglia
    file_origin_microglia = \
        '/local/zy/PEI/compare_DLPFC/microglia/microglia_interactions.txt'
    uniform_microglia = \
        '/local/zy/PEI/compare_DLPFC/microglia/interactions.txt'
    interaction_annotation_microglia = \
        '/local/zy/PEI/compare_DLPFC/microglia/interactions.cRE.txt'
    transform_plac(file_origin_microglia, uniform_microglia)
    annotate_hic(file_cre_braincc, uniform_microglia,
                 interaction_annotation_microglia)

    # PLAC-seq neuron
    file_origin_neuron = \
        '/local/zy/PEI/compare_DLPFC/neuron/Neuronal_interactions.txt'
    uniform_neuron = \
        '/local/zy/PEI/compare_DLPFC/neuron/interactions.txt'
    interaction_annotation_neuron = \
        '/local/zy/PEI/compare_DLPFC/neuron/interactions.cRE.txt'
    transform_plac(file_origin_neuron, uniform_neuron)
    annotate_hic(file_cre_braincc, uniform_neuron,
                 interaction_annotation_neuron)

    # PLAC-seq oligodendrocyte
    file_origin_oligodendrocyte = \
        '/local/zy/PEI/compare_DLPFC/oligodendrocyte/' \
        'oligodendrocyte_interactions.txt'
    uniform_oligodendrocyte = \
        '/local/zy/PEI/compare_DLPFC/oligodendrocyte/interactions.txt'
    interaction_annotation_oligodendrocyte = \
        '/local/zy/PEI/compare_DLPFC/oligodendrocyte/interactions.cRE.txt'
    transform_plac(file_origin_oligodendrocyte, uniform_oligodendrocyte)
    annotate_hic(file_cre_braincc, uniform_oligodendrocyte,
                 interaction_annotation_oligodendrocyte)

    # PLAC-seq
    file_origin_plac = \
        '/local/zy/PEI/compare_DLPFC/PLACseq/PLACseq.uniq.txt'
    uniform_plac = \
        '/local/zy/PEI/compare_DLPFC/PLACseq/interactions.txt'
    interaction_annotation_plac = \
        '/local/zy/PEI/compare_DLPFC/PLACseq/interactions.cRE.txt'
    transform_plac(file_origin_plac, uniform_plac)
    annotate_hic(file_cre_braincc, uniform_plac, interaction_annotation_plac)

    # eQTL
    # eGenes
    file_egenes = '/local/zy/PEI/compare_DLPFC/eQTL/Cortex/' \
                  'Brain_Cortex.v7.egenes.txt'
    protein_gene = '/local/zy/PEI/data/gene/genes.protein.gencode.v19.bed'
    df_protein_genes = pd.read_csv(protein_gene, sep='\t', header=None)
    df_protein_genes.columns = ['chrom', 'start', 'end',
                                'gene_name', 'gene_id', 'strand']
    df_egenes = pd.read_csv(file_egenes, sep='\t')
    df_egenes = df_egenes.loc[
        df_egenes['qval'] <= 0.05,
        ['gene_id', 'gene_name', 'gene_start', 'gene_end']
    ]
    df_egenes_pos = pd.merge(
        df_protein_genes, df_egenes, on=['gene_id', 'gene_name']
    )
    df_egenes_pos = df_egenes_pos.loc[
                    :, ['chrom', 'gene_start', 'gene_end', 'gene_name',
                        'gene_id', 'strand']
                    ]

    file_eqtl = '/local/zy/PEI/compare_DLPFC/eQTL/Cortex/' \
                'Brain_Cortex.v7.signif_variant_gene_pairs.txt'
    file_uniform = '/local/zy/PEI/compare_DLPFC/eQTL/Cortex/' \
                   'uniform_500.txt'
    pairs_eqtl_egene = \
        '/local/zy/PEI/compare_DLPFC/eQTL/Cortex/' \
        'eQTL_cRE_pairs_500.txt'
    pairs_dhs_egene = '/local/zy/PEI/compare_DLPFC/eQTL/Cortex/' \
                      'cRE_pairs_500.txt'
    transform_eqtl(
        file_eqtl, file_uniform, pairs_eqtl_egene, pairs_dhs_egene,
        file_cre_braincc, df_egenes_pos)

    # effective interactions and get gene-cre pairs
    list_files = [interaction_annotation_ng, interaction_annotation_3div,
                  interaction_annotation_3dgb, interaction_annotation_psych,
                  interaction_annotation_microglia,
                  interaction_annotation_neuron,
                  interaction_annotation_oligodendrocyte,
                  interaction_annotation_plac]
    list_labels = ['pcHi-C', '3DIV', '3DGB', 'psychENCODE', 'microglia',
                   'neuron', 'oligodendrocyte', 'PLAC-seq']
    pairs_ng = \
        '/local/zy/PEI/compare_DLPFC/NGpcHiC/cRE_pairs_0.01.txt'
    pairs_3div = \
        '/local/zy/PEI/compare_DLPFC/3DIV/cRE_pairs_14.txt'
    pairs_3dgb = \
        '/local/zy/PEI/compare_DLPFC/3DGB/cRE_pairs.txt'
    pairs_psych = \
        '/local/zy/PEI/compare_DLPFC/psychENCODE/cRE_pairs.txt'
    pairs_microglia = \
        '/local/zy/PEI/compare_DLPFC/microglia/cRE_pairs.txt'
    pairs_neuron = \
        '/local/zy/PEI/compare_DLPFC/neuron/cRE_pairs.txt'
    pairs_oligodendrocyte = \
        '/local/zy/PEI/compare_DLPFC/oligodendrocyte/cRE_pairs.txt'
    pairs_plac = \
        '/local/zy/PEI/compare_DLPFC/PLACseq/cRE_pairs.txt'
    list_out = [pairs_ng, pairs_3div, pairs_3dgb, pairs_psych, pairs_microglia,
                pairs_neuron, pairs_oligodendrocyte, pairs_plac]
    egene_pairs_ng = \
        '/local/zy/PEI/compare_DLPFC/NGpcHiC/egene_cRE_pairs_0.05.txt'
    egene_pairs_3div = \
        '/local/zy/PEI/compare_DLPFC/3DIV/egene_cRE_pairs_11.txt'
    egene_pairs_3dgb = \
        '/local/zy/PEI/compare_DLPFC/3DGB/egene_cRE_pairs.txt'
    egene_pairs_psych = \
        '/local/zy/PEI/compare_DLPFC/psychENCODE/egene_cRE_pairs.txt'
    egene_pairs_microglia = \
        '/local/zy/PEI/compare_DLPFC/microglia/egene_cRE_pairs.txt'
    egene_pairs_neuron = \
        '/local/zy/PEI/compare_DLPFC/neuron/egene_cRE_pairs.txt'
    egene_pairs_oligodendrocyte = \
        '/local/zy/PEI/compare_DLPFC/oligodendrocyte/egene_cRE_pairs.txt'
    egene_pairs_plac = \
        '/local/zy/PEI/compare_DLPFC/PLACseq/egene_cRE_pairs.txt'
    list_egene = [egene_pairs_ng, egene_pairs_3div, egene_pairs_3dgb,
                  egene_pairs_psych, egene_pairs_microglia, egene_pairs_neuron,
                  egene_pairs_oligodendrocyte, egene_pairs_plac]
    overlap_pairs_ng = \
        '/local/zy/PEI/compare_DLPFC/NGpcHiC/' \
        'overlap_egene_cRE_pairs_0.05.txt'
    overlap_pairs_3div = \
        '/local/zy/PEI/compare_DLPFC/3DIV/overlap_egene_cRE_pairs_11.txt'
    overlap_pairs_3dgb = \
        '/local/zy/PEI/compare_DLPFC/3DGB/overlap_egene_cRE_pairs.txt'
    overlap_pairs_psych = \
        '/local/zy/PEI/compare_DLPFC/psychENCODE/overlap_egene_cRE_pairs.txt'
    overlap_pairs_microglia = \
        '/local/zy/PEI/compare_DLPFC/microglia/overlap_egene_cRE_pairs.txt'
    overlap_pairs_neuron = \
        '/local/zy/PEI/compare_DLPFC/neuron/overlap_egene_cRE_pairs.txt'
    overlap_pairs_oligodendrocyte = \
        '/local/zy/PEI/compare_DLPFC/oligodendrocyte/' \
        'overlap_egene_cRE_pairs.txt'
    overlap_pairs_plac = \
        '/local/zy/PEI/compare_DLPFC/PLACseq/overlap_egene_cRE_pairs.txt'
    list_overlap = [overlap_pairs_ng, overlap_pairs_3div, overlap_pairs_3dgb,
                    overlap_pairs_psych, overlap_pairs_microglia,
                    overlap_pairs_neuron, overlap_pairs_oligodendrocyte,
                    overlap_pairs_plac]
    set_egenes = set(df_egenes_pos['gene_name'].tolist())
    file_compare = '/local/zy/PEI/compare_DLPFC/compare_results.txt'
    stat_effective_interaction(
        list_files, list_labels, list_out, set_egenes, pairs_dhs_egene,
        list_egene, list_overlap, file_compare)

    # overlap of Hi-C from different source
    file_mat = '/local/zy/PEI/compare_DLPFC/overlap.mtx'
    overlap_matrix(list_labels + ['eQTL'], list_egene + [pairs_dhs_egene],
                   file_mat)

    # different cutoff
    range_cutoff_ng = [1, 1.3, 1.6, 2, 2.5, 3, 3.5, 4, 5]
    path_cutoff_ng = '/local/zy/PEI/compare_DLPFC/NGpcHiC/cutoff'
    pool = Pool(processes=10)
    func_cutoff = partial(
        sub_cutoff_ng, ng_pp, ng_po, file_cre_braincc, set_egenes,
        pairs_eqtl, path_cutoff_ng)
    list_dict = pool.map(func_cutoff, range_cutoff_ng)
    pool.close()
    df_ng = pd.DataFrame(list_dict)
    df_ng.to_csv(os.path.join(path_cutoff_ng, 'result.txt'),
                 sep='\t', index=None)

    range_cutoff_3div = [11, 12, 14, 16, 18, 20, 23, 26, 30]
    path_cutoff_3div = '/local/zy/PEI/compare_DLPFC/3DIV/cutoff'
    pool = Pool(processes=10)
    func_cutoff = partial(
        sub_cutoff_3div, file_origin_3div, file_cre_braincc, set_egenes,
        pairs_eqtl, path_cutoff_3div)
    list_dict = pool.map(func_cutoff, range_cutoff_3div)
    pool.close()
    df_3div = pd.DataFrame(list_dict)
    df_3div.to_csv(os.path.join(path_cutoff_3div, 'result.txt'),
                   sep='\t', index=None)

    time_end = time()
    print(time_end - time_start)
