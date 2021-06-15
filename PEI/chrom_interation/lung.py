#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: lung.py
# @time: 5/24/20 2:38 PM

from time import time
import re
import os
import pandas as pd
import numpy as np
from multiprocessing import Pool


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
                if chrom1 == chrom2:
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
                if chrom1 == chrom2:
                    w_out.write(fmt.format(**locals()))

    os.system(f"sort -k 1,1 -k2,2n {file_tmp} > {file_out}")
    os.remove(file_tmp)

    return


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
                # if ('.' in set_cres1) | ('.' in set_cres2):
                #     continue
                dhs_ids1 = list_line[3].split(',')
                dhs_ids2 = list_line[9].split(',')
                promoters1 = set(list_line[5].split(','))
                promoters2 = set(list_line[11].split(','))
                if loop_id[:2] == 'pp':
                    genes1 = [pattern_gene.search(val).group()[:-2]
                              for val in promoters1 if val != '.']
                    genes2 = [pattern_gene.search(val).group()[:-2]
                              for val in promoters2 if val != '.']
                    if genes2:
                        for i, dhs_id in enumerate(dhs_ids1):
                            if dhs_id == '.':
                                continue
                            for gene2 in genes2:
                                dict_pro1 = dict(gene=gene2, dhs_id=dhs_id,
                                                 cre_type=cres1[i],
                                                 loop_score=loop_score)
                                w_out.write(fmt.format(**dict_pro1))
                    if genes1:
                        for i, dhs_id in enumerate(dhs_ids2):
                            if dhs_id == '.':
                                continue
                            for gene1 in genes1:
                                dict_pro2 = dict(gene=gene1, dhs_id=dhs_id,
                                                 cre_type=cres2[i],
                                                 loop_score=loop_score)
                                w_out.write(fmt.format(**dict_pro2))
                elif loop_id[:2] == 'po':
                    genes1 = [pattern_gene.search(val).group()[:-2]
                              for val in promoters1 if val != '.']
                    if genes1:
                        for i, dhs_id in enumerate(dhs_ids2):
                            if dhs_id == '.':
                                continue
                            for gene1 in genes1:
                                dict_pro2 = dict(gene=gene1, dhs_id=dhs_id,
                                                 cre_type=cres2[i],
                                                 loop_score=loop_score)
                                w_out.write(fmt.format(**dict_pro2))
                else:
                    genes1 = [pattern_gene.search(val).group()[:-2]
                              for val in promoters1 if val != '.']
                    genes2 = [pattern_gene.search(val).group()[:-2]
                              for val in promoters2 if val != '.']
                    if genes2:
                        for i, dhs_id in enumerate(dhs_ids1):
                            if dhs_id == '.':
                                continue
                            for gene2 in genes2:
                                dict_pro1 = dict(gene=gene2, dhs_id=dhs_id,
                                                 cre_type=cres1[i],
                                                 loop_score=loop_score)
                                w_out.write(fmt.format(**dict_pro1))
                    if genes1:
                        for i, dhs_id in enumerate(dhs_ids2):
                            if dhs_id == '.':
                                continue
                            for gene1 in genes1:
                                dict_pro2 = dict(gene=gene1, dhs_id=dhs_id,
                                                 cre_type=cres2[i],
                                                 loop_score=loop_score)
                                w_out.write(fmt.format(**dict_pro2))

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
              # f"grep -w 'Promoter\\|Enhancer' | "
              f"bedtools sort -i > {file_intersect_bin1}")
    os.system(f"bedtools intersect -a {file_bin2} -b {file_cre} -loj | "
              f"cut -f 1,2,3,4,5,9,10,12 | "
              # f"grep -w 'Promoter\\|Enhancer' | "
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

    file_pair_pre = file_pair + '.pre'
    # get_pairs(file_out, file_pair_pre)
    # os.system(f"grep -w 'Enhancer' {file_pair_pre} > {file_pair}")
    # os.remove(file_pair_pre)
    get_pairs(file_out, file_pair)

    return


def uniform_snp(file_in, file_out, file_cre, file_pair, flank=500):
    file_in_bed = file_out + '.bed'
    with open(file_in, 'r') as r_f:
        with open(file_in_bed, 'w') as w_f:
            fmt = "{chrom}\t{start}\t{end}\t{snp_id}\n"
            for line in r_f:
                tmp_line = line.strip().split('\t')
                chrom = 'chr' + tmp_line[0]
                start = int(tmp_line[1]) - flank
                end = int(tmp_line[1]) + flank
                snp_id = tmp_line[2]
                w_f.write(fmt.format(**locals()))

    file_snp_cre = file_out + '.snp.cRE'
    os.system(f"bedtools intersect -a {file_in_bed} -b {file_cre} -wao | "
              f"cut -f 4,8 > {file_snp_cre}")
    df_snp_cre = pd.read_csv(file_snp_cre, sep='\t', header=None,
                             names=['snp_id', 'dhs_id'])
    df_pair = pd.read_csv(file_pair, sep='\t', header=None,
                          names=['gene', 'dhs_id', 'type_cre', 'loop_score'])
    df_out = pd.merge(df_snp_cre, df_pair, on='dhs_id', how='left')
    df_out.to_csv(file_out, sep='\t', index=None)
    #
    # file_pro = '/local/zy/PEI/lung_snps/pro.txt'
    # file_enh = '/local/zy/PEI/lung_snps/enh.txt'
    # file_snp_hic = '/local/zy/PEI/lung_snps/snp_enh.txt'
    # os.system(f"cut -f 1,2,3,7,8 {file_uniform_lung} > {file_pro}")
    # os.system(f"cut -f 4,5,6,7,8 {file_uniform_lung} > {file_enh}")
    # os.system(f"bedtools intersect -a {file_in_bed} -b {file_enh} -wao | "
    #           f"cut -f 4,8 > {file_snp_hic}")

    return


def map_tissue(dict_in):
    file_pp = dict_in['file_pp']
    file_po = dict_in['file_po']
    file_cre = dict_in['file_cre']
    file_uniform = dict_in['file_uniform']
    file_out = dict_in['file_out']
    file_pair = dict_in['file_pair']
    tissue = dict_in['tissue']
    file_out_snp = dict_in['file_out_snp']
    transform_ng2019(file_pp, file_po, file_uniform, 1.3)
    annotate_hic(file_cre, file_uniform, file_out, file_pair)

    file_snp_cre = file_out_snp + '.snp.cRE'
    os.system(f"bedtools intersect -a {file_snp_uniform} -b {file_cre} -wao | "
              f"cut -f 4,8 > {file_snp_cre}")
    df_snp_cre = pd.read_csv(file_snp_cre, sep='\t', header=None,
                             names=['snp_id', 'dhs_id'])
    df_pair = pd.read_csv(file_pair, sep='\t', header=None,
                          names=['gene', 'dhs_id', 'type_cre', 'loop_score'])
    df_out = pd.merge(df_snp_cre, df_pair, on='dhs_id', how='left')
    df_out.to_csv(file_out_snp, sep='\t', index=None)
    df_out['tissue'] = [tissue for _ in range(df_out.shape[0])]

    return df_out


def map_label(dict_in):
    file_cre = dict_in['file_cre']
    file_uniform = dict_in['file_uniform']
    file_out = dict_in['file_out']
    file_pair = dict_in['file_pair']
    label = dict_in['label']
    file_out_snp = dict_in['file_out_snp']
    annotate_hic(file_cre, file_uniform, file_out, file_pair)

    file_snp_cre = file_out_snp + '.snp.cRE'
    os.system(f"bedtools intersect -a {file_snp_uniform} -b {file_cre} -wao | "
              f"cut -f 4,8 > {file_snp_cre}")
    df_snp_cre = pd.read_csv(file_snp_cre, sep='\t', header=None,
                             names=['snp_id', 'dhs_id'])
    df_pair = pd.read_csv(file_pair, sep='\t', header=None,
                          names=['gene', 'dhs_id', 'type_cre', 'loop_score'])
    df_out = pd.merge(df_snp_cre, df_pair, on='dhs_id', how='left')
    df_out.to_csv(file_out_snp, sep='\t', index=None)
    df_out['tissue'] = [label for _ in range(df_out.shape[0])]

    return df_out


def snp_map_to_multiple_tissue():
    with open(file_snp, 'r') as r_f:
        with open(file_snp_uniform, 'w') as w_f:
            fmt = "{chrom}\t{start}\t{end}\t{snp_id}\n"
            for line in r_f:
                tmp_line = line.strip().split('\t')
                chrom = 'chr' + tmp_line[0]
                start = int(tmp_line[1]) - flank
                end = int(tmp_line[1]) + flank
                snp_id = tmp_line[2]
                w_f.write(fmt.format(**locals()))

    df_meta = pd.read_csv(meta_file_pchic, sep='\t')
    df_meta = df_meta.dropna()
    list_input = []
    for sub_dict in df_meta.to_dict('records'):
        path_tissue = os.path.join(path_lung, sub_dict['tissue'])
        if not os.path.exists(path_tissue):
            os.mkdir(path_tissue)
        file_pp = os.path.join(path_hic, f"{sub_dict['tissue']}.pp.txt")
        file_po = os.path.join(path_hic, f"{sub_dict['tissue']}.po.txt")
        file_cre = os.path.join(sub_dict['file_cre'], 'cRE.txt')
        file_uniform = os.path.join(path_tissue, 'uniform.txt')
        file_out = os.path.join(path_tissue, 'interactions.cRE.txt')
        file_pair = os.path.join(path_tissue, 'pairs.gene.cRE.txt')
        file_out_snp = os.path.join(path_tissue, 'SNPs.gene.cRE.txt')
        list_input.append(
            {'file_pp': file_pp, 'file_po': file_po, 'file_cre': file_cre,
             'file_uniform': file_uniform, 'file_out': file_out,
             'file_pair': file_pair, 'tissue': sub_dict['tissue'],
             'file_out_snp': file_out_snp})

    # pool = Pool(processes=num_cpu)
    # list_out = pool.map(map_tissue, list_input)
    # pool.close()
    # df_out = pd.concat(list_out)
    # df_out.to_csv(
    #     os.path.join(path_lung, 'SNPs_gene_cRE.pcHiC.all.txt'),
    #     sep='\t', index=None)
    # df_out = df_out.loc[df_out['dhs_id'] != '.', :]
    # df_out.to_csv(
    #     os.path.join(path_lung, 'SNPs_gene_cRE.pcHiC.txt'),
    #     sep='\t', index=None)

    df_meta = pd.read_csv(flie_meta_label, sep='\t')
    list_input = []
    for sub_dict in df_meta.to_dict('records'):
        term = sub_dict['Biosample term name']
        method = sub_dict['Method']
        source = sub_dict['Source']
        filename = sub_dict['Filename']
        path_term = os.path.join(path_lung, term)
        if not os.path.exists(path_term):
            os.mkdir(path_term)
        file_cre = os.path.join(path_ref_cellline, f"{term}/cRE.txt")
        file_uniform = os.path.join(
            path_label, f'{term}/{method}__{source}__{filename[:-4]}.uniform')
        file_out = os.path.join(
            path_term,
            f'{method}__{source}__{filename[:-4]}.interactions.cRE.txt')
        file_pair = os.path.join(
            path_term,
            f'{method}__{source}__{filename[:-4]}.pairs.gene.cRE.txt')
        file_out_snp = os.path.join(
            path_term,
            f'{method}__{source}__{filename[:-4]}.SNPs.gene.cRE.txt')
        list_input.append(
            {'file_cre': file_cre,
             'label': f"{method}__{source}__{filename[:-4]}",
             'file_uniform': file_uniform, 'file_out': file_out,
             'file_pair': file_pair, 'file_out_snp': file_out_snp})

    pool = Pool(processes=num_cpu)
    list_out = pool.map(map_label, list_input)
    pool.close()
    df_out = pd.concat(list_out)
    df_out.to_csv(
        os.path.join(path_lung, 'SNPs_gene_cRE.label.all.txt'),
        sep='\t', index=None)
    df_out = df_out.loc[df_out['dhs_id'] != '.', :]
    df_out.to_csv(
        os.path.join(path_lung, 'SNPs_gene_cRE.label.txt'),
        sep='\t', index=None)

    return


if __name__ == '__main__':
    time_start = time()
    num_cpu = 40
    flank = 1000
    path_hic = \
        '/local/zy/PEI/origin_data/Chromatin_interactions/pcHi-C/Jung_NG_2019'
    path_lung = '/local/zy/PEI/lung_snps/multi_tissue'

    # example
    file_pp_lung = \
        '/local/zy/PEI/origin_data/Chromatin_interactions/pcHi-C/' \
        'Jung_NG_2019/LG.pp.txt'
    file_po_lung = \
        '/local/zy/PEI/origin_data/Chromatin_interactions/pcHi-C/' \
        'Jung_NG_2019/LG.po.txt'
    file_uniform_lung = '/local/zy/PEI/lung_snps/lung.uniform'
    file_cre_lung = \
        '/local/zy/PEI/mid_data/tissue/DHS/cRE_annotation/' \
        'adult_lung/lung/cRE.txt'
    file_out_lung = '/local/zy/PEI/lung_snps/lung.interactions.cRE.txt'
    file_pair_lung = '/local/zy/PEI/lung_snps/lung.pairs.gene.cRE.txt'
    # transform_ng2019(file_pp_lung, file_po_lung, file_uniform_lung, 1.3)
    # annotate_hic(
    #     file_cre_lung, file_uniform_lung, file_out_lung, file_pair_lung)

    file_snp = '/local/zy/PEI/lung_snps/top_SNPs.txt'
    file_snp_cre_gene = '/local/zy/PEI/lung_snps/SNPs_cRE_gene.txt'
    # uniform_snp(file_snp, file_snp_cre_gene, file_cre_lung, file_pair_lung,
    #             flank=1000)

    # multiple tissue
    file_snp_uniform = '/local/zy/PEI/lung_snps/top_SNPs_uniform.txt'
    meta_file_pchic = '/local/zy/PEI/origin_data/meta_file/meta_pcHi-C.txt'
    path_label = \
        '/local/zy/PEI/mid_data/training_label/label_interactions_V1'
    flie_meta_label = os.path.join(path_label, 'meta_label.txt')
    path_ref_cellline = '/local/zy/PEI/mid_data/cell_line/DHS/cRE_annotation'
    snp_map_to_multiple_tissue()

    time_end = time()
    print(time_end - time_start)
