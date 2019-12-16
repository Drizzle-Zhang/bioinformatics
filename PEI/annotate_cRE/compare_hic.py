#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: compare_hic.py
# @time: 12/16/19 4:08 PM

from time import time
import re
import os
import pandas as pd


def transform_ng2019(file_pp, file_po, file_out, cutoff):
    with open(file_out, 'w') as w_out:
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
    file_merge_bin1 = file_out + '.merge.bin1'
    file_merge_bin2 = file_out + '.merge.bin2'
    os.system(f"bedtools merge -i {file_intersect_bin1} -c 4,5,6,7 "
              f"-o distinct,collapse,collapse,collapse > {file_merge_bin1}")
    os.system(f"bedtools merge -i {file_intersect_bin2} -c 4,5,6,7 "
              f"-o distinct,collapse,collapse,collapse > {file_merge_bin2}")
    file_res_bin1 = file_out + '.result.bin1'
    file_res_bin2 = file_out + '.result.bin2'
    os.system(f"bedtools intersect -a {file_bin1} -b {file_merge_bin1} -loj | "
              f"cut -f 1,2,3,4,9,10,11 > {file_res_bin1}")
    os.system(f"bedtools intersect -a {file_bin2} -b {file_merge_bin2} -loj | "
              f"cut -f 1,2,3,4,9,10,11 > {file_res_bin2}")
    df_bin1 = pd.read_csv(file_res_bin1, sep='\t', header=None)
    df_bin2 = pd.read_csv(file_res_bin2, sep='\t', header=None)
    df_out = pd.merge(df_bin1, df_bin2, on=3)
    df_out[7] = df_out[3]
    df_out = df_out.drop(3, axis=1)
    df_out.to_csv(file_out, sep='\t', header=None, index=None)

    os.remove(file_bin1)
    os.remove(file_bin2)
    os.remove(file_intersect_bin1)
    os.remove(file_intersect_bin2)
    os.remove(file_merge_bin1)
    os.remove(file_merge_bin2)
    os.remove(file_res_bin1)
    os.remove(file_res_bin2)

    return


def stat_effective_interaction(list_files):

    return


if __name__ == '__main__':
    time_start = time()
    # unify data format and cRE annotation
    file_cre_braincc = '/local/zy/PEI/data/DHS/combine_H3K27ac/brain/' \
                       'adult_cerebral_cortex/cRE.txt'
    # promoter capture Hi-C
    ng_pp = '/local/zy/PEI/compare_DLPFC/NGpcHiC/X5628FC.pp.txt'
    ng_po = '/local/zy/PEI/compare_DLPFC/NGpcHiC/X5628FC.po.txt'

    # ng_uniform = '/local/zy/PEI/compare_DLPFC/NGpcHiC/interactions_0.05.txt'
    # interaction_annotation = \
    #     '/local/zy/PEI/compare_DLPFC/NGpcHiC/interactions_0.05.cRE.txt'
    # transform_ng2019(ng_pp, ng_po, ng_uniform, 1.3)
    # annotate_hic(file_cre_braincc, ng_uniform, interaction_annotation)
    #
    # ng_uniform = '/local/zy/PEI/compare_DLPFC/NGpcHiC/interactions_0.01.txt'
    # interaction_annotation = \
    #     '/local/zy/PEI/compare_DLPFC/NGpcHiC/interactions_0.01.cRE.txt'
    # transform_ng2019(ng_pp, ng_po, ng_uniform, 2)
    # annotate_hic(file_cre_braincc, ng_uniform, interaction_annotation)
    #
    # ng_uniform = '/local/zy/PEI/compare_DLPFC/NGpcHiC/interactions_0.001.txt'
    # interaction_annotation = \
    #     '/local/zy/PEI/compare_DLPFC/NGpcHiC/interactions_0.001.cRE.txt'
    # transform_ng2019(ng_pp, ng_po, ng_uniform, 3)
    # annotate_hic(file_cre_braincc, ng_uniform, interaction_annotation)
    #
    # ng_uniform = '/local/zy/PEI/compare_DLPFC/NGpcHiC/interactions_0.0001.txt'
    # interaction_annotation = \
    #     '/local/zy/PEI/compare_DLPFC/NGpcHiC/interactions_0.0001.cRE.txt'
    # transform_ng2019(ng_pp, ng_po, ng_uniform, 4)
    # annotate_hic(file_cre_braincc, ng_uniform, interaction_annotation)

    ng_uniform = '/local/zy/PEI/compare_DLPFC/NGpcHiC/interactions_0.00001.txt'
    interaction_annotation = \
        '/local/zy/PEI/compare_DLPFC/NGpcHiC/interactions_0.00001.cRE.txt'
    transform_ng2019(ng_pp, ng_po, ng_uniform, 5)
    annotate_hic(file_cre_braincc, ng_uniform, interaction_annotation)

    # 3DIV
    file_origin = \
        '/local/zy/PEI/compare_DLPFC/3DIV/HiCaptureSeq_DL_cutoff_10.txt'

    # file_uniform = '/local/zy/PEI/compare_DLPFC/3DIV/interactions_20.txt'
    # interaction_annotation = \
    #     '/local/zy/PEI/compare_DLPFC/3DIV/interactions_20.cRE.txt'
    # transform_3div(file_origin, file_uniform, 20)
    # annotate_hic(file_cre_braincc, file_uniform, interaction_annotation)

    file_uniform = '/local/zy/PEI/compare_DLPFC/3DIV/interactions_30.txt'
    interaction_annotation = \
        '/local/zy/PEI/compare_DLPFC/3DIV/interactions_30.cRE.txt'
    transform_3div(file_origin, file_uniform, 30)
    annotate_hic(file_cre_braincc, file_uniform, interaction_annotation)

    # 3DGB
    file_origin = '/local/zy/PEI/compare_DLPFC/3DGB/' \
                  'Schmitt_2016.Cortex_DLPFC.hg19.peakachu-merged.loops'
    file_uniform = '/local/zy/PEI/compare_DLPFC/3DGB/interactions.txt'
    interaction_annotation = \
        '/local/zy/PEI/compare_DLPFC/3DGB/interactions.cRE.txt'
    transform_3dgb(file_origin, file_uniform)
    annotate_hic(file_cre_braincc, file_uniform, interaction_annotation)

    time_end = time()
    print(time_end - time_start)
