#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: annotate_DHS.py
# @time: 10/28/19 4:08 PM

from time import time
import os
from multiprocessing import Pool
from functools import partial
import pandas as pd
from prepare_bed_file import merge_bed
from subprocess import check_output


def generate_file_list(path_in, path_out):
    folder_1 = os.listdir(path_in)
    list_input = []
    for element_1 in folder_1:
        if element_1[:4] == 'meta':
            continue
        path_1 = os.path.join(path_in, element_1)
        if not os.path.isdir(path_1):
            list_input.append(dict(path_out=path_out,
                                   file=element_1,
                                   path_in=path_in,
                                   organ='.',
                                   life_stage='.',
                                   term='.'))
        else:
            path_1_out = os.path.join(path_out, element_1)
            if not os.path.exists(path_1_out):
                os.makedirs(path_1_out)
            folder_2 = os.listdir(path_1)
            for element_2 in folder_2:
                path_2 = os.path.join(path_1, element_2)
                if not os.path.isdir(path_2):
                    list_input.append(dict(path_out=path_1_out,
                                           file=element_2,
                                           path_in=path_1,
                                           organ=element_1,
                                           life_stage='.',
                                           term='.'))
                else:
                    path_2_out = os.path.join(path_1_out, element_2)
                    if not os.path.exists(path_2_out):
                        os.makedirs(path_2_out)
                    folder_3 = os.listdir(path_2)
                    for element_3 in folder_3:
                        path_3 = os.path.join(path_2, element_3)
                        if not os.path.isdir(path_3):
                            list_input.append(dict(path_out=path_2_out,
                                                   file=element_3,
                                                   path_in=path_2,
                                                   organ=element_1,
                                                   life_stage=element_2,
                                                   term='.'))
                        else:
                            path_3_out = os.path.join(path_2_out, element_3)
                            if not os.path.exists(path_3_out):
                                os.makedirs(path_3_out)
                            folder_4 = os.listdir(path_3)
                            for element_4 in folder_4:
                                path_4 = os.path.join(path_3, element_4)
                                if not os.path.isdir(path_4):
                                    list_input.append(dict(path_out=path_3_out,
                                                           file=element_4,
                                                           path_in=path_3,
                                                           organ=element_1,
                                                           life_stage=
                                                           element_2,
                                                           term=element_3))

    return list_input


def sub_stan(type_bed, col_score, dict_in):
    path_out = dict_in['path_out']
    path_in = dict_in['path_in']
    file = dict_in['file']
    organ = dict_in['organ']
    life_stage = dict_in['life_stage']
    term = dict_in['term']
    file_label = f"{organ}|{life_stage}|{term}"
    file_in = os.path.join(path_in, file)
    with open(file_in, 'r') as r_f:
        with open(os.path.join(path_out, file), 'w') as w_f:
            fmt_dhs = "{chrom}\t{start}\t{end}\t{label}\t{file_label}\n"
            fmt_histone = "{chrom}\t{start}\t{end}\t{label}\t" \
                          "{score}\t{file_label}\n"
            for line in r_f:
                list_line = line.strip().split('\t')
                chrom = list_line[0]
                start = list_line[1]
                end = list_line[2]
                label = f"{type_bed}<-{chrom}:{start}-{end}"
                if type_bed == 'DHS':
                    w_f.write(fmt_dhs.format(**locals()))
                else:
                    score = max(list_line[col_score].strip().split(','))
                    w_f.write(fmt_histone.format(**locals()))

    return


def standardize_bed(path_in, path_out, type_bed):
    os.system(f"rm -rf {path_out}")
    os.mkdir(path_out)

    list_input = generate_file_list(path_in, path_out)
    pool = Pool(processes=40)
    func_stan = partial(sub_stan, type_bed, 8)
    pool.map(func_stan, list_input)
    pool.close()

    return


def split_ref_bed(ref_file, dict_in):
    organ = dict_in['organ']
    life_stage = dict_in['life_stage']
    term = dict_in['term']
    label = [val for val in [organ, life_stage, term] if val != '.']
    len_label = len(label)
    label = '|'.join(label)
    with open(ref_file, 'r') as r_ref:
        with open(os.path.join(dict_in['path_out'], dict_in['file']), 'w') \
                as w_f:
            fmt_dhs = "{chrom}\t{start}\t{end}\t{dhs_id}\t{label}\n"
            for line in r_ref:
                list_line = line.strip().split('\t')
                chrom = list_line[0]
                start = list_line[1]
                end = list_line[2]
                dhs_id = list_line[3]
                list_label = ['|'.join(val.strip().split('|')[:len_label])
                              for val in list_line[4].strip().split(',')]
                if label in list_label:
                    w_f.write(fmt_dhs.format(**locals()))

    return


def merge_split_bed(path_in, path_out):
    os.system(f"rm -rf {path_out}")
    os.mkdir(path_out)
    os.system(f"cp {os.path.join(path_in, 'metadata.tsv')} "
              f"{os.path.join(path_out, 'metadata.tsv')}")

    list_input = generate_file_list(path_in, path_out)
    df_list = pd.DataFrame(list_input)
    # merge
    df_subs = df_list.loc[df_list['life_stage'] != '.',
                          ['organ', 'life_stage', 'file']]
    dict_merge = dict(term_name='Reference_DHS', path=path_out,
                      accession_ids=[line['organ'] + '/' + line['life_stage']
                                     + '/' + line['file'][:-4]
                                     for line in df_subs.to_dict('records')])
    merge_bed(path_in, '5', dict_merge)

    # add uniform label
    all_ref = os.path.join(path_out, 'Reference_DHS.plus.bed')
    with open(os.path.join(path_out, 'Reference_DHS.bed'), 'r') as r_f:
        with open(all_ref, 'w') \
                as w_f:
            fmt_dhs = "{chrom}\t{start}\t{end}\t{label}\t{file_label}\n"
            for line in r_f:
                list_line = line.strip().split('\t')
                chrom = list_line[0]
                start = list_line[1]
                end = list_line[2]
                label = f"DHS<-{chrom}:{start}-{end}"
                file_label = list_line[3]
                w_f.write(fmt_dhs.format(**locals()))

    # split
    pool = Pool(processes=40)
    func_split = partial(split_ref_bed, all_ref)
    pool.map(func_split, list_input)
    pool.close()

    return


def annotate_dhs_promoter(path_promoter_in, dict_in):
    path_out = dict_in['path_out']
    path_in = dict_in['path_in']
    file = dict_in['file']
    file_in = os.path.join(path_in, file)
    bedtools_out = os.path.join(path_out, f"{file}.bedtools.out")
    os.system(f"bedtools intersect -a {file_in} -b {path_promoter_in} -loj "
              f"> {bedtools_out}")
    col_num = int(check_output("head -n 1 " + file_in + " | awk '{print NF}'",
                               shell=True).strip())
    use_col_list = list(range(1, col_num + 1))
    use_col_list.extend([col_num + 4, col_num + 5])
    use_col = ','.join([str(num) for num in use_col_list])
    os.system(f"cut -f {use_col} {bedtools_out} > "
              f"{os.path.join(path_out, file)}")
    os.remove(bedtools_out)

    return


def annotate_dhs_histone(dict_in):
    path_out = dict_in['path_out']
    path_ref = dict_in['path_ref']
    path_in = dict_in['path_in']
    file = dict_in['file']
    file_in = os.path.join(path_in, file)
    bedtools_out = os.path.join(path_out, f"{file}.bedtools.out")
    os.system(f"bedtools intersect -a {path_ref} -b {file_in} -loj "
              f"> {bedtools_out}")
    col_num = int(check_output(
        "head -n 1 " + path_ref + " | awk '{print NF}'",
        shell=True).strip())
    use_col_list = list(range(1, col_num + 1))
    use_col_list.extend([col_num + 4, col_num + 5])
    use_col = ','.join([str(num) for num in use_col_list])
    os.system(f"cut -f {use_col} {bedtools_out} > "
              f"{os.path.join(path_out, file)}")
    # os.remove(bedtools_out)

    return


def annotate_dhs(path_dhs_in, path_promoter_in, path_h3k27ac_in, path_dhs_out):
    path_out_pro = path_dhs_out + '_pro'
    os.system(f"rm -rf {path_out_pro}")
    os.mkdir(path_out_pro)

    list_dhs = generate_file_list(path_dhs_in, path_out_pro)
    pool = Pool(processes=40)
    func_pro = partial(annotate_dhs_promoter, path_promoter_in)
    pool.map(func_pro, list_dhs)
    pool.close()

    list_h3k27ac = generate_file_list(path_h3k27ac_in, path_dhs_out)
    list_h3k27ac_input = []
    for sub_dict in list_h3k27ac:
        if (sub_dict['organ'] != '.') & (sub_dict['life_stage'] != '.'):
            sub_dict['path_ref'] = os.path.join(
                path_out_pro, f"{sub_dict['organ']}/{sub_dict['life_stage']}/"
                              f"{sub_dict['life_stage']}.bed")
            list_h3k27ac_input.append(sub_dict)
    # df_dhs_pro = pd.DataFrame(list_dhs_pro,
    #                           columns=['path_out', 'file', 'path_in',
    #                                    'organ', 'life_stage', 'term'])
    # df_dhs_pro.columns = ['path_out', 'file', 'path_in_1',
    #                       'organ', 'life_stage', 'term']
    # df_h3k27ac = pd.DataFrame(list_h3k27ac,
    #                           columns=['path_out', 'file', 'path_in',
    #                                    'organ', 'life_stage', 'term'])
    # df_h3k27ac.columns = ['path_out', 'file', 'path_in_2',
    #                       'organ', 'life_stage', 'term']
    # df_merge = pd.merge(df_dhs_pro, df_h3k27ac,
    #                     on=['path_out', 'file',
    #                         'organ', 'life_stage', 'term'], how='inner')
    # list_dict = df_merge.to_dict('records')

    pool = Pool(processes=40)
    pool.map(annotate_dhs_histone, list_h3k27ac_input)
    pool.close()

    return


if __name__ == '__main__':
    time_start = time()
    # build DHS reference by organ

    # standardization
    # DHS
    path_dhs = '/home/zy/driver_mutation/data/DHS/hg19'
    path_dhs_stan = '/home/zy/driver_mutation/data/DHS/hg19_standard'
    standardize_bed(path_dhs, path_dhs_stan, 'DHS')

    # H3K27ac
    path_h3k27ac = \
        '/home/zy/driver_mutation/data/ENCODE/histone_ChIP-seq/hg19/H3K27ac'
    path_h3k27ac_stan = \
        '/home/zy/driver_mutation/data/ENCODE/histone_ChIP-seq/' \
        'hg19/H3K27ac_standard'
    standardize_bed(path_h3k27ac, path_h3k27ac_stan, 'H3K27ac')

    # unify DHS labels
    path_dhs_uniform = '/home/zy/driver_mutation/data/DHS/hg19_uniform'
    merge_split_bed(path_dhs_stan, path_dhs_uniform)

    # annotate DHS
    path_promoter = '/home/zy/driver_mutation/data/gene/' \
                    'promoters.up2k.protein.gencode.v19.bed'
    path_anno = '/home/zy/driver_mutation/data/DHS/hg19_annotation'
    annotate_dhs(path_dhs_uniform, path_promoter, path_h3k27ac_stan, path_anno)

    time_end = time()
    print(time_end - time_start)
