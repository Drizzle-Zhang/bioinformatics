#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: annotate_DHS.V0.py
# @time: 10/28/19 4:08 PM

from time import time
import os
from multiprocessing import Pool
from functools import partial
import pandas as pd
import numpy as np
from subprocess import check_output


def generate_file_list(path_in, path_out):
    folder_1 = os.listdir(path_in)
    list_input = []
    for element_1 in folder_1:
        if element_1[:4] == 'meta':
            continue
        path_1 = os.path.join(path_in, element_1)
        if not os.path.isdir(path_1):
            continue
        else:
            path_1_out = os.path.join(path_out, element_1)
            if not os.path.exists(path_1_out):
                os.makedirs(path_1_out)
            folder_2 = os.listdir(path_1)
            for element_2 in folder_2:
                path_2 = os.path.join(path_1, element_2)
                if not os.path.isdir(path_2):
                    continue
                else:
                    path_2_out = os.path.join(path_1_out, element_2)
                    if not os.path.exists(path_2_out):
                        os.makedirs(path_2_out)
                    folder_3 = os.listdir(path_2)
                    for element_3 in folder_3:
                        path_3 = os.path.join(path_2, element_3)
                        if not os.path.isdir(path_3):
                            continue
                        else:
                            path_3_out = os.path.join(path_2_out, element_3)
                            if not os.path.exists(path_3_out):
                                os.makedirs(path_3_out)
                            folder_4 = os.listdir(path_3)
                            for element_4 in folder_4:
                                path_4 = os.path.join(path_3, element_4)
                                if not os.path.isdir(path_4):
                                    if element_4[-4:] != '.bed':
                                        continue
                                    list_input.append(
                                        dict(path_out=path_3_out,
                                             file=element_4, path_in=path_3,
                                             organ=element_1,
                                             life_stage=element_2,
                                             term=element_3))

    return list_input


def sub_stan(type_bed, df_meta, dict_in):
    path_out = dict_in['path_out']
    path_in = dict_in['path_in']
    file = dict_in['file']
    organ = dict_in['organ'].replace('_', ' ')
    life_stage = dict_in['life_stage'].replace('_', ' ')
    term = dict_in['term'].replace('_', ' ').replace(
                             '+', '/').replace("--", "'")
    file_label = f"{organ}|{life_stage}|{term}"
    file_in = os.path.join(path_in, file)
    file_out = os.path.join(path_out, file)
    if type_bed == 'DHS':
        file_tmp = file_in
    else:
        file_tmp = file_out + '.tmp'
        df_sub = df_meta.loc[
            (df_meta['Biosample organ'].apply(
                lambda x: organ in x.strip().split(','))) &
            (df_meta['Biosample life stage'] == life_stage) &
            (df_meta['Biosample term name'] == term), :]
        total_num = np.sum(df_sub['Total peak number']).tolist()
        file_num = df_sub.shape[0]
        os.system(f"Rscript adjust_p_value.R "
                  f"{file_in} {file_tmp} {total_num} {file_num}")

    with open(file_tmp, 'r') as r_f:
        with open(file_out, 'w') as w_f:
            fmt_dhs = "{chrom}\t{start}\t{end}\t{label}\t.\t.\t{file_label}\n"
            fmt_histone = "{chrom}\t{start}\t{end}\t{label}\t" \
                          "{score}\t.\t{file_label}\n"
            for line in r_f:
                list_line = line.strip().split('\t')
                chrom = list_line[0]
                start = list_line[1]
                end = list_line[2]
                label = f"{type_bed}<-{chrom}:{start}-{end}"
                if type_bed == 'DHS':
                    w_f.write(fmt_dhs.format(**locals()))
                else:
                    score = list_line[3]
                    w_f.write(fmt_histone.format(**locals()))

    if type_bed != 'DHS':
        os.system(f"rm {file_tmp}")

    return


def standardize_bed(path_in, path_out, type_bed, num_process):
    if os.path.exists(path_out):
        os.system(f"rm -rf {path_out}")
    os.mkdir(path_out)
    os.system(f"cp {os.path.join(path_in, 'metadata.simple.tsv')} "
              f"{os.path.join(path_out, 'metadata.simple.tsv')}")
    os.system(f"cp {os.path.join(path_in, 'meta.reference.tsv')} "
              f"{os.path.join(path_out, 'meta.reference.tsv')}")

    df_meta = pd.read_csv(
        os.path.join(path_in, 'metadata.simple.tsv'), sep='\t'
    )

    list_input = generate_file_list(path_in, path_out)
    pool = Pool(processes=num_process)
    func_stan = partial(sub_stan, type_bed, df_meta)
    pool.map(func_stan, list_input)
    pool.close()

    return


def split_ref_bed(ref_file, dict_in):
    organ = dict_in['organ']
    life_stage = dict_in['life_stage']
    term = dict_in['term']
    label = '|'.join([organ, life_stage, term])
    with open(ref_file, 'r') as r_ref:
        with open(os.path.join(dict_in['path_out'], dict_in['file']), 'w') \
                as w_f:
            fmt_dhs = "{chrom}\t{start}\t{end}\t{dhs_id}\t.\t.\t{label}\n"
            for line in r_ref:
                list_line = line.strip().split('\t')
                chrom = list_line[0]
                start = list_line[1]
                end = list_line[2]
                dhs_id = list_line[3]
                list_label = list_line[6].strip().split(',')
                if label in list_label:
                    w_f.write(fmt_dhs.format(**locals()))

    return


def merge_bed(path_bed, dict_in):
    flank_percent = dict_in['flank_percent']
    term_name = dict_in['term_name']
    path_out = dict_in['path']
    col_collapse = '2,3,5,7'
    str_collapse = \
        ','.join([val for val in ['collapse']
                  for i in range(len(col_collapse.split(',')))])
    cat_out = os.path.join(path_out, f"{term_name}.bed.cat")
    sort_out = os.path.join(path_out, f"{term_name}.bed.sort")
    merge_out = os.path.join(path_out, f"{term_name}.bed.merge")
    cat_in = ' '.join([os.path.join(path_bed, acce_id + '.bed')
                       for acce_id in dict_in['accession_ids']])
    os.system(f"cat {cat_in} > {cat_out}")
    os.system(f"bedtools sort -i {cat_out} > {sort_out}")
    # os.system(f"sort -k 1,1 -k2,2n {cat_out} > {sort_out}")
    os.system(f"bedtools merge -i {sort_out} "
              f"-c {col_collapse} -o {str_collapse} > {merge_out}")

    # split merge file
    split_out = os.path.join(path_out, f"{term_name}.bed")
    with open(merge_out, 'r') as r_f:
        with open(split_out, 'w') as w_f:
            fmt = "{chrom}\t{start}\t{end}\t{label}\t{p_value}\n"
            for line in r_f:
                list_line = line.strip().split('\t')
                chrom = list_line[0]
                array_start = np.array(
                    [int(num) for num in list_line[3].strip().split(',')])
                array_end = np.array(
                    [int(num) for num in list_line[4].strip().split(',')])
                array_label = np.array(
                    [num for num in list_line[6].strip().split(',')])
                array_p_value = np.array(
                    [num for num in list_line[5].strip().split(',')])
                array_length = array_end - array_start
                while array_length.shape[0] > 0:
                    idx = np.argmax(array_length)
                    start_idx = array_start[idx]
                    end_idx = array_end[idx]
                    flank = array_length[idx] * flank_percent
                    select_bool = (array_start >= start_idx - flank) & \
                                  (array_end <= end_idx + flank)
                    select_start = array_start[select_bool]
                    select_end = array_end[select_bool]
                    select_label = array_label[select_bool]
                    select_p_value = array_p_value[select_bool]
                    # reset arrays
                    array_start = array_start[~select_bool]
                    array_end = array_end[~select_bool]
                    array_label = array_label[~select_bool]
                    array_p_value = array_p_value[~select_bool]
                    array_length = array_end - array_start
                    # write new rows
                    start = str(np.min(select_start))
                    end = str(np.max(select_end))
                    label = \
                        ','.join(map(str, select_label.tolist()))
                    p_value = ','.join(map(str, select_p_value.tolist()))
                    w_f.write(fmt.format(**locals()))

    os.remove(cat_out)
    os.remove(sort_out)

    return


def sub_merge(dict_in):
    sub_meta = dict_in['sub_meta']
    sub_path_in = dict_in['path_in']
    sub_path_out = dict_in['path_out']
    file_out = \
        os.path.join(sub_path_out, dict_in['organ'].replace(' ', '_') + '.bed')
    if not os.path.exists(sub_path_out):
        os.mkdir(sub_path_out)
    list_meta = sub_meta.to_dict("records")
    accession_ids = [
        '/'.join([sub_dict['Biosample life stage'].replace(' ', '_'),
                  sub_dict['Biosample term name'].replace(
                        ' ', '_').replace('/', '+').replace("'", '--'),
                  sub_dict['Biosample term name'].replace(
                        ' ', '_').replace('/', '+').replace("'", '--')])
        for sub_dict in list_meta]
    if sub_meta.shape[0] == 1:
        os.system(
            f"cp {os.path.join(sub_path_in, accession_ids[0] + '.bed')} "
            f"{file_out}"
        )
    elif sub_meta.shape[0] > 1:
        dict_merge = dict(
            path=sub_path_out,
            term_name=dict_in['organ'].replace(' ', '_'),
            accession_ids=accession_ids,
            flank_percent=0.3)
        merge_bed(sub_path_in, dict_merge)
        labels = [f"{dict_in['organ']}|{sub_dict['Biosample life stage']}|"
                  f"{sub_dict['Biosample term name']}"
                  for sub_dict in list_meta]
        list_bed = \
            (pd.read_csv(file_out, sep='\t', header=None)).to_dict('records')
        list_dict = []
        for sub_dict in list_bed:
            out_dict = dict()
            out_dict['peak_id'] = \
                f"{sub_dict[0]}:{str(sub_dict[1])}-{str(sub_dict[2])}"
            list_label = sub_dict[3].strip().split(',')
            list_lgp = \
                [num for num in sub_dict[4].strip().split(',')]
            for label in labels:
                try:
                    idx = list_label.index(label)
                except ValueError:
                    out_dict[label] = 0
                    continue
                lgp = list_lgp[idx]
                if lgp != '.':
                    out_dict[label] = float(lgp)
                else:
                    out_dict[label] = -1
            list_dict.append(out_dict)
        df_label_peak = pd.DataFrame(list_dict, columns=['peak_id'] + labels)
        df_label_peak.index = df_label_peak['peak_id']
        df_label_peak = df_label_peak.drop('peak_id', 1)
        mat_peak = os.path.join(sub_path_out, 'label_peak.txt')
        df_label_peak.to_csv(mat_peak, sep='\t')
        os.system(f"Rscript scatter.plot.R {mat_peak} "
                  f"{os.path.join(sub_path_out, 'scatter.pdf')}")

    return


def merge_organ_cluster(path_in, path_out, num_process):
    if os.path.exists(path_out):
        os.system(f"rm -rf {path_out}")
    os.mkdir(path_out)
    os.system(f"cp {os.path.join(path_in, 'meta.reference.tsv')} "
              f"{os.path.join(path_out, 'meta.reference.tsv')}")

    df_meta = \
        pd.read_csv(os.path.join(path_in, 'meta.reference.tsv'), sep='\t')
    organs = []
    for line in df_meta['Biosample organ'].tolist():
        organs.extend(line.strip().split(','))
    organs = set(organs)
    list_input = []
    for organ in organs:
        sub_meta = df_meta.loc[
                   df_meta['Biosample organ'].apply(
                       lambda x: organ in x.strip().split(',')), :]
        list_input.append(
            dict(sub_meta=sub_meta, organ=organ,
                 path_out=os.path.join(path_out, organ.replace(' ', '_')),
                 path_in=os.path.join(path_in, organ.replace(' ', '_'))))

    pool = Pool(processes=num_process)
    pool.map(sub_merge, list_input)
    pool.close()

    return


def annotate_dhs_promoter(path_promoter_in, dict_in):
    path_out = dict_in['path_out']
    path_in = dict_in['path_in']
    file = dict_in['file']
    file_in = os.path.join(path_in, file)
    bedtools_out = os.path.join(path_out, f"{file}.bedtools.out")
    col_num = int(check_output("head -n 1 " + file_in + " | awk '{print NF}'",
                               shell=True).strip())
    use_col_list = list(range(1, col_num + 1))
    use_col_list.extend([col_num + 4, col_num + 5])
    use_col = ','.join([str(num) for num in use_col_list])
    os.system(f"bedtools intersect -a {file_in} -b {path_promoter_in} -loj "
              f"> {bedtools_out}")
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
    col_num = int(check_output(
        "head -n 1 " + path_ref + " | awk '{print NF}'",
        shell=True).strip())
    use_col_list = list(range(1, col_num + 1))
    use_col_list.extend([col_num + 4, col_num + 5, col_num + 7])
    use_col = ','.join([str(num) for num in use_col_list])
    bedtools_out = os.path.join(path_out, f"{file}.bedtools.out")
    os.system(f"bedtools intersect -a {path_ref} -b {file_in} -loj "
              f"> {bedtools_out}")
    os.system(f"cut -f {use_col} {bedtools_out} > "
              f"{os.path.join(path_out, file)}")
    os.remove(bedtools_out)

    return


def match_dhs_file(list_ref, list_histone):
    meta_dhs = pd.DataFrame(list_ref)
    df_all = meta_dhs.loc[meta_dhs['organ'] == '.', :]
    list_out = []
    for dict_line in list_histone:
        organ = dict_line['organ']
        life_stage = dict_line['life_stage']
        term = dict_line['term']
        df_organ = meta_dhs.loc[meta_dhs['organ'] == organ, :]
        if df_organ.shape[0] == 0:
            dict_line['path_ref'] = os.path.join(
                df_all['path_in'].tolist()[0], df_all['file'].tolist()[0]
            )
            list_out.append(dict_line)
        else:
            df_life = df_organ.loc[df_organ['life_stage'] == life_stage, :]
            if df_life.shape[0] == 0:
                df_all_life = df_organ.loc[df_organ['life_stage'] == '.', :]
                dict_line['path_ref'] = os.path.join(
                    df_all_life['path_in'].tolist()[0],
                    df_all_life['file'].tolist()[0]
                )
                list_out.append(dict_line)
            else:
                df_term = df_life.loc[df_life['term'] == term, :]
                if df_term.shape[0] == 0:
                    df_all_term = df_life.loc[df_life['term'] == '.', :]
                    dict_line['path_ref'] = os.path.join(
                        df_all_term['path_in'].tolist()[0],
                        df_all_term['file'].tolist()[0]
                    )
                    list_out.append(dict_line)
                else:
                    dict_line['path_ref'] = os.path.join(
                        df_term['path_in'].tolist()[0],
                        df_term['file'].tolist()[0]
                    )
                    list_out.append(dict_line)

    return list_out


def annotate_dhs(path_dhs_in, path_promoter_in, path_h3k27ac_in,
                 path_h3k4me3_in, path_dhs_out, num_process):
    path_out_pro = path_dhs_out + '_pro'
    if os.path.exists(path_out_pro):
        os.system(f"rm -rf {path_out_pro}")
    os.mkdir(path_out_pro)

    list_dhs = generate_file_list(path_dhs_in, path_out_pro)
    pool = Pool(processes=40)
    func_pro = partial(annotate_dhs_promoter, path_promoter_in)
    pool.map(func_pro, list_dhs)
    pool.close()

    path_out_tmp = path_dhs_out + '_tmp'
    if os.path.exists(path_out_tmp):
        os.system(f"rm -rf {path_out_tmp}")
    os.mkdir(path_out_tmp)
    list_dhs_pro = generate_file_list(path_out_pro, path_out_tmp)
    list_h3k4me3 = generate_file_list(path_h3k4me3_in, path_out_tmp)
    list_h3k4me3_input = match_dhs_file(list_dhs_pro, list_h3k4me3)
    pool = Pool(processes=num_process)
    pool.map(annotate_dhs_histone, list_h3k4me3_input)
    pool.close()

    if os.path.exists(path_dhs_out):
        os.system(f"rm -rf {path_dhs_out}")
    os.mkdir(path_dhs_out)
    list_dhs_tmp = generate_file_list(path_out_tmp, path_dhs_out)
    list_h3k27ac = generate_file_list(path_h3k27ac_in, path_dhs_out)
    list_h3k27ac_input = match_dhs_file(list_dhs_tmp, list_h3k27ac)
    pool = Pool(processes=num_process)
    pool.map(annotate_dhs_histone, list_h3k27ac_input)
    pool.close()

    return


def sub_anno_cre(dict_in):
    path_out = dict_in['path_out']
    path_in = dict_in['path_in']
    file = dict_in['file']
    file_in = os.path.join(path_in, file)
    file_out = os.path.join(path_out, file)
    file_out_tmp = file_out + '.tmp'
    file_out_origin = file_out + '.origin'
    with open(file_in, 'r') as r_f:
        with open(file_out_tmp, 'w') as w_f:
            fmt_dhs = "{chrom}\t{start}\t{end}\t{dhs_id}\t{score}\t.\t" \
                      "{cre}\n"
            with open(file_out_origin, 'w') as w_ori:
                for line in r_f:
                    list_line = line.strip().split('\t')
                    chrom = list_line[0]
                    start = list_line[1]
                    end = list_line[2]
                    dhs_id = list_line[3]
                    promoter = list_line[7]
                    h3k4me3 = list_line[9]
                    h3k27ac = list_line[12]
                    if (promoter != '.') & (h3k4me3 != '.'):
                        cre = 'Promoter'
                        score = list_line[10]
                    elif h3k27ac != '.':
                        cre = 'Enhancer'
                        score = list_line[13]
                    else:
                        cre = 'Other'
                        score = '.'
                    list_line.insert(6, cre)
                    list_line[4] = score
                    w_ori.write('\t'.join(list_line) + '\n')
                    w_f.write(fmt_dhs.format(**locals()))

    os.system(f"uniq {file_out_tmp} > {file_out}")

    return


def annotate_cre(path_in, path_out, num_process):
    if os.path.exists(path_out):
        os.system(f"rm -rf {path_out}")
    os.mkdir(path_out)

    list_dhs_anno = generate_file_list(path_in, path_out)
    pool = Pool(processes=num_process)
    pool.map(sub_anno_cre, list_dhs_anno)
    pool.close()

    return


if __name__ == '__main__':
    time_start = time()
    num_cpu = 40
    # build DHS reference by organ

    # standardization
    # DHS
    path_dhs = '/home/zy/driver_mutation/data/DHS/GRCh38tohg19'
    path_dhs_stan = '/home/zy/driver_mutation/data/' \
                    'DHS/GRCh38tohg19_standard'
    standardize_bed(path_dhs, path_dhs_stan, 'DHS', num_cpu)
    print('Standardization of DHS completed!')

    # H3K27ac
    path_h3k27ac = \
        '/home/zy/driver_mutation/data/ENCODE/' \
        'histone_ChIP-seq/GRCh38tohg19/H3K27ac_merge'
    path_h3k27ac_stan = \
        '/home/zy/driver_mutation/data/ENCODE/' \
        'histone_ChIP-seq/GRCh38tohg19/H3K27ac_standard'
    standardize_bed(path_h3k27ac, path_h3k27ac_stan, 'H3K27ac', num_cpu)
    print('Standardization of H3K27ac completed!')

    # H3K4me3
    path_h3k4me3 = \
        '/home/zy/driver_mutation/data/ENCODE/' \
        'histone_ChIP-seq/GRCh38tohg19/H3K4me3_merge'
    path_h3k4me3_stan = \
        '/home/zy/driver_mutation/data/ENCODE/' \
        'histone_ChIP-seq/GRCh38tohg19/H3K4me3_standard'
    standardize_bed(path_h3k4me3, path_h3k4me3_stan, 'H3K4me3', num_cpu)
    print('Standardization of H3K4me3 completed!')

    # merge and cluster
    # DHS
    path_dhs_cluster = \
        '/home/zy/driver_mutation/data/DHS/GRCh38tohg19_cluster'
    merge_organ_cluster(path_dhs_stan, path_dhs_cluster, num_cpu)

    # H3K27ac
    path_h3k27ac_cluster = \
        '/home/zy/driver_mutation/data/ENCODE/' \
        'histone_ChIP-seq/GRCh38tohg19/H3K27ac_cluster'
    merge_organ_cluster(path_h3k27ac_stan, path_h3k27ac_cluster, num_cpu)

    # H3K4me3
    path_h3k4me3_cluster = \
        '/home/zy/driver_mutation/data/ENCODE/' \
        'histone_ChIP-seq/GRCh38tohg19/H3K4me3_cluster'
    merge_organ_cluster(path_h3k4me3_stan, path_h3k4me3_cluster, num_cpu)

    # # unify DHS labels
    # path_dhs_uniform = '/home/zy/driver_mutation/data/' \
    #                    'DHS/GRCh38tohg19_uniform'
    # # merge_split_bed(path_dhs_stan, path_dhs_uniform, num_cpu)
    # print('Uniform of DHS completed!')
    #
    # # annotate DHS
    # path_promoter = '/home/zy/driver_mutation/data/gene/' \
    #                 'promoters.up2k.protein.gencode.v19.bed'
    # path_anno = '/home/zy/driver_mutation/data/DHS/' \
    #             'GRCh38tohg19_annotation'
    # annotate_dhs(path_dhs_uniform, path_promoter, path_h3k27ac_stan,
    #              path_h3k4me3_stan, path_anno, num_cpu)
    # print('Annotation of DHS completed!')
    #
    # # annotate cRE
    # path_cre = '/home/zy/driver_mutation/data/DHS/' \
    #            'GRCh38tohg19_cRE'
    # annotate_cre(path_anno, path_cre, num_cpu)
    # print('Annotation completed!')

    time_end = time()
    print(time_end - time_start)
