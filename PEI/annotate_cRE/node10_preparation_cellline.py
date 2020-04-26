#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: node10_preparation.py
# @time: 12/3/19 10:52 AM

from time import time
import pandas as pd
import numpy as np
import os
from collections import defaultdict
from multiprocessing import Pool
from functools import partial
from subprocess import check_output
from itertools import combinations
from scipy.spatial.distance import pdist
from sklearn.preprocessing import StandardScaler


def filter_meta(meta_in):
    # reserve released data and original(untreated) peak file
    # simplify meta file
    df_meta = pd.read_csv(meta_in, sep='\t')
    df_meta_released = df_meta.loc[
        (df_meta['File Status'] == 'released') &
        (df_meta['Output type'] == 'peaks') &
        (df_meta['Biological replicate(s)'].apply(
            lambda x: len(str(x).split(', ')) == 1
        )) &
        (df_meta['Biosample treatments'].apply(
            lambda x: np.isnan(x) if isinstance(x, float) else False)) &
        (df_meta['Biosample genetic modifications methods'].apply(
            lambda x: np.isnan(x) if isinstance(x, float) else False)),
        ['File accession', 'Experiment accession', 'Biosample term id',
         'Biosample term name', 'Biosample type', 'Biosample treatments',
         'Biosample genetic modifications methods', 'Output type', 'Assembly',
         'Biological replicate(s)', 'File Status']]
    df_meta_released.index = df_meta_released['File accession']

    return df_meta_released


def build_dict_attr(path_meta):
    # build a dictionary to add organ/life stage/cell labels
    dict_accession_attr = defaultdict(list)
    meta_files = os.listdir(path_meta)
    for file in meta_files:
        if file[:4] == 'meta':
            attr = file[9:-4].replace('_', ' ')
            with open(os.path.join(path_meta, file), 'r') as r_sub:
                for line in r_sub:
                    list_line = line.strip().split('\t')
                    if list_line[0][:3] == 'ENC':
                        dict_accession_attr[list_line[0]].append(attr)

    return dict_accession_attr


def add_attr(df_meta, dict_attr, column_name):
    # add labels to meta file
    df_access = df_meta['File accession']
    df_attr = df_access.apply(lambda x: ','.join(dict_attr[x]))
    df_attr.name = column_name
    df_out = pd.concat([df_meta, df_attr], axis=1, sort=False)

    return df_out


def modify_meta(df_meta):
    # only select cell line data
    df_meta = df_meta.loc[df_meta['Biosample type'] == 'cell line', :]

    return df_meta


def sub_hg38tohg19(path_hg38, path_hg19, dict_in):
    file_hg38 = os.path.join(path_hg38, dict_in['File accession'] + '.bed')
    file_hg19 = os.path.join(path_hg19, dict_in['File accession'] + '.bed')
    # Download from UCSC
    file_chain = \
        '/local/zy/tools/files_liftOver/hg38ToHg19.over.chain.gz'
    file_ummap = os.path.join(
        path_hg19, dict_in['File accession'] + '.bed.unmap')
    if dict_in['Assembly'] == 'hg19':
        with open(file_hg19, 'w') as w_f:
            fmt = "{chrom}\t{start}\t{end}\t{peak_id}\t{score}\t{strand}\t" \
                  "{fold_change}\t{p_value}\t{q_value}\t{peak_location}\n"
            with open(file_hg38, 'r') as r_hg38:
                for line in r_hg38:
                    list_line = line.strip().split('\t')
                    dict_hg19 = dict(
                        chrom=list_line[0], start=list_line[1],
                        end=list_line[2], peak_id=dict_in['File accession'],
                        score=list_line[4], strand=list_line[5],
                        fold_change=list_line[6], p_value=list_line[7],
                        q_value=list_line[8], peak_location=list_line[9]
                    )
                    w_f.write(fmt.format(**dict_hg19))
        df_bed = pd.read_csv(file_hg19, sep='\t', header=None)
        scores = np.array(df_bed.iloc[:, 6]).reshape(-1, 1)
        scale_scores = StandardScaler().fit_transform(scores)
        df_bed.iloc[:, 6] = scale_scores
        df_bed.to_csv(file_hg19, sep='\t', index=None, header=None)

    if dict_in['Assembly'] == 'GRCh38':
        file_hg38_labeled = file_hg38 + '.labeled'
        label_assess = str(
            check_output(f"head -1 {file_hg38}", shell=True).strip()
        ).split('\\t')[3]
        if label_assess == '.':
            with open(file_hg38_labeled, 'w') as w_f:
                with open(file_hg38, 'r') as r_f:
                    for i, line in enumerate(r_f):
                        list_line = line.strip().split('\t')
                        list_line[3] = 'peak_' + str(i)
                        w_f.write('\t'.join(list_line) + '\n')
        else:
            os.system(f"cp {file_hg38} {file_hg38_labeled}")

        file_prefix = file_hg38 + '.prefix'
        file_suffix = file_hg38 + '.suffix'
        file_hg19_prefix = file_hg19 + '.prefix'
        file_hg19_format = file_hg19 + '.format'
        os.system(f"cut -f 1,2,3,4 {file_hg38_labeled} > {file_prefix}")
        os.system(f"cut -f 4,5,6,7,8,9,10 {file_hg38_labeled} > {file_suffix}")
        os.system(f"/local/zy/tools/liftOver {file_prefix} {file_chain} "
                  f"{file_hg19_prefix} {file_ummap}")
        dict_peak_score = defaultdict(list)
        with open(file_suffix, 'r') as r_f:
            for line in r_f:
                list_line = line.strip().split('\t')
                dict_peak_score[list_line[0]].append(list_line[1:])
        with open(file_hg19_format, 'w') as w_f:
            fmt = "{chrom}\t{start}\t{end}\t{peak_id}\t{score}\t{strand}\t" \
                  "{fold_change}\t{p_value}\t{q_value}\t{peak_location}\n"
            with open(file_hg19_prefix, 'r') as r_hg19:
                for line in r_hg19:
                    list_line = line.strip().split('\t')
                    list_suffix = dict_peak_score[list_line[3]][0]
                    dict_hg19 = dict(
                        chrom=list_line[0], start=list_line[1],
                        end=list_line[2], peak_id=dict_in['File accession'],
                        score=list_suffix[0], strand=list_suffix[1],
                        fold_change=list_suffix[2], p_value=list_suffix[3],
                        q_value=list_suffix[4], peak_location=list_suffix[5]
                    )
                    w_f.write(fmt.format(**dict_hg19))

        df_old = pd.read_csv(file_prefix, sep='\t', header=None)
        length_old = df_old.iloc[:, 2] - df_old.iloc[:, 1]
        up_limit = min(np.max(length_old) + 10, 15000)
        down_limit = max(np.min(length_old) - 10, 140)
        df_bed = pd.read_csv(file_hg19_format, sep='\t', header=None)
        length = df_bed.iloc[:, 2] - df_bed.iloc[:, 1]
        df_bed = df_bed.loc[(length < up_limit) & (length > down_limit), :]
        scores = np.array(df_bed.iloc[:, 6]).reshape(-1, 1)
        scale_scores = StandardScaler().fit_transform(scores)
        df_bed.iloc[:, 6] = scale_scores
        df_bed.to_csv(file_hg19, sep='\t', index=None, header=None)

        os.remove(file_hg38_labeled)
        os.remove(file_ummap)
        os.remove(file_prefix)
        os.remove(file_suffix)
        os.remove(file_hg19_prefix)
        os.remove(file_hg19_format)

    return


def calculate_peak_numbers(path_in, dict_in):
    file_in = os.path.join(path_in, dict_in['File accession'] + '.bed')
    file_tail = check_output(f"sort -k 8 -n {file_in} | head -1", shell=True)
    list_tail = str(file_tail).strip().split('\\t')
    order = int(str(check_output(f"wc -l {file_in}",
                                 shell=True).strip()).split(' ')[0][2:])
    p_value = 10 ** (-float(list_tail[7]))
    q_value = 10 ** (-float(list_tail[8]))
    if list_tail[6] == -1:
        total = '.'
    else:
        total = int(q_value/p_value * order)

    return {'File accession': dict_in['File accession'],
            'Peak number': order, 'Inferred peak number': total}


def hg38tohg19(path_hg38, path_hg19, path_meta, num_process):
    if os.path.exists(path_hg19):
        os.system(f"rm -rf {path_hg19}")
    os.mkdir(path_hg19)

    df_meta = pd.read_csv(path_meta, sep='\t')
    list_meta = []
    experiments = set(df_meta['Experiment accession'].tolist())
    for exp in experiments:
        df_exp = df_meta.loc[df_meta['Experiment accession'] == exp, :]
        hg19 = df_exp.loc[df_exp['Assembly'] == 'hg19', :]
        hg38 = df_exp.loc[df_exp['Assembly'] == 'GRCh38', :]
        if hg19.shape[0] == hg38.shape[0]:
            list_meta.append(hg19)
        elif hg19.shape[0] > hg38.shape[0]:
            list_meta.append(hg19)
        elif hg19.shape[0] < hg38.shape[0]:
            list_meta.append(hg38)
    new_meta = pd.concat(list_meta, sort=False)
    list_dict = new_meta.to_dict('records')

    pool = Pool(processes=num_process)
    func_hg38tohg19 = partial(sub_hg38tohg19, path_hg38, path_hg19)
    pool.map(func_hg38tohg19, list_dict)
    pool.close()

    # count peak numbers and infer total peak numbers
    pool = Pool(processes=num_process)
    func_calc = partial(calculate_peak_numbers, path_hg19)
    result = pool.map(func_calc, list_dict)
    pool.close()
    df_res = pd.DataFrame(result)
    df_meta_merge = pd.merge(new_meta, df_res, on='File accession')

    df_meta_merge.to_csv(os.path.join(path_hg19, 'metadata.simple.tsv'),
                         sep='\t', index=None)

    return


def merge_bed(path_bed, dict_in):
    flank_percent = dict_in['flank_percent']
    term_name = dict_in['term_name']
    path_out = dict_in['path']
    accession_ids = dict_in['accession_ids']
    if len(accession_ids) == 1:
        file_in = os.path.join(path_bed, accession_ids[0] + '.bed')
        file_out = os.path.join(path_out, f"{term_name}.bed")
        os.system(f"cp {file_in} {file_out}")

        return

    # merge bed files from same term
    col_collapse = '2,3,4,5,6,7,8,9'
    str_collapse = \
        ','.join([val for val in ['collapse']
                  for i in range(len(col_collapse.split(',')))])
    cat_out = os.path.join(path_out, f"{term_name}.bed.cat")
    sort_out = os.path.join(path_out, f"{term_name}.bed.sort")
    merge_out = os.path.join(path_out, f"{term_name}.bed.merge")
    cat_in = ' '.join([os.path.join(path_bed, acce_id + '.bed')
                       for acce_id in accession_ids])
    os.system(f"cat {cat_in} > {cat_out}")
    os.system(f"bedtools sort -i {cat_out} > {sort_out}")
    os.system(f"bedtools merge -i {sort_out} "
              f"-c {col_collapse} -o {str_collapse} > {merge_out}")

    # split merge file
    split_out = os.path.join(path_out, f"{term_name}.bed.unsort")
    split_sort_out = os.path.join(path_out, f"{term_name}.bed.split.sort")
    final_out = os.path.join(path_out, f"{term_name}.bed")
    with open(merge_out, 'r') as r_f:
        with open(split_out, 'w') as w_f:
            fmt = "{chrom}\t{start}\t{end}\t{access}\t0\t.\t" \
                  "{fold_change}\t{p_value}\t-1\t0\n"
            for line in r_f:
                list_line = line.strip().split('\t')
                chrom = list_line[0]
                array_start = np.array(
                    [int(num) for num in list_line[3].strip().split(',')])
                array_end = np.array(
                    [int(num) for num in list_line[4].strip().split(',')])
                array_access = np.array(
                    [num for num in list_line[5].strip().split(',')])
                array_fold_change = np.array(
                    [num for num in list_line[8].strip().split(',')])
                array_p_value = np.array(
                    [num for num in list_line[9].strip().split(',')])
                array_length = array_end - array_start
                while array_length.shape[0] > 0:
                    idx = np.argmax(array_length)
                    start_idx = array_start[idx]
                    end_idx = array_end[idx]
                    # merge similar peaks
                    flank = array_length[idx] * flank_percent
                    select_bool = (array_start >= start_idx - flank) & \
                                  (array_end <= end_idx + flank)
                    select_start = array_start[select_bool]
                    select_end = array_end[select_bool]
                    select_access = array_access[select_bool]
                    select_fold_change = array_fold_change[select_bool]
                    select_p_value = array_p_value[select_bool]
                    # reset arrays
                    array_start = array_start[~select_bool]
                    array_end = array_end[~select_bool]
                    array_access = array_access[~select_bool]
                    array_fold_change = array_fold_change[~select_bool]
                    array_p_value = array_p_value[~select_bool]
                    array_length = array_end - array_start
                    # write new rows
                    start = str(np.min(select_start))
                    end = str(np.max(select_end))
                    access = '|'.join(map(str, select_access.tolist()))
                    fold_change = \
                        '|'.join(map(str, select_fold_change.tolist()))
                    p_value = '|'.join(map(str, select_p_value.tolist()))
                    w_f.write(fmt.format(**locals()))

    os.system(f"bedtools sort -i {split_out} > {split_sort_out}")
    with open(final_out, 'w') as w_final:
        old_end = 0
        old_chrom = '0'
        with open(split_sort_out, 'r') as r_sort:
            for line in r_sort:
                list_line = line.strip().split('\t')
                chrom = list_line[0]
                if chrom != old_chrom:
                    old_end = 0
                start = int(list_line[1])
                if start > old_end:
                    w_final.write(line)
                else:
                    list_line[1] = str(old_end + 1)
                    w_final.write('\t'.join(list_line) + '\n')
                old_end = int(list_line[2])
                old_chrom = chrom

    os.remove(cat_out)
    os.remove(sort_out)
    os.remove(merge_out)
    os.remove(split_out)
    os.remove(split_sort_out)

    return


def overlap_matrix(path_in, dict_in):
    term_name = dict_in['term_name']
    path_out = dict_in['path']
    accession_ids = dict_in['accession_ids']
    df_meta = pd.read_csv(
        os.path.join(path_in, 'metadata.simple.tsv'), sep='\t')

    if len(accession_ids) <= 1:
        return
    else:
        file_merge = os.path.join(path_out, term_name + '.bed')
        list_bed = \
            (pd.read_csv(file_merge, sep='\t', header=None)).to_dict('records')
        list_dict = []
        for sub_dict in list_bed:
            out_dict = dict()
            out_dict['peak_id'] = \
                f"{sub_dict[0]}:{str(sub_dict[1])}-{str(sub_dict[2])}"
            set_access = set(sub_dict[3].strip().split('|'))
            if accession_ids[0][:5] == 'ENCSR':
                for access in accession_ids:
                    access_files = (df_meta.loc[
                        df_meta['Experiment accession'] == access,
                        'File accession']).tolist()
                    if set(access_files).intersection(set_access):
                        out_dict[access] = 1
                    else:
                        out_dict[access] = 0
                    list_dict.append(out_dict)
            else:
                for access in accession_ids:
                    if access in set_access:
                        out_dict[access] = 1
                    else:
                        out_dict[access] = 0
                list_dict.append(out_dict)
        df_ref = pd.DataFrame(
            list_dict, columns=['peak_id'] + accession_ids
        )
        df_ref.index = df_ref['peak_id']
        df_ref = df_ref.drop('peak_id', 1)
        # write score matrix to txt file
        file_ref = os.path.join(path_out, term_name + '.ref')
        df_ref.to_csv(file_ref, sep='\t')

        list_out = []
        list_com = combinations(df_ref.columns, 2)
        for com in list_com:
            jdist = pdist(
                (df_ref.loc[
                    (df_ref.loc[:, com[0]] != 0) |
                    (df_ref.loc[:, com[1]] != 0), com]).T,
                'jaccard')[0]
            list_out.append({'Name': term_name, 'Combination': com,
                             'Jaccard distance': jdist})

        df_out = pd.DataFrame(list_out)

        return df_out


def merge_experiment(path_in, path_out, flank_percent, num_process):
    # integrate accession files from same experiment to a single file
    df_meta = pd.read_csv(
        os.path.join(path_in, 'metadata.simple.tsv'), sep='\t')

    if os.path.exists(path_out):
        os.system(f"rm -rf {path_out}")
    os.mkdir(path_out)
    os.system(f"cp {os.path.join(path_in, 'metadata.simple.tsv')} "
              f"{os.path.join(path_out, 'metadata.simple.tsv')}")

    list_input = []
    experiments = set(df_meta['Experiment accession'].tolist())
    for experiment in experiments:
        df_exp = df_meta.loc[df_meta['Experiment accession'] == experiment, :]
        accession_ids = list(set(df_exp['File accession'].tolist()))
        list_input.append(
            dict(path=path_out,
                 term_name=experiment.replace(' ', '_').replace(
                     '/', '+').replace("'", '--'),
                 accession_ids=accession_ids,
                 flank_percent=flank_percent))

    pool = Pool(processes=num_process)
    func_merge = partial(merge_bed, path_in)
    pool.map(func_merge, list_input)
    pool.close()

    pool = Pool(processes=num_process)
    func_overlap = partial(overlap_matrix, path_in)
    list_df = pool.map(func_overlap, list_input)
    pool.close()

    df_overlap = pd.concat(list_df, sort=False)
    df_overlap.to_csv(
        os.path.join(path_out, 'overlap.txt'), sep='\t', index=None
    )

    return


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
                 accession_ids=accession_ids,
                 flank_percent=flank_percent))

    pool = Pool(processes=num_process)
    func_merge = partial(merge_bed, path_in)
    pool.map(func_merge, list_input)
    pool.close()

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
        file_out = os.path.join(path_out, file)
        # make folder
        if not os.path.exists(folder3):
            os.mkdir(folder3)
    else:
        file_in = os.path.join(path_in, dict_in['File accession'] + '.bed')
        file_out = os.path.join(path_out, dict_in['File accession'] + '.bed')

    with open(file_in, 'r') as r_f:
        with open(file_out, 'w') as w_f:
            fmt_dhs = "{chrom}\t{start}\t{end}\t{label}\t{score}\t.\t" \
                      "{file_label}\t{accessions}\n"
            fmt_chip = "{chrom}\t{start}\t{end}\t{label}\t" \
                       "{score}\t{pvalue}\t{accessions}\n"
            for line in r_f:
                list_line = line.strip().split('\t')
                chrom = list_line[0]
                start = list_line[1]
                end = list_line[2]
                accessions = list_line[3]
                label = f"{type_bed}<-{chrom}:{start}-{end}"
                score = max([float(num)
                             for num in list_line[6].strip().split('|')])
                if type_bed == 'DHS':
                    w_f.write(fmt_dhs.format(**locals()))
                elif type_bed in {'H3K4me3', 'H3K27ac'}:
                    pvalue = list_line[7]
                    w_f.write(fmt_chip.format(**locals()))
                elif type_bed == 'CTCF':
                    pvalue = list_line[8]
                    w_f.write(fmt_chip.format(**locals()))

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


if __name__ == '__main__':
    time_start = time()
    # parameters
    num_cpu = 40

    # build life stage dictionary
    path_lifestage = '/local/zy/PEI/origin_data/ENCODE/metadata/life_stage'
    dict_lifestage = build_dict_attr(path_lifestage)

    # build organ dictionary
    path_meta_organ = '/local/zy/PEI/origin_data/ENCODE/metadata/organ'
    dict_organ = build_dict_attr(path_meta_organ)

    # build organ dictionary
    path_meta_cell = '/local/zy/PEI/origin_data/ENCODE/metadata/cell'
    dict_cell = build_dict_attr(path_meta_cell)

    # build organ dictionary
    path_meta_lab = '/local/zy/PEI/origin_data/ENCODE/metadata/lab'
    dict_lab = build_dict_attr(path_meta_lab)
    print("Preparation of dictionary files and reference files is completed")

    # DHS reference
    # metafile
    path_dhs = '/local/zy/PEI/origin_data/ENCODE/DNase-seq/all'
    ori_meta_dhs = os.path.join(path_dhs, 'metadata.tsv')
    df_meta_dhs = filter_meta(ori_meta_dhs)
    df_meta_dhs = add_attr(df_meta_dhs, dict_lifestage, 'Biosample life stage')
    df_meta_dhs = add_attr(df_meta_dhs, dict_organ, 'Biosample organ')
    df_meta_dhs = add_attr(df_meta_dhs, dict_cell, 'Biosample cell')
    df_meta_dhs = add_attr(df_meta_dhs, dict_lab, 'Lab')
    df_meta_dhs = modify_meta(df_meta_dhs)
    meta_dhs = '/local/zy/PEI/mid_data/cell_line/metadata.simple.tsv'
    df_meta_dhs.to_csv(meta_dhs, sep='\t', index=None)
    print("DHS metadata ---- completed")

    # hg38 to hg19
    path_hg38tohg19 = \
        '/local/zy/PEI/mid_data/cell_line/ENCODE/DNase-seq/GRCh38tohg19'
    hg38tohg19(path_dhs, path_hg38tohg19, meta_dhs, num_cpu)
    print("Format conversion: hg38 -> hg19 ---- completed")

    # integrate files from same experiment
    path_exp_dhs = \
        '/local/zy/PEI/mid_data/cell_line/ENCODE/' \
        'DNase-seq/GRCh38tohg19_experiment'
    merge_experiment(path_hg38tohg19, path_exp_dhs, 0.4, num_cpu)
    print("Integration of files from same experiment ---- completed")

    # build DHS reference
    path_dhs_hg38tohg19 = '/local/zy/PEI/mid_data/cell_line/DHS/GRCh38tohg19/'
    unique_bed_files(path_exp_dhs, path_dhs_hg38tohg19, 0.5, num_cpu)
    print("Integration of files from same term ---- completed")

    # standardization
    path_dhs_stan = \
        '/local/zy/PEI/mid_data/cell_line/DHS/GRCh38tohg19_standard'
    standardize_bed(path_dhs_hg38tohg19, path_dhs_stan, 'DHS', num_cpu)
    print('Standardization of DHS completed!')

    # preparation of bed files of histone and TF
    # H3K4me3
    path_h3k4me3 = \
        '/local/zy/PEI/origin_data/ENCODE/histone_ChIP-seq/H3K4me3'
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
    meta_h3k4me3 = '/local/zy/PEI/mid_data/cell_line/ENCODE/' \
                   'histone_ChIP-seq/metadata.simple.H3K4me3.tsv'
    df_meta_h3k4me3.to_csv(meta_h3k4me3, sep='\t', index=None)
    print("H3K4me3 metadata ---- completed")

    # hg38 to hg19
    path_hg38tohg19 = \
        '/local/zy/PEI/mid_data/cell_line/ENCODE/histone_ChIP-seq/H3K4me3'
    hg38tohg19(path_h3k4me3, path_hg38tohg19, meta_h3k4me3, num_cpu)
    print("Format conversion: hg38 -> hg19 ---- completed")

    # standardization
    path_h3k4me3_stan = \
        '/local/zy/PEI/mid_data/cell_line/ENCODE/histone_ChIP-seq/' \
        'H3K4me3_standard'
    standardize_bed(path_hg38tohg19, path_h3k4me3_stan, 'H3K4me3', num_cpu)
    print('Standardization of H3K4me3 completed!')

    # H3K27ac
    path_h3k27ac = \
        '/local/zy/PEI/origin_data/ENCODE/histone_ChIP-seq/H3K27ac'
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
    meta_h3k27ac = '/local/zy/PEI/mid_data/cell_line/ENCODE/' \
                   'histone_ChIP-seq/metadata.simple.H3K27ac.tsv'
    df_meta_h3k27ac.to_csv(meta_h3k27ac, sep='\t', index=None)
    print("H3K27ac metadata ---- completed")

    # hg38 to hg19
    path_hg38tohg19 = \
        '/local/zy/PEI/mid_data/cell_line/ENCODE/histone_ChIP-seq/H3K27ac'
    hg38tohg19(path_h3k27ac, path_hg38tohg19, meta_h3k27ac, num_cpu)
    print("Format conversion: hg38 -> hg19 ---- completed")

    # standardization
    path_h3k27ac_stan = \
        '/local/zy/PEI/mid_data/cell_line/ENCODE/histone_ChIP-seq/' \
        'H3K27ac_standard'
    standardize_bed(path_hg38tohg19, path_h3k27ac_stan, 'H3K27ac', num_cpu)
    print('Standardization of H3K27ac completed!')

    # CTCF
    path_ctcf = \
        '/local/zy/PEI/origin_data/ENCODE/TF_ChIP-seq/CTCF'
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
    meta_ctcf = '/local/zy/PEI/mid_data/cell_line/ENCODE/' \
                'TF_ChIP-seq/metadata.simple.CTCF.tsv'
    df_meta_ctcf.to_csv(meta_ctcf, sep='\t', index=None)
    print("CTCF metadata ---- completed")

    # hg38 to hg19
    path_hg38tohg19 = \
        '/local/zy/PEI/mid_data/cell_line/ENCODE/TF_ChIP-seq/CTCF'
    hg38tohg19(path_ctcf, path_hg38tohg19, meta_ctcf, num_cpu)
    print("Format conversion: hg38 -> hg19 ---- completed")

    # standardization
    path_ctcf_stan = \
        '/local/zy/PEI/mid_data/cell_line/ENCODE/TF_ChIP-seq/CTCF_standard'
    standardize_bed(path_hg38tohg19, path_ctcf_stan, 'CTCF', num_cpu)
    print('Standardization of CTCF completed!')

    time_end = time()
    print(time_end - time_start)
