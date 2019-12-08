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


def generate_gene_file(gtf_file, protein_file, promoter_file):
    with open(protein_file, 'w') as w_gene:
        with open(promoter_file, 'w') as w_pro:
            fmt_gene = \
                "{chrom}\t{start}\t{end}\t{symbol}\t{ensg_id}\t{strand}\n"
            fmt_promoter = \
                "{chrom}\t{start}\t{end}\t{symbol}\t{ensg_id}\t{strand}\n"
            with open(gtf_file, 'r') as r_gtf:
                for line_gene in r_gtf:
                    if line_gene[0] == '#':
                        continue
                    list_line_gene = line_gene.strip().split('\t')
                    if list_line_gene[2] != 'gene':
                        continue
                    list_attr = list_line_gene[8].strip().split('; ')
                    gene_type = list_attr[2][11:-1]
                    if gene_type != "protein_coding":
                        continue
                    gene_name = list_attr[4][11:-1]
                    ensg_id = list_attr[0][9:-1]
                    strand = list_line_gene[6]
                    dict_gene = dict(chrom=list_line_gene[0],
                                     start=list_line_gene[3],
                                     end=list_line_gene[4], symbol=gene_name,
                                     ensg_id=ensg_id,
                                     strand=strand)
                    w_gene.write(fmt_gene.format(**dict_gene))
                    if strand == '+':
                        pro_start = str(int(list_line_gene[3]) - 2000)
                        pro_end = str(int(list_line_gene[3]) + 2000)
                    elif strand == '-':
                        pro_start = str(int(list_line_gene[4]) - 2000)
                        pro_end = str(int(list_line_gene[4]) + 2000)
                    else:
                        print('Error')
                        break
                    dict_promoter = dict(chrom=list_line_gene[0],
                                         start=pro_start,
                                         end=pro_end,
                                         symbol=f"{gene_name}<-"
                                                f"{list_line_gene[0]}:"
                                                f"{pro_start}-{pro_end}",
                                         ensg_id=ensg_id,
                                         strand=strand)
                    w_pro.write(fmt_promoter.format(**dict_promoter))

    return


def filter_meta(meta_in):
    # reserve released data and original(untreated) peak file
    # simplify meta file
    df_meta = pd.read_csv(meta_in, sep='\t')
    df_meta_released = df_meta.loc[
        (df_meta['File Status'] == 'released') &
        (df_meta['Output type'] == 'peaks') &
        (df_meta['Biological replicate(s)'].apply(
            lambda x: len(str(x).split(', ')) == 1
        )),
        ['File accession', 'Experiment accession', 'Biosample term id',
         'Biosample term name', 'Biosample type', 'Biosample treatments',
         'Biosample genetic modifications methods', 'Assembly',
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


def modify_meta(df_meta, set_ref, df_com):
    # only select tissue data
    df_meta = df_meta.loc[df_meta['Biosample type'] == 'tissue', :]
    rows1 = df_meta.shape[0]

    # reference organs
    df_meta_nan = df_meta.loc[df_meta['Biosample organ'] == '', :]
    df_meta_nan = df_meta_nan.drop(['Biosample organ'], 1)
    df_meta_nan = pd.merge(df_meta_nan, df_com,
                           on=['Biosample term id', 'Biosample term name'])
    df_meta = df_meta.loc[df_meta['Biosample organ'] != '', :]
    df_meta = pd.concat([df_meta, df_meta_nan], sort=False)
    rows2 = df_meta.shape[0]
    assert rows1 == rows2
    organs = df_meta['Biosample organ'].tolist()
    new_organs = []
    for organ in organs:
        list_organ = \
            [val for val in organ.strip().split(',') if val in set_ref]
        new_organs.append(','.join(list_organ))

    df_meta = df_meta.drop('Biosample organ', 1)
    df_meta['Biosample organ'] = new_organs

    # correct fuzzy life stage labels to 'unknown'
    life_stages = df_meta['Biosample life stage'].tolist()
    new_life_stages = []
    for life_stage in life_stages:
        if life_stage in ['adult', 'embryonic', 'unknown', 'newborn', 'child']:
            new_life_stages.append(life_stage)
        else:
            new_life_stages.append('unknown')
    df_meta = df_meta.drop('Biosample life stage', 1)
    df_meta['Biosample life stage'] = new_life_stages

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
    if dict_in['Assembly'] == 'GRCh38':
        file_prefix = file_hg38 + '.prefix'
        file_suffix = file_hg38 + '.suffix'
        file_hg19_prefix = file_hg19 + '.prefix'
        file_hg19_format = file_hg19 + '.format'
        os.system(f"cut -f 1,2,3,4 {file_hg38} > {file_prefix}")
        os.system(f"cut -f 4,5,6,7,8,9,10 {file_hg38} > {file_suffix}")
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
        up_limit = np.max(length_old) + 10
        down_limit = np.min(length_old) - 10
        df_bed = pd.read_csv(file_hg19_format, sep='\t', header=None)
        length = df_bed.iloc[:, 2] - df_bed.iloc[:, 1]
        df_bed = df_bed.loc[(length < up_limit) & (length > down_limit), :]
        df_bed.to_csv(file_hg19, sep='\t', index=None, header=None)
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


def hg38tohg19(path_hg38, path_hg19, num_process):
    if os.path.exists(path_hg19):
        os.system(f"rm -rf {path_hg19}")
    os.mkdir(path_hg19)

    df_meta = pd.read_csv(
        os.path.join(path_hg38, 'metadata.simple.tsv'), sep='\t'
    )
    list_meta = []
    experiments = set(df_meta['Experiment accession'].tolist())
    for exp in experiments:
        df_exp = df_meta.loc[df_meta['Experiment accession'] == exp, :]
        hg19 = df_exp.loc[df_exp['Assembly'] == 'hg19', :]
        hg38 = df_exp.loc[df_exp['Assembly'] == 'GRCh38', :]
        if hg19.shape[0] == hg38.shape[0]:
            list_meta.append(hg19)
        else:
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

    os.system(f"bedtools sort -i {split_out} > {final_out}")

    os.remove(cat_out)
    os.remove(sort_out)
    os.remove(merge_out)
    os.remove(split_out)

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
            len_list1 = (df_ref.loc[df_ref.loc[:, com[0]] != 0, :]).shape[0]
            len_list2 = (df_ref.loc[df_ref.loc[:, com[1]] != 0, :]).shape[0]
            total = min(len_list1, len_list2)
            len_overlap = (df_ref.loc[
                           (df_ref.loc[:, com[0]] != 0) &
                           (df_ref.loc[:, com[1]] != 0), :]).shape[0]
            list_out.append({'Name': term_name, 'Combination': com,
                             'Overlap ratio': len_overlap/total})

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
                   'Biosample organ']]
    meta_out = meta_out.drop_duplicates()
    dict_meta = []
    for sub_dict in meta_out.to_dict('records'):
        organs = sub_dict['Biosample organ'].split(',')
        for organ in organs:
            dict_meta.append(
                {'Biosample term name': sub_dict['Biosample term name'],
                 'Biosample life stage': sub_dict['Biosample life stage'],
                 'Biosample organ': organ})
    meta_out = pd.DataFrame(dict_meta)
    meta_out.to_csv(os.path.join(path_out, 'meta.reference.tsv'),
                    sep='\t', index=None)

    list_input = []
    organs = []
    for line in df_meta['Biosample organ'].tolist():
        organs.extend(line.strip().split(','))
    organs = set(organs)
    for organ in organs:
        organ_meta = df_meta.loc[
                     df_meta['Biosample organ'].apply(
                         lambda x: organ in x.strip().split(',')), :]
        organ_path = \
            os.path.join(path_out, organ.replace(' ', '_'))
        if not os.path.exists(organ_path):
            os.makedirs(organ_path)
        life_stages = []
        for line in organ_meta['Biosample life stage'].tolist():
            life_stages.extend(line.strip().split(','))
        life_stages = set(life_stages)
        for life_stage in life_stages:
            life_meta = \
                organ_meta.loc[
                 organ_meta['Biosample life stage'].apply(
                     lambda x: life_stage in x.strip().split(',')), :]
            path_life_stage = \
                os.path.join(organ_path, life_stage.replace(' ', '_'))
            if not os.path.exists(path_life_stage):
                os.makedirs(path_life_stage)
            terms = set(life_meta['Biosample term name'].tolist())
            for term in terms:
                term_meta = \
                    life_meta.loc[life_meta['Biosample term name'] == term, :]
                accession_ids = \
                    list(set(term_meta['Experiment accession'].tolist()))
                path_term = \
                    os.path.join(path_life_stage, term.replace(
                        ' ', '_').replace('/', '+').replace("'", '--'))
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
        organ = dict_in['Biosample organ']
        life_stage = dict_in['Biosample life stage']
        term = dict_in['Biosample term name']
        term_name = term.replace(' ', '_').replace('/', '+').replace("'", "--")
        file_label = f"{organ}|{life_stage}|{term}"
        file = f"{organ.replace(' ', '_')}/{life_stage.replace(' ', '_')}/" \
               f"{term_name}/{term_name}.bed"
        folder1 = os.path.join(path_out, f"{organ.replace(' ', '_')}")
        folder2 = os.path.join(
            path_out,
            f"{organ.replace(' ', '_')}/{life_stage.replace(' ', '_')}"
        )
        folder3 = os.path.join(
            path_out,
            f"{organ.replace(' ', '_')}/{life_stage.replace(' ', '_')}/"
            f"{term_name}"
        )
        file_in = os.path.join(path_in, file)
        file_out = os.path.join(path_out, file)
        # make folder
        if not os.path.exists(folder1):
            os.mkdir(folder1)
        if not os.path.exists(folder2):
            os.mkdir(folder2)
        if not os.path.exists(folder3):
            os.mkdir(folder3)
    else:
        file_in = os.path.join(path_in, dict_in['File accession'] + '.bed')
        file_out = os.path.join(path_out, dict_in['File accession'] + '.bed')

    with open(file_in, 'r') as r_f:
        with open(file_out, 'w') as w_f:
            fmt_dhs = "{chrom}\t{start}\t{end}\t{label}\t.\t.\t{file_label}\n"
            fmt_histone = "{chrom}\t{start}\t{end}\t{label}\t" \
                          "{score}\t.\t{file_label}\t{accessions}\n"
            for line in r_f:
                list_line = line.strip().split('\t')
                chrom = list_line[0]
                start = list_line[1]
                end = list_line[2]
                accessions = list_line[3]
                label = f"{type_bed}<-{chrom}:{start}-{end}"
                if type_bed == 'DHS':
                    w_f.write(fmt_dhs.format(**locals()))
                else:
                    score = list_line[7]
                    w_f.write(fmt_histone.format(**locals()))

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

    if type_bed == 'DHS':
        df_meta = pd.read_csv(
            os.path.join(path_in, 'meta.reference.tsv'), sep='\t'
        )
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


def merge_standard_bed(path_bed, dict_in):
    flank_percent = dict_in['flank_percent']
    term_name = dict_in['term_name']
    path_out = dict_in['path']
    col_collapse = '2,3,5,7,8'
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
            fmt = "{chrom}\t{start}\t{end}\t.\t{p_value}\t.\t{label}\t" \
                  "{accessions}\n"
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
                array_accessions = np.array(
                    [num for num in list_line[7].strip().split(',')])
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
                    select_accessions = array_accessions[select_bool]
                    # reset arrays
                    array_start = array_start[~select_bool]
                    array_end = array_end[~select_bool]
                    array_label = array_label[~select_bool]
                    array_p_value = array_p_value[~select_bool]
                    array_accessions = array_accessions[~select_bool]
                    array_length = array_end - array_start
                    # write new rows
                    start = str(np.min(select_start))
                    end = str(np.max(select_end))
                    label = \
                        '/'.join(map(str, select_label.tolist()))
                    p_value = '/'.join(map(str, select_p_value.tolist()))
                    accessions = '|'.join(map(str, select_accessions.tolist()))
                    w_f.write(fmt.format(**locals()))

    os.remove(cat_out)
    os.remove(sort_out)

    return


def sub_merge(dict_in):
    sub_ref_meta = dict_in['sub_ref_meta']
    sub_meta = dict_in['sub_meta']
    sub_path_in = dict_in['path_in']
    sub_path_out = dict_in['path_out']
    organ = dict_in['organ']
    file_out = os.path.join(sub_path_out, organ.replace(' ', '_') + '.bed')
    if not os.path.exists(sub_path_out):
        os.mkdir(sub_path_out)
    list_meta = sub_ref_meta.to_dict("records")
    accession_ids = [
        '/'.join([sub_dict['Biosample life stage'].replace(' ', '_'),
                  sub_dict['Biosample term name'].replace(
                        ' ', '_').replace('/', '+').replace("'", '--'),
                  sub_dict['Biosample term name'].replace(
                        ' ', '_').replace('/', '+').replace("'", '--')])
        for sub_dict in list_meta]
    if sub_ref_meta.shape[0] == 1:
        os.system(
            f"cp {os.path.join(sub_path_in, accession_ids[0] + '.bed')} "
            f"{file_out}"
        )
        with open(file_out, 'w') as w_f:
            fmt = "{chrom}\t{start}\t{end}\t.\t{score}\t.\t{label}\n"
            with open(os.path.join(sub_path_in, accession_ids[0] + '.bed'),
                      'r') as r_f:
                for line in r_f:
                    list_line = line.strip().split('\t')
                    dict_out = dict(
                        chrom=list_line[0], start=list_line[1],
                        end=list_line[2], label=list_line[6],
                        score=list_line[4]
                    )
                    w_f.write(fmt.format(**dict_out))
        return
    elif sub_ref_meta.shape[0] > 1:
        dict_merge = dict(
            path=sub_path_out,
            term_name=organ.replace(' ', '_'),
            accession_ids=accession_ids,
            flank_percent=0.6)
        merge_standard_bed(sub_path_in, dict_merge)
        labels = [f"{organ}|{sub_dict['Biosample life stage']}|"
                  f"{sub_dict['Biosample term name']}"
                  for sub_dict in list_meta]
        accessions = [sub_dict['File accession']
                      for sub_dict in sub_meta.to_dict("records")]
        list_bed = \
            (pd.read_csv(file_out, sep='\t', header=None)).to_dict('records')
        list_label = []
        list_exp = []
        for sub_dict in list_bed:
            dict_label = dict()
            dict_exp = dict()
            dict_label['peak_id'] = \
                f"{sub_dict[0]}:{str(sub_dict[1])}-{str(sub_dict[2])}"
            dict_exp['peak_id'] = \
                f"{sub_dict[0]}:{str(sub_dict[1])}-{str(sub_dict[2])}"
            set_label = set(sub_dict[6].strip().split('/'))
            for label in labels:
                if label in set_label:
                    dict_label[label] = 1
                else:
                    dict_label[label] = 0
            set_accession = set(sub_dict[7].strip().split('|'))
            for accession in accessions:
                if accession in set_accession:
                    dict_exp[accession] = 1
                else:
                    dict_exp[accession] = 0
            list_label.append(dict_label)
            list_exp.append(dict_exp)

        # label matrix
        df_label = pd.DataFrame(list_label, columns=['peak_id'] + labels)
        df_label.index = df_label['peak_id']
        df_label = df_label.drop('peak_id', 1)
        # write score matrix to txt file
        mat_label = os.path.join(sub_path_out, 'label_matrix.txt')
        df_label.to_csv(mat_label, sep='\t')

        # experiment matrix
        df_exp = \
            pd.DataFrame(list_exp, columns=['peak_id'] + accessions)
        df_exp.index = df_exp['peak_id']
        df_exp = df_exp.drop('peak_id', 1)
        # write score matrix to txt file
        mat_exp = os.path.join(sub_path_out, 'experiment_matrix.txt')
        df_exp.to_csv(mat_exp, sep='\t')

        # calculate overlap ratio
        list_out = []
        list_com = combinations(df_label.columns, 2)
        for com in list_com:
            len_list1 = (df_label.loc[
                         df_label.loc[:, com[0]] != 0, :]).shape[0]
            len_list2 = (df_label.loc[
                         df_label.loc[:, com[1]] != 0, :]).shape[0]
            total = min(len_list1, len_list2)
            len_overlap = (df_label.loc[
                           (df_label.loc[:, com[0]] != 0) &
                           (df_label.loc[:, com[1]] != 0), :]).shape[0]
            list_out.append({'Name': organ, 'Combination': com,
                             'Overlap ratio': len_overlap/total})

        df_out = pd.DataFrame(list_out)

        # scatter plot
        # os.system(f"Rscript scatter.plot.organ.R {mat_label} "
        #           f"{os.path.join(sub_path_out, 'scatter.pdf')}")
        # os.system(f"Rscript scatter.plot.organ.accession.R {mat_label} "
        #           f"{os.path.join(sub_path_out, 'scatter.pdf')}")

        return df_out


def merge_organ_cluster(path_in, path_out, num_process):
    if os.path.exists(path_out):
        os.system(f"rm -rf {path_out}")
    os.mkdir(path_out)
    os.system(f"cp {os.path.join(path_in, 'meta.reference.tsv')} "
              f"{os.path.join(path_out, 'meta.reference.tsv')}")
    os.system(f"cp {os.path.join(path_in, 'metadata.simple.tsv')} "
              f"{os.path.join(path_out, 'metadata.simple.tsv')}")

    df_meta_ref = \
        pd.read_csv(os.path.join(path_in, 'meta.reference.tsv'), sep='\t')
    df_meta = \
        pd.read_csv(os.path.join(path_in, 'metadata.simple.tsv'), sep='\t')
    organs = []
    for line in df_meta_ref['Biosample organ'].tolist():
        organs.extend(line.strip().split(','))
    organs = set(organs)
    list_input = []
    for organ in organs:
        sub_ref_meta = df_meta_ref.loc[
                   df_meta_ref['Biosample organ'].apply(
                       lambda x: organ in x.strip().split(',')), :]
        sub_meta = df_meta.loc[
                   df_meta['Biosample organ'].apply(
                       lambda x: organ in x.strip().split(',')), :]
        list_input.append(
            dict(sub_ref_meta=sub_ref_meta, organ=organ, sub_meta=sub_meta,
                 path_out=os.path.join(path_out, organ.replace(' ', '_')),
                 path_in=os.path.join(path_in, organ.replace(' ', '_'))))

    pool = Pool(processes=num_process)
    list_df = pool.map(sub_merge, list_input)
    pool.close()

    df_overlap = pd.concat(list_df, sort=False)
    df_overlap.to_csv(os.path.join(path_out, 'overlap.txt'), sep='\t')

    # merge all organ
    organs = os.listdir(path_out)
    accession_ids = []
    labels = []
    for organ in organs:
        path_organ = os.path.join(path_out, organ)
        if os.path.isdir(path_organ):
            labels.append(organ.replace('_', ' '))
            accession_ids.append(os.path.join(path_organ, organ))

    dict_merge = dict(
        path=path_out,
        term_name='all_organs',
        accession_ids=accession_ids,
        flank_percent=0.7)
    merge_standard_bed(path_out, dict_merge)
    list_bed = \
        (pd.read_csv(os.path.join(path_out, 'all_organs.bed'),
                     sep='\t', header=None)).to_dict('records')
    list_dict = []
    for sub_dict in list_bed:
        out_dict = dict()
        out_dict['peak_id'] = \
            f"{sub_dict[0]}:{str(sub_dict[1])}-{str(sub_dict[2])}"
        list_label = set([val.split('|')[0]
                          for val in sub_dict[6].strip().split('/')])
        for label in labels:
            if label in list_label:
                out_dict[label] = 1
            else:
                out_dict[label] = 0
        list_dict.append(out_dict)
    df_label_peak = pd.DataFrame(list_dict, columns=['peak_id'] + labels)
    df_label_peak.index = df_label_peak['peak_id']
    df_label_peak = df_label_peak.drop('peak_id', 1)
    # write score matrix to txt file
    mat_peak = os.path.join(path_out, 'all_organs_label_peak.txt')
    df_label_peak.to_csv(mat_peak, sep='\t')
    # calculate overlap ratio
    list_out = []
    list_com = combinations(df_label_peak.columns, 2)
    for com in list_com:
        len_list1 = (df_label_peak.loc[
                     df_label_peak.loc[:, com[0]] != 0, :]).shape[0]
        len_list2 = (df_label_peak.loc[
                     df_label_peak.loc[:, com[1]] != 0, :]).shape[0]
        total = (len_list1 + len_list2) / 2
        len_overlap = (df_label_peak.loc[
                       (df_label_peak.loc[:, com[0]] != 0) &
                       (df_label_peak.loc[:, com[1]] != 0), :]).shape[0]
        list_out.append({'Combination': com,
                         'Overlap ratio': len_overlap / total})

    df_overlap = pd.DataFrame(list_out)
    df_overlap.to_csv(
        os.path.join(path_out, 'all_organs_overlap.txt'), sep='\t', index=None
    )

    # scatter plot
    # os.system(f"Rscript scatter.plot.R {mat_peak} "
    #           f"{os.path.join(path_out, 'all_organs_scatter.pdf')}")

    return


if __name__ == '__main__':
    time_start = time()
    # parameters
    num_cpu = 40
    # get bed file annotating protein-coding genes
    gtf_file_hg19 = \
        '/local/zy/PEI/data/ENCODE/gencode.v19.annotation.gtf'
    protein_file_hg19 = \
        '/local/zy/PEI/data/gene/genes.protein.gencode.v19.bed'
    promoter_file_hg19 = \
        '/local/zy/PEI/data/gene/' \
        'promoters.up2k.protein.gencode.v19.bed'
    generate_gene_file(gtf_file_hg19, protein_file_hg19, promoter_file_hg19)

    # build life stage dictionary
    path_lifestage = '/local/zy/PEI/data/ENCODE/metadata/life_stage'
    dict_lifestage = build_dict_attr(path_lifestage)

    # build organ dictionary
    path_meta_organ = '/local/zy/PEI/data/ENCODE/metadata/organ'
    dict_organ = build_dict_attr(path_meta_organ)

    # build organ dictionary
    path_meta_cell = '/local/zy/PEI/data/ENCODE/metadata/cell'
    dict_cell = build_dict_attr(path_meta_cell)

    # build organ dictionary
    path_meta_lab = '/local/zy/PEI/data/ENCODE/metadata/lab'
    dict_lab = build_dict_attr(path_meta_lab)

    # read reference organ
    ref_organ = '/local/zy/PEI/data/ENCODE/metadata/organ_ref.txt'
    with open(ref_organ, 'r') as r_ref:
        set_organs = set([organ.strip() for organ in r_ref])
    # organ complement
    file_complement = \
        '/local/zy/PEI/data/ENCODE/metadata/complement_organ.txt'
    df_complement = pd.read_csv(file_complement, sep='\t')
    print("Preparation of dictionary files and reference files is completed")

    # DHS reference
    # metafile
    path_dhs = '/local/zy/PEI/data/ENCODE/DNase-seq/all'
    ori_meta_dhs = os.path.join(path_dhs, 'metadata.tsv')
    df_meta_dhs = filter_meta(ori_meta_dhs)
    df_meta_dhs = add_attr(df_meta_dhs, dict_lifestage, 'Biosample life stage')
    df_meta_dhs = add_attr(df_meta_dhs, dict_organ, 'Biosample organ')
    df_meta_dhs = add_attr(df_meta_dhs, dict_cell, 'Biosample cell')
    df_meta_dhs = add_attr(df_meta_dhs, dict_lab, 'Lab')
    df_meta_dhs = modify_meta(df_meta_dhs, set_organs, df_complement)
    meta_dhs = os.path.join(path_dhs, 'metadata.simple.tsv')
    df_meta_dhs.to_csv(meta_dhs, sep='\t', index=None)
    print("DHS metadata ---- completed")

    # hg38 to hg19
    path_hg38tohg19 = \
        '/local/zy/PEI/data/ENCODE/DNase-seq/GRCh38tohg19'
    hg38tohg19(path_dhs, path_hg38tohg19, num_cpu)
    print("Format conversion: hg38 -> hg19 ---- completed")

    # integrate files from same experiment
    path_exp_dhs = \
        '/local/zy/PEI/data/ENCODE/DNase-seq/GRCh38tohg19_experiment'
    merge_experiment(path_hg38tohg19, path_exp_dhs, 0.4, num_cpu)
    print("Integration of files from same experiment ---- completed")

    # build DHS reference
    path_dhs_hg38tohg19 = '/local/zy/PEI/data/DHS/GRCh38tohg19/'
    unique_bed_files(path_exp_dhs, path_dhs_hg38tohg19, 0.5, num_cpu)
    print("Integration of files from same term ---- completed")

    # standardization
    path_dhs_stan = '/local/zy/PEI/data/DHS/GRCh38tohg19_standard'
    standardize_bed(path_dhs_hg38tohg19, path_dhs_stan, 'DHS', num_cpu)
    print('Standardization of DHS completed!')

    # merge and cluster
    path_dhs_cluster = '/local/zy/PEI/data/DHS/GRCh38tohg19_cluster'
    merge_organ_cluster(path_dhs_stan, path_dhs_cluster, num_cpu)
    print('Cluster and merge of DHS completed!')

    # preparation of bed files of histone and TF
    # H3K4me3
    path_h3k4me3 = \
        '/local/zy/PEI/data/ENCODE/histone_ChIP-seq/H3K4me3'
    ori_meta_h3k4me3 = os.path.join(path_h3k4me3, 'metadata.tsv')
    df_meta_h3k4me3 = filter_meta(ori_meta_h3k4me3)
    df_meta_h3k4me3 = \
        add_attr(df_meta_h3k4me3, dict_lifestage, 'Biosample life stage')
    df_meta_h3k4me3 = add_attr(df_meta_h3k4me3, dict_organ, 'Biosample organ')
    df_meta_h3k4me3 = add_attr(df_meta_h3k4me3, dict_cell, 'Biosample cell')
    df_meta_h3k4me3 = add_attr(df_meta_h3k4me3, dict_lab, 'Lab')
    df_meta_h3k4me3 = modify_meta(df_meta_h3k4me3, set_organs, df_complement)
    meta_h3k4me3 = os.path.join(path_h3k4me3, 'metadata.simple.tsv')
    df_meta_h3k4me3.to_csv(meta_h3k4me3, sep='\t', index=None)

    # hg38 to hg19
    path_hg38tohg19 = \
        '/local/zy/PEI/data/ENCODE/histone_ChIP-seq/' \
        'GRCh38tohg19/H3K4me3'
    hg38tohg19(path_h3k4me3, path_hg38tohg19, num_cpu)

    # H3K27ac
    path_h3k27ac = \
        '/local/zy/PEI/data/ENCODE/histone_ChIP-seq/H3K27ac'
    ori_meta_h3k27ac = os.path.join(path_h3k27ac, 'metadata.tsv')
    df_meta_h3k27ac = filter_meta(ori_meta_h3k27ac)
    df_meta_h3k27ac = \
        add_attr(df_meta_h3k27ac, dict_lifestage, 'Biosample life stage')
    df_meta_h3k27ac = add_attr(df_meta_h3k27ac, dict_organ, 'Biosample organ')
    df_meta_h3k27ac = add_attr(df_meta_h3k27ac, dict_cell, 'Biosample cell')
    df_meta_h3k27ac = add_attr(df_meta_h3k27ac, dict_lab, 'Lab')
    df_meta_h3k27ac = modify_meta(df_meta_h3k27ac, set_organs, df_complement)
    meta_h3k27ac = os.path.join(path_h3k27ac, 'metadata.simple.tsv')
    df_meta_h3k27ac.to_csv(meta_h3k27ac, sep='\t', index=None)

    # hg38 to hg19
    path_hg38tohg19 = \
        '/local/zy/PEI/data/ENCODE/histone_ChIP-seq/' \
        'GRCh38tohg19/H3K27ac'
    hg38tohg19(path_h3k27ac, path_hg38tohg19, num_cpu)

    time_end = time()
    print(time_end - time_start)
