#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: preparation.py
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
import glob
import sys
# sys.path.append('/local/zy/my_git/bioinformatics/PEI/annotate_cRE')
# file_chain = \
#     '/local/zy/tools/files_liftOver/hg38ToHg19.over.chain.gz'
# liftover = '/local/zy/tools/liftOver'
file_chain = '/lustre/tianlab/tools/files_liftOver/hg38ToHg19.over.chain.gz'
liftover = '/lustre/tianlab/tools/liftOver'
root_path = '/local/tianlab/zhangyu/my_git/bioinformatics/PEI/annotate_cRE'


def generate_gene_file(gtf_file, protein_file, promoter_file, promoter_merge,
                       exon_file):
    exon_file_tmp = exon_file + '.tmp'
    with open(protein_file, 'w') as w_gene:
        with open(promoter_file, 'w') as w_pro:
            with open(exon_file_tmp, 'w') as w_exon:
                fmt_gene = \
                    "{chrom}\t{start}\t{end}\t{symbol}\t{ensg_id}\t{strand}\n"
                with open(gtf_file, 'r') as r_gtf:
                    for line_gene in r_gtf:
                        if line_gene[0] == '#':
                            continue
                        list_line_gene = line_gene.strip().split('\t')
                        list_attr = list_line_gene[8].strip().split('; ')
                        gene_name = list_attr[4][11:-1]
                        ensg_id = list_attr[0][9:-1]
                        strand = list_line_gene[6]
                        gene_type = list_attr[2][11:-1]
                        if gene_type != "protein_coding":
                            continue
                        if list_line_gene[2] == 'gene':
                            dict_gene = dict(chrom=list_line_gene[0],
                                             start=list_line_gene[3],
                                             end=list_line_gene[4],
                                             symbol=gene_name,
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
                                                        f"{pro_start}-"
                                                        f"{pro_end}",
                                                 ensg_id=ensg_id,
                                                 strand=strand)
                            w_pro.write(fmt_gene.format(**dict_promoter))
                        elif list_line_gene[2] == 'exon':
                            dict_exon = dict(chrom=list_line_gene[0],
                                             start=list_line_gene[3],
                                             end=list_line_gene[4],
                                             symbol=gene_name,
                                             ensg_id=ensg_id,
                                             strand=strand)
                            w_exon.write(fmt_gene.format(**dict_exon))

    promoter_sort = promoter_file + '.sort'
    os.system(f"bedtools sort -i {promoter_file} > {promoter_sort}")
    os.system(f"bedtools merge -i {promoter_sort} "
              f"-c 4,5,6 -o collapse,collapse,collapse > {promoter_merge}")
    os.system(f"bedtools sort -i {exon_file_tmp} > {exon_file}")
    os.remove(promoter_sort)
    os.remove(exon_file_tmp)

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


def modify_meta(df_meta, set_ref, df_com):
    # only select tissue data
    df_meta = df_meta.loc[df_meta['Biosample type'] == 'tissue', :]
    # rows1 = df_meta.shape[0]

    # reference organs
    df_meta_nan = df_meta.loc[df_meta['Biosample organ'] == '', :]
    df_meta_nan = df_meta_nan.drop(['Biosample organ'], 1)
    df_meta_nan = pd.merge(df_meta_nan, df_com,
                           on=['Biosample term id', 'Biosample term name'])
    df_meta = df_meta.loc[df_meta['Biosample organ'] != '', :]
    df_meta = pd.concat([df_meta, df_meta_nan], sort=False)
    # rows2 = df_meta.shape[0]
    # assert rows1 == rows2
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
    df_meta.loc[df_meta['Biosample organ'] == 'extraembryonic component',
                'Biosample life stage'] = 'embryonic'
    df_meta = df_meta.loc[df_meta['Biosample life stage'] != 'unknown', :]

    # delete 'embryo'
    df_meta = df_meta.loc[df_meta['Biosample organ'] != 'embryo', :]

    # combine life and organ
    # df_meta['Biosample life_organ'] = df_meta.apply(
    #     lambda x: x['Biosample life stage'] + '_' + x['Biosample organ'],
    #     axis=1
    # )

    return df_meta


def sub_hg38tohg19(path_hg38, path_hg19, dict_in):
    file_hg38 = os.path.join(path_hg38, dict_in['File accession'] + '.bed')
    file_hg19_unsort = \
        os.path.join(path_hg19, dict_in['File accession'] + '.unsort.bed')
    file_hg19 = os.path.join(path_hg19, dict_in['File accession'] + '.bed')
    set_chroms = {'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
                  'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                  'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
                  'chr20', 'chr21', 'chr22', 'chrX', 'chrY'}
    # Download from UCSC
    file_ummap = os.path.join(
        path_hg19, dict_in['File accession'] + '.bed.unmap')
    if dict_in['Assembly'] == 'hg19':
        with open(file_hg19_unsort, 'w') as w_f:
            fmt = "{chrom}\t{start}\t{end}\t{peak_id}\t{score}\t{strand}\t" \
                  "{fold_change}\t{p_value}\t{q_value}\t{peak_location}\n"
            with open(file_hg38, 'r') as r_hg38:
                for line in r_hg38:
                    list_line = line.strip().split('\t')
                    chrom = list_line[0]
                    if chrom not in set_chroms:
                        continue
                    dict_hg19 = dict(
                        chrom=list_line[0], start=list_line[1],
                        end=list_line[2], peak_id=dict_in['File accession'],
                        score=list_line[4], strand=list_line[5],
                        fold_change=list_line[6], p_value=list_line[7],
                        q_value=list_line[8], peak_location=list_line[9]
                    )
                    w_f.write(fmt.format(**dict_hg19))
        df_bed = pd.read_csv(file_hg19_unsort, sep='\t', header=None)
        # scores = np.array(df_bed.iloc[:, 6]).reshape(-1, 1)
        #
        # # scale_scores = StandardScaler().fit_transform(scores)
        # scale_scores = scores / np.mean(scores)
        #
        # df_bed.iloc[:, 6] = scale_scores
        df_bed = df_bed.drop_duplicates([0, 1, 2])
        df_bed.to_csv(file_hg19_unsort, sep='\t', index=None, header=None)

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
        os.system(f"{liftover} {file_prefix} {file_chain} "
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
                    chrom = list_line[0]
                    if chrom not in set_chroms:
                        continue
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
        # scores = np.array(df_bed.iloc[:, 6]).reshape(-1, 1)
        #
        # # scale_scores = StandardScaler().fit_transform(scores)
        #
        # scale_scores = scores / np.mean(scores)
        #
        # df_bed.iloc[:, 6] = scale_scores
        # if dict_in['File accession'] == 'ENCFF742USA':
        #     print(df_bed.shape[0])
        df_bed = df_bed.drop_duplicates([0, 1, 2])
        # if dict_in['File accession'] == 'ENCFF742USA':
        #     print(df_bed.shape[0])
        df_bed.to_csv(file_hg19_unsort, sep='\t', index=None, header=None)

        os.remove(file_hg38_labeled)
        os.remove(file_ummap)
        os.remove(file_prefix)
        os.remove(file_suffix)
        os.remove(file_hg19_prefix)
        os.remove(file_hg19_format)

    os.system(f"bedtools sort -i {file_hg19_unsort} > {file_hg19}")
    os.remove(file_hg19_unsort)

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
    os.system(f"cp {path_meta} "
              f"{os.path.join(path_hg19, 'metadata.simple.tsv')}")

    df_meta = pd.read_csv(path_meta, sep='\t')
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


def split_merge_bed(sub_path_out, flank_percent, subfile):
    # split merge file
    _, subfile_name = os.path.split(subfile)
    sub_split_out = os.path.join(sub_path_out, f"{subfile_name}.bed.unsort")
    sub_split_sort_out = \
        os.path.join(sub_path_out, f"{subfile_name}.bed.split.sort")
    sub_final_out = os.path.join(sub_path_out, f"{subfile_name}.bed")
    with open(subfile, 'r') as r_f:
        with open(sub_split_out, 'w') as w_f:
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

    os.system(f"bedtools sort -i {sub_split_out} > {sub_split_sort_out}")
    with open(sub_final_out, 'w') as w_final:
        old_end = 0
        old_chrom = '0'
        with open(sub_split_sort_out, 'r') as r_sort:
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

    df_final = pd.read_csv(sub_final_out, sep='\t', header=None,
                           names=['chrom', 'start', 'end', 'access',
                                  'score', 'strand', 'fold_change',
                                  'p_value', 'q_value', 'center'])
    len_df = df_final.shape[0]
    region_length = df_final['end'] - df_final['start']
    df_narrow = df_final.loc[region_length <= 80, :]
    set_index = set(df_narrow.index)
    for i in df_narrow.index:
        if (i - 1) in set_index:
            continue
        if i + 1 == len_df:
            if df_final.loc[i, 'start'] == df_final.loc[i-1, 'end'] + 1:
                df_final.loc[i-1, 'end'] = df_final.loc[i, 'end']
                df_final.loc[i-1, 'access'] = \
                    '|'.join((df_final.loc[i-1:i+1, 'access']).tolist())
                df_final.loc[i-1, 'fold_change'] = \
                    '|'.join((df_final.loc[i-1:i+1, 'fold_change']).tolist())
                df_final.loc[i-1, 'p_value'] = \
                    '|'.join((df_final.loc[i-1:i+1, 'p_value']).tolist())
                df_final = df_final.drop(i)
            continue
        # try:
        #     a = (df_final.loc[i, 'end'] == df_final.loc[i+1, 'start'] - 1)
        # except KeyError:
        #     print(sub_final_out)
        #     print(i + 1 == df_final.shape[0])
        #     print(i + 1)
        #     print(df_final.shape[0])
        if i == 0:
            if df_final.loc[i, 'end'] == df_final.loc[i+1, 'start'] - 1:
                df_final.loc[i+1, 'start'] = df_final.loc[i, 'start']
                df_final.loc[i+1, 'access'] = \
                    '|'.join((df_final.loc[i:i+2, 'access']).tolist())
                df_final.loc[i+1, 'fold_change'] = \
                    '|'.join((df_final.loc[i:i+2, 'fold_change']).tolist())
                df_final.loc[i+1, 'p_value'] = \
                    '|'.join((df_final.loc[i:i+2, 'p_value']).tolist())
                df_final = df_final.drop(i)
            continue
        if df_final.loc[i, 'end'] == df_final.loc[i+1, 'start'] - 1:
            df_final.loc[i+1, 'start'] = df_final.loc[i, 'start']
            df_final.loc[i+1, 'access'] = \
                '|'.join((df_final.loc[i:i+2, 'access']).tolist())
            df_final.loc[i+1, 'fold_change'] = \
                '|'.join((df_final.loc[i:i+2, 'fold_change']).tolist())
            df_final.loc[i+1, 'p_value'] = \
                '|'.join((df_final.loc[i:i+2, 'p_value']).tolist())
            df_final = df_final.drop(i)
            continue
        elif df_final.loc[i, 'start'] == df_final.loc[i-1, 'end'] + 1:
            df_final.loc[i-1, 'end'] = df_final.loc[i, 'end']
            df_final.loc[i-1, 'access'] = \
                '|'.join((df_final.loc[i-1:i+1, 'access']).tolist())
            df_final.loc[i-1, 'fold_change'] = \
                '|'.join((df_final.loc[i-1:i+1, 'fold_change']).tolist())
            df_final.loc[i-1, 'p_value'] = \
                '|'.join((df_final.loc[i-1:i+1, 'p_value']).tolist())
            df_final = df_final.drop(i)

    df_final.to_csv(sub_final_out, sep='\t', index=None, header=None)

    return


def pre_merge_bed(path_bed, dict_in):
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

    os.remove(cat_out)
    os.remove(sort_out)

    return


def concat_subfiles(dict_in):
    term_name = dict_in['term_name']
    sub_path_out = dict_in['path']
    path_subfile_in = dict_in['path_subfile_in']
    path_subfile_out = dict_in['path_subfile_out']
    split_out = os.path.join(sub_path_out, f"{term_name}.bed.unsort")
    split_sort_out = os.path.join(sub_path_out, f"{term_name}.bed")
    cat_in = ' '.join([subfile_out for subfile_out in
                       glob.glob(path_subfile_out + '/*.bed')])
    os.system(f"cat {cat_in} > {split_out}")
    os.system(f"bedtools sort -i {split_out} > {split_sort_out}")

    os.remove(split_out)
    os.system(f"rm -rf {path_subfile_in} {path_subfile_out}")

    return


def merge_peak_bed(path_in, list_input, num_process):
    pool = Pool(processes=num_process)
    func_merge = partial(pre_merge_bed, path_in)
    pool.map(func_merge, list_input)
    pool.close()
    print("Preparing merge bed is completed!")

    list_path_subout = []
    for dict_in in list_input:
        accession_ids = dict_in['accession_ids']
        if len(accession_ids) == 1:
            continue
        flank_percent = dict_in['flank_percent']
        term_name = dict_in['term_name']
        sub_path_out = dict_in['path']
        merge_out = os.path.join(sub_path_out, f"{term_name}.bed.merge")

        path_subfile_in = \
            os.path.join(sub_path_out, f'{term_name}_subfiles_in')
        if not os.path.exists(path_subfile_in):
            os.mkdir(path_subfile_in)
        path_subfile_out = \
            os.path.join(sub_path_out, f'{term_name}_subfiles_out')
        if not os.path.exists(path_subfile_out):
            os.mkdir(path_subfile_out)
        dict_in['path_subfile_in'] = path_subfile_in
        dict_in['path_subfile_out'] = path_subfile_out
        list_path_subout.append(dict_in)
        os.system(f"split -a 4 -d -l 10000 {merge_out} "
                  f"{path_subfile_in}/subfile")
        os.remove(merge_out)

        subfiles_in = glob.glob(path_subfile_in + '/*')
        pool = Pool(processes=num_process)
        func_split = partial(split_merge_bed, path_subfile_out, flank_percent)
        pool.map(func_split, subfiles_in)
        pool.close()
    print("Spliting merge bed is completed!")

    pool = Pool(processes=num_process)
    pool.map(concat_subfiles, list_path_subout)
    pool.close()
    print("Concating is completed!")

    return


def overlap_matrix(path_in, dict_in):
    term_name = dict_in['term_name']
    path_out = dict_in['path']
    accession_ids = dict_in['accession_ids']
    df_meta = pd.read_csv(
        os.path.join(path_in, 'metadata.simple.tsv'), sep='\t')

    file_merge = os.path.join(path_out, term_name + '.bed')
    list_bed = \
        (pd.read_csv(file_merge, sep='\t', header=None)).to_dict('records')
    term = term_name.replace('_', ' ').replace('+', '/').replace("--", "'")
    list_term_access = df_meta.loc[df_meta['Biosample term name'] == term,
                                   'File accession'].tolist()
    if accession_ids[0][:5] == 'ENCSR':
        list_score = []
        for sub_dict in list_bed:
            score_dict = dict()
            score_dict['peak_id'] = \
                f"{sub_dict[0]}:{str(sub_dict[1])}-{str(sub_dict[2])}"
            list_access = sub_dict[3].strip().split('|')
            if len(list_access) == 1:
                score_dict[list_access[0]] = sub_dict[6]
            else:
                set_access = set(list_access)
                list_row_score = sub_dict[6].strip().split('|')
                for access in list_term_access:
                    if access in set_access:
                        idx_access = list_access.index(access)
                        score_dict[access] = list_row_score[idx_access]
                    else:
                        score_dict[access] = np.nan
            list_score.append(score_dict)
        df_score = pd.DataFrame(
            list_score, columns=['peak_id'] + list_term_access
        )
        df_score.index = df_score['peak_id']
        df_score = df_score.drop('peak_id', 1)
        # write score matrix to txt file
        file_score = os.path.join(path_out, 'score.txt')
        # df_ref = df_ref.drop_duplicates()
        df_score.to_csv(file_score, sep='\t', na_rep='NA')

    if len(accession_ids) <= 1:
        return
    else:
        list_dict = []
        for sub_dict in list_bed:
            out_dict = dict()
            out_dict['peak_id'] = \
                f"{sub_dict[0]}:{str(sub_dict[1])}-{str(sub_dict[2])}"
            score_dict = dict()
            score_dict['peak_id'] = \
                f"{sub_dict[0]}:{str(sub_dict[1])}-{str(sub_dict[2])}"
            list_access = sub_dict[3].strip().split('|')
            set_access = set(list_access)
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
        # df_ref = df_ref.drop_duplicates()
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
            term_name = sub_dict['Biosample term name']
            life_stage = sub_dict['Biosample life stage']
            dict_meta.append(
                {'Biosample term name': term_name,
                 'Biosample life stage': life_stage,
                 'Biosample organ': organ,
                 'Biosample life_organ': life_stage + '_' + organ})
    meta_out = pd.DataFrame(dict_meta)
    meta_out.to_csv(os.path.join(path_out, 'meta.reference.tsv'),
                    sep='\t', index=None)

    list_input = []
    life_organs = set(meta_out['Biosample life_organ'].tolist())
    for life_organ in life_organs:
        organ = pd.unique(
            meta_out.loc[
                meta_out['Biosample life_organ'] == life_organ,
                'Biosample organ'])[0]
        life_stage = pd.unique(
            meta_out.loc[
                meta_out['Biosample life_organ'] == life_organ,
                'Biosample life stage'])[0]
        organ_meta = df_meta.loc[
                     (df_meta['Biosample organ'].apply(
                         lambda x: organ in x.strip().split(','))) &
                     (df_meta['Biosample life stage'] == life_stage), :]
        organ_path = \
            os.path.join(path_out, life_organ.replace(' ', '_'))
        if not os.path.exists(organ_path):
            os.makedirs(organ_path)
        terms = set(organ_meta['Biosample term name'].tolist())
        for term in terms:
            term_meta = \
                organ_meta.loc[organ_meta['Biosample term name'] == term, :]
            accession_ids = \
                list(set(term_meta['Experiment accession'].tolist()))
            path_term = \
                os.path.join(organ_path, term.replace(
                    ' ', '_').replace('/', '+').replace("'", '--'))
            if not os.path.exists(path_term):
                os.makedirs(path_term)
            list_input.append(
                dict(path=path_term,
                     term_name=term.replace(' ', '_').replace(
                         '/', '+').replace("'", '--'),
                     accession_ids=accession_ids,
                     flank_percent=flank_percent))

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
        life_organ = dict_in['Biosample life_organ']
        term = dict_in['Biosample term name']
        term_name = term.replace(' ', '_').replace('/', '+').replace("'", "--")
        file_label = f"{life_organ}|{term}"
        file = f"{life_organ.replace(' ', '_')}/{term_name}/{term_name}.bed"
        folder3 = os.path.join(
            path_out, f"{life_organ.replace(' ', '_')}/{term_name}"
        )
        file_in = os.path.join(path_in, file)
        file_out_unsort = os.path.join(path_out, file + '.unsort')
        file_out = os.path.join(path_out, file)
        # make folder
        if not os.path.exists(folder3):
            os.mkdir(folder3)
        file_score = os.path.join(
            path_in, f"{life_organ.replace(' ', '_')}/{term_name}/score.txt")
        file_quantile = os.path.join(folder3, 'quantile.txt')
        os.system(
            f"Rscript {os.path.join(root_path, 'correct_DHS_score.R')} "
            f"{file_score} {file_quantile}")
        df_quantile = pd.read_csv(
            file_quantile, sep='\t', header=None, index_col=0,
            names=['dhs_id', 'quantile'])
    else:
        file_in = os.path.join(path_in, dict_in['File accession'] + '.bed')
        file_out_unsort = \
            os.path.join(path_out, dict_in['File accession'] + '.unsort.bed')
        file_out = os.path.join(path_out, dict_in['File accession'] + '.bed')

    with open(file_in, 'r') as r_f:
        with open(file_out_unsort, 'w') as w_f:
            fmt_dhs = "{chrom}\t{start}\t{end}\t{label}\t{score}\t.\t" \
                      "{file_label}\t{accessions}\n"
            fmt_histone = "{chrom}\t{start}\t{end}\t{label}\t" \
                          "{score}\t{pvalue}\t{accessions}\n"
            for line in r_f:
                list_line = line.strip().split('\t')
                chrom = list_line[0]
                start = list_line[1]
                end = list_line[2]
                accessions = list_line[3]
                label = f"{type_bed}<-{chrom}:{start}-{end}"
                if type_bed == 'DHS':
                    dhs_id = f"{chrom}:{start}-{end}"
                    score = df_quantile.loc[dhs_id, 'quantile']
                    w_f.write(fmt_dhs.format(**locals()))
                else:
                    score = list_line[6]
                    pvalue = list_line[7]
                    w_f.write(fmt_histone.format(**locals()))

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
                       'Biosample organ']]
        meta_out = meta_out.drop_duplicates()
        dict_meta = []
        for sub_dict in meta_out.to_dict('records'):
            organs = sub_dict['Biosample organ'].split(',')
            for organ in organs:
                term_name = sub_dict['Biosample term name']
                life_stage = sub_dict['Biosample life stage']
                dict_meta.append(
                    {'Biosample term name': term_name,
                     'Biosample life stage': life_stage,
                     'Biosample organ': organ,
                     'Biosample life_organ': life_stage + '_' + organ})
        meta_out = pd.DataFrame(dict_meta)
        meta_out.to_csv(os.path.join(path_out, 'meta.reference.tsv'),
                        sep='\t', index=None)

    if type_bed == 'DHS':
        df_meta = pd.read_csv(
            os.path.join(path_in, 'meta.reference.tsv'), sep='\t'
        )
        list_input = df_meta.to_dict('records')
        for sub_dict in list_input:
            life_organ = sub_dict['Biosample life_organ']
            folder1 = os.path.join(path_out, f"{life_organ.replace(' ', '_')}")
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


def pre_merge_stan_bed(path_bed, dict_in):
    term_name = dict_in['term_name']
    path_out = dict_in['path']
    accession_ids = dict_in['accession_ids']

    # merge bed files from same term
    col_collapse = '2,3,5,7,8'
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
    # os.system(f"sort -k 1,1 -k2,2n {cat_out} > {sort_out}")
    os.system(f"bedtools merge -i {sort_out} "
              f"-c {col_collapse} -o {str_collapse} > {merge_out}")

    os.remove(cat_out)
    os.remove(sort_out)

    return


def split_merge_stan_bed(sub_path_out, flank_percent, subfile):
    # split merge file
    _, subfile_name = os.path.split(subfile)
    sub_split_out = os.path.join(sub_path_out, f"{subfile_name}.bed.unsort")
    sub_split_sort_out = \
        os.path.join(sub_path_out, f"{subfile_name}.bed.split.sort")
    sub_final_out = os.path.join(sub_path_out, f"{subfile_name}.bed")
    with open(subfile, 'r') as r_f:
        with open(sub_split_out, 'w') as w_f:
            fmt = "{chrom}\t{start}\t{end}\t{dhs_id}\t{score}\t.\t" \
                  "{label}\t{accessions}\n"
            for line in r_f:
                list_line = line.strip().split('\t')
                chrom = list_line[0]
                array_start = np.array(
                    [int(num) for num in list_line[3].strip().split(',')])
                array_end = np.array(
                    [int(num) for num in list_line[4].strip().split(',')])
                array_label = np.array(
                    [num for num in list_line[6].strip().split(',')])
                array_score = np.array(
                    [float(num) for num in list_line[5].strip().split(',')])
                array_accessions = np.array(
                    [num for num in list_line[7].strip().split(',')])
                array_length = array_end - array_start
                list_dict = []
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
                    select_p_value = array_score[select_bool]
                    select_accessions = array_accessions[select_bool]
                    # reset arrays
                    array_start = array_start[~select_bool]
                    array_end = array_end[~select_bool]
                    array_label = array_label[~select_bool]
                    array_score = array_score[~select_bool]
                    array_accessions = array_accessions[~select_bool]
                    array_length = array_end - array_start
                    # write new rows
                    start = str(np.min(select_start))
                    end = str(np.max(select_end))
                    dhs_id = f"DHS<-{chrom}:{start}-{end}"
                    label = '/'.join(map(str, select_label.tolist()))
                    score = str(np.max(select_p_value.tolist()))
                    accessions = '|'.join(map(str, select_accessions.tolist()))
                    w_f.write(fmt.format(**locals()))

    os.system(f"bedtools sort -i {sub_split_out} > {sub_split_sort_out}")

    with open(sub_final_out, 'w') as w_final:
        old_end = 0
        old_chrom = '0'
        with open(sub_split_sort_out, 'r') as r_sort:
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
                    list_line[3] = \
                        f"DHS<-{chrom}:{list_line[1]}-{list_line[2]}"
                    w_final.write('\t'.join(list_line) + '\n')
                old_end = int(list_line[2])
                old_chrom = chrom

    df_final = pd.read_csv(sub_final_out, sep='\t', header=None,
                           names=['chrom', 'start', 'end', 'dhs_id', 'score',
                                  'strand', 'label', 'accessions'])
    len_df = df_final.shape[0]
    region_length = df_final['end'] - df_final['start']
    df_narrow = df_final.loc[region_length <= 120, :]
    set_index = set(df_narrow.index)
    for i in df_narrow.index:
        if (i - 1) in set_index:
            continue
        if i + 1 == len_df:
            if df_final.loc[i, 'start'] == df_final.loc[i-1, 'end'] + 1:
                df_final.loc[i-1, 'end'] = df_final.loc[i, 'end']
                df_final.loc[i-1, 'dhs_id'] = \
                    f"DHS<-{df_final.loc[i, 'chrom']}:" \
                    f"{df_final.loc[i-1, 'start']}-{df_final.loc[i-1, 'end']}"
                df_final.loc[i-1, 'score'] = np.max(df_final.loc[i-1:i+1, 'score'])
                df_final.loc[i-1, 'label'] = \
                    '/'.join((df_final.loc[i-1:i+1, 'label']).tolist())
                df_final.loc[i-1, 'accessions'] = \
                    '|'.join((df_final.loc[i-1:i+1, 'accessions']).tolist())
                df_final = df_final.drop(i)
            continue
        # try:
        #     a = (df_final.loc[i, 'end'] == df_final.loc[i+1, 'start'] - 1)
        # except KeyError:
        #     print(sub_final_out)
        #     print(i + 1 == df_final.shape[0])
        #     print(i + 1)
        #     print(df_final.shape[0])
        if i == 0:
            try:
                a = (df_final.loc[i, 'end'] == df_final.loc[i+1, 'start'] - 1)
            except KeyError:
                print(sub_final_out)
                print(i)
                print(i + 1)
                # print(df_final.shape[0])
            if df_final.loc[i, 'end'] == df_final.loc[i+1, 'start'] - 1:
                df_final.loc[i+1, 'start'] = df_final.loc[i, 'start']
                df_final.loc[i+1, 'dhs_id'] = \
                    f"DHS<-{df_final.loc[i, 'chrom']}:" \
                    f"{df_final.loc[i+1, 'start']}-{df_final.loc[i+1, 'end']}"
                df_final.loc[i+1, 'score'] = np.max(df_final.loc[i:i+2, 'score'])
                df_final.loc[i+1, 'label'] = \
                    '/'.join((df_final.loc[i:i+2, 'label']).tolist())
                df_final.loc[i+1, 'accessions'] = \
                    '|'.join((df_final.loc[i:i+2, 'accessions']).tolist())
                df_final = df_final.drop(i)
            continue
        if df_final.loc[i, 'end'] == df_final.loc[i+1, 'start'] - 1:
            df_final.loc[i+1, 'start'] = df_final.loc[i, 'start']
            df_final.loc[i+1, 'dhs_id'] = \
                f"DHS<-{df_final.loc[i, 'chrom']}:" \
                f"{df_final.loc[i+1, 'start']}-{df_final.loc[i+1, 'end']}"
            df_final.loc[i+1, 'score'] = np.max(df_final.loc[i:i+2, 'score'])
            df_final.loc[i+1, 'label'] = \
                '/'.join((df_final.loc[i:i+2, 'label']).tolist())
            df_final.loc[i+1, 'accessions'] = \
                '|'.join((df_final.loc[i:i+2, 'accessions']).tolist())
            df_final = df_final.drop(i)
            continue
        elif df_final.loc[i, 'start'] == df_final.loc[i-1, 'end'] + 1:
            df_final.loc[i-1, 'end'] = df_final.loc[i, 'end']
            df_final.loc[i-1, 'dhs_id'] = \
                f"DHS<-{df_final.loc[i, 'chrom']}:" \
                f"{df_final.loc[i-1, 'start']}-{df_final.loc[i-1, 'end']}"
            df_final.loc[i-1, 'score'] = np.max(df_final.loc[i-1:i+1, 'score'])
            df_final.loc[i-1, 'label'] = \
                '/'.join((df_final.loc[i-1:i+1, 'label']).tolist())
            df_final.loc[i-1, 'accessions'] = \
                '|'.join((df_final.loc[i-1:i+1, 'accessions']).tolist())
            df_final = df_final.drop(i)

    df_final.to_csv(sub_final_out, sep='\t', index=None, header=None)

    return


def merge_standard_bed(path_in, list_input, num_process):
    pool = Pool(processes=num_process)
    func_merge = partial(pre_merge_stan_bed, path_in)
    pool.map(func_merge, list_input)
    pool.close()
    print("Preparing merge bed is completed!")

    list_path_subout = []
    for dict_in in list_input:
        accession_ids = dict_in['accession_ids']
        if len(accession_ids) == 1:
            continue
        flank_percent = dict_in['flank_percent']
        term_name = dict_in['term_name']
        sub_path_out = dict_in['path']
        merge_out = os.path.join(sub_path_out, f"{term_name}.bed.merge")

        path_subfile_in = \
            os.path.join(sub_path_out, f'{term_name}_subfiles_in')
        if not os.path.exists(path_subfile_in):
            os.mkdir(path_subfile_in)
        path_subfile_out = \
            os.path.join(sub_path_out, f'{term_name}_subfiles_out')
        if not os.path.exists(path_subfile_out):
            os.mkdir(path_subfile_out)
        dict_in['path_subfile_in'] = path_subfile_in
        dict_in['path_subfile_out'] = path_subfile_out
        list_path_subout.append(dict_in)
        os.system(f"split -a 4 -d -l 10000 {merge_out} "
                  f"{path_subfile_in}/subfile")
        os.remove(merge_out)

        subfiles_in = glob.glob(path_subfile_in + '/*')
        pool = Pool(processes=num_process)
        func_split = partial(
            split_merge_stan_bed, path_subfile_out, flank_percent)
        pool.map(func_split, subfiles_in)
        pool.close()
    print("Spliting merge bed is completed!")

    pool = Pool(processes=num_process)
    pool.map(concat_subfiles, list_path_subout)
    pool.close()
    print("Concating is completed!")

    return


def label_mat(labels, dict_in):
    dict_label = dict()
    dict_label['peak_id'] = \
        f"{dict_in[0]}:{str(dict_in[1])}-{str(dict_in[2])}"
    set_label = set(dict_in[6].strip().split('/'))
    for label in labels:
        if label in set_label:
            dict_label[label] = 1
        else:
            dict_label[label] = 0

    return dict_label


def accession_mat(accessions, dict_in):
    dict_accession = dict()
    dict_accession['peak_id'] = \
        f"{dict_in[0]}:{str(dict_in[1])}-{str(dict_in[2])}"
    set_accession = set(dict_in[7].strip().split('|'))
    for accession in accessions:
        if accession in set_accession:
            dict_accession[accession] = 1
        else:
            dict_accession[accession] = 0

    return dict_accession


def sub_merge(dict_in):
    sub_ref_meta = dict_in['sub_ref_meta']
    sub_meta = dict_in['sub_meta']
    sub_path_in = dict_in['path_in']
    sub_path_out = dict_in['path_out']
    life_organ = dict_in['life_organ']
    bool_accession = dict_in['bool_accession']
    bool_plot = dict_in['bool_plot']
    num_process = dict_in['num_process']
    file_out = \
        os.path.join(sub_path_out, life_organ.replace(' ', '_') + '.bed')
    if os.path.exists(sub_path_out):
        os.system(f"rm -rf {sub_path_out}")
    os.mkdir(sub_path_out)
    list_meta = sub_ref_meta.to_dict("records")
    accession_ids = [
        '/'.join([sub_dict['Biosample term name'].replace(
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
            fmt = "{chrom}\t{start}\t{end}\t{dhs_id}\t{score}\t.\t{label}\t" \
                  "{accessions}\n"
            with open(os.path.join(sub_path_in, accession_ids[0] + '.bed'),
                      'r') as r_f:
                for line in r_f:
                    list_line = line.strip().split('\t')
                    dict_out = dict(
                        chrom=list_line[0], start=list_line[1],
                        end=list_line[2], label=list_line[6],
                        dhs_id=f"DHS<-{list_line[0]}:"
                               f"{list_line[1]}-{list_line[2]}",
                        score=list_line[4], accessions=list_line[7]
                    )
                    w_f.write(fmt.format(**dict_out))
        return
    elif sub_ref_meta.shape[0] > 1:
        dict_merge = [dict(
            path=sub_path_out,
            term_name=life_organ.replace(' ', '_'),
            accession_ids=accession_ids,
            flank_percent=0.5)]
        merge_standard_bed(sub_path_in, dict_merge, num_process)
        labels = [f"{sub_dict['Biosample organ']}|"
                  f"{sub_dict['Biosample life stage']}|"
                  f"{sub_dict['Biosample term name']}"
                  for sub_dict in list_meta]
        list_bed = \
            (pd.read_csv(file_out, sep='\t', header=None)).to_dict('records')

        # label matrix
        func_label = partial(label_mat, labels)
        list_label = [func_label(sub_dict) for sub_dict in list_bed]
        df_label = pd.DataFrame(list_label, columns=['peak_id'] + labels)
        df_label.index = df_label['peak_id']
        df_label = df_label.drop('peak_id', 1)
        # write score matrix to txt file
        mat_label = os.path.join(sub_path_out, 'label_matrix.txt')
        df_label.to_csv(mat_label, sep='\t')

        # Jaccard distance
        list_out = []
        list_com = combinations(df_label.columns, 2)
        for com in list_com:
            jdist = pdist(
                (df_label.loc[
                    (df_label.loc[:, com[0]] != 0) |
                    (df_label.loc[:, com[1]] != 0), com]).T,
                'jaccard')[0]
            list_out.append({'Name': life_organ, 'Combination': com,
                             'Jaccard distance': jdist})

        df_out = pd.DataFrame(list_out)

        if bool_plot:
            # scatter plot
            os.system(
                f"Rscript scatter.plot.organ.R {mat_label} {sub_path_out}")

        # accession matrix
        if bool_accession:
            accessions = [sub_dict['File accession']
                          for sub_dict in sub_meta.to_dict("records")]
            func_accession = partial(accession_mat, accessions)
            list_accession = \
                [func_accession(sub_dict) for sub_dict in list_bed]
            df_accession = \
                pd.DataFrame(list_accession, columns=['peak_id'] + accessions)
            df_accession.index = df_accession['peak_id']
            df_accession = df_accession.drop('peak_id', 1)
            # write score matrix to txt file
            mat_accession = os.path.join(sub_path_out, 'accession_matrix.txt')
            df_accession.to_csv(mat_accession, sep='\t')

            if bool_plot:
                # scatter plot
                path_in = '/'.join(sub_path_in.split('/')[:-1])
                os.system(
                    f"Rscript scatter.plot.organ.accession.R {mat_accession} "
                    f"{os.path.join(path_in, 'metadata.simple.tsv')} "
                    f"{sub_path_out}")

        return df_out


def organ_mat(organs, dict_in):
    organ_dict = dict()
    organ_dict['peak_id'] = \
        f"{dict_in[0]}:{str(dict_in[1])}-{str(dict_in[2])}"
    list_organ = set([val.split('|')[0]
                      for val in dict_in[6].strip().split('/')])
    for organ in organs:
        if organ in list_organ:
            organ_dict[organ] = 1
        else:
            organ_dict[organ] = 0

    return organ_dict


def calculate_overlap(df_organ, com):
    # Overlap ratio
    # len_list1 = (df_organ.loc[
    #              df_organ.loc[:, com[0]] != 0, :]).shape[0]
    # len_list2 = (df_organ.loc[
    #              df_organ.loc[:, com[1]] != 0, :]).shape[0]
    # total = min(len_list1, len_list2)
    # len_overlap = (df_organ.loc[
    #                (df_organ.loc[:, com[0]] != 0) &
    #                (df_organ.loc[:, com[1]] != 0), :]).shape[0]

    # Jaccard distance
    jdist = pdist(
        (df_organ.loc[
            (df_organ.loc[:, com[0]] != 0) |
            (df_organ.loc[:, com[1]] != 0), com]).T,
        'jaccard')[0]

    return {'Combination': com, 'Jaccard distance': jdist}


def merge_organ_cluster(path_in, path_out, num_process,
                        bool_accession=True, bool_plot=True):
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

    life_organs = list(set(df_meta_ref['Biosample life_organ'].tolist()))

    list_df = []
    for life_organ in life_organs:
        sub_ref_meta = df_meta_ref.loc[
                   df_meta_ref['Biosample life_organ'] == life_organ, :]
        organ = pd.unique(sub_ref_meta['Biosample organ'])[0]
        life_stage = pd.unique(sub_ref_meta['Biosample life stage'])[0]
        sub_meta = df_meta.loc[
                   (df_meta['Biosample organ'].apply(
                       lambda x: organ in x.strip().split(','))) &
                   (df_meta['Biosample life stage'] == life_stage), :]
        dict_in = \
            dict(sub_ref_meta=sub_ref_meta, life_organ=life_organ,
                 sub_meta=sub_meta, num_process=num_process,
                 path_out=os.path.join(path_out, life_organ.replace(' ', '_')),
                 path_in=os.path.join(path_in, life_organ.replace(' ', '_')),
                 bool_accession=bool_accession, bool_plot=bool_plot)
        sub_df = sub_merge(dict_in)
        list_df.append(sub_df)

    df_overlap = pd.concat(list_df, sort=False)
    df_overlap.to_csv(os.path.join(path_out, 'overlap.txt'), sep='\t')
    print("Integration of organ ---- completed!")

    # merge all organ
    organ_folders = os.listdir(path_out)
    accession_ids = []
    for organ_folder in organ_folders:
        path_organ = os.path.join(path_out, organ_folder)
        if (os.path.isdir(path_organ)) & (organ_folder != 'flank'):
            accession_ids.append(os.path.join(path_organ, organ_folder))

    dict_merge = [dict(
        path=path_out,
        term_name='all_organs',
        accession_ids=accession_ids,
        flank_percent=0.5)]
    merge_standard_bed(path_out, dict_merge, num_process=20)
    list_bed = \
        (pd.read_csv(os.path.join(path_out, 'all_organs.bed'),
                     sep='\t', header=None)).to_dict('records')
    print("Integration of all data ---- completed!")

    # label matrix
    labels = [f"{sub_dict['Biosample life stage']}_"
              f"{sub_dict['Biosample organ']}|"
              f"{sub_dict['Biosample term name']}"
              for sub_dict in df_meta_ref.to_dict("records")]
    pool = Pool(processes=num_process)
    func_label = partial(label_mat, labels)
    list_label = pool.map(func_label, list_bed)
    pool.close()

    df_label = pd.DataFrame(list_label, columns=['peak_id'] + labels)
    df_label.index = df_label['peak_id']
    df_label = df_label.drop('peak_id', 1)
    # write score matrix to txt file
    mat_label = os.path.join(path_out, 'label_matrix.txt')
    df_label.to_csv(mat_label, sep='\t')
    print("Label matrix ---- completed!")

    # experiment matrix
    accessions = [sub_dict['File accession']
                  for sub_dict in df_meta.to_dict("records")]
    pool = Pool(processes=num_process)
    func_accession = partial(accession_mat, accessions)
    list_accession = pool.map(func_accession, list_bed)
    pool.close()

    df_accession = \
        pd.DataFrame(list_accession, columns=['peak_id'] + accessions)
    df_accession.index = df_accession['peak_id']
    df_accession = df_accession.drop('peak_id', 1)
    # write score matrix to txt file
    mat_accession = os.path.join(path_out, 'accession_matrix.txt')
    df_accession.to_csv(mat_accession, sep='\t')
    print("Experiment matrix ---- completed!")

    # organ matrix
    pool = Pool(processes=num_process)
    func_organ = partial(organ_mat, life_organs)
    list_organ = pool.map(func_organ, list_bed)
    pool.close()

    df_organ = pd.DataFrame(list_organ, columns=['peak_id'] + life_organs)
    df_organ.index = df_organ['peak_id']
    df_organ = df_organ.drop('peak_id', 1)
    # write score matrix to txt file
    mat_organ = os.path.join(path_out, 'organ_matrix.txt')
    df_organ.to_csv(mat_organ, sep='\t')
    print("Organ matrix ---- completed!")

    # calculate overlap ratio(有问题)
    # list_com = combinations(df_organ.columns, 2)
    # pool = Pool(processes=num_process)
    # func_overlap = partial(calculate_overlap, df_organ)
    # list_out = pool.map(func_overlap, list_com)
    # pool.close()
    #
    # df_overlap = pd.DataFrame(list_out)
    # df_overlap.to_csv(
    #     os.path.join(path_out, 'all_organs_overlap.txt'),
    #     sep='\t', index=None
    # )
    # print("Overlap ratio ---- completed!")

    # scatter plot
    if bool_plot:
        os.system(f"Rscript scatter.plot.all.R {mat_label} "
                  f"{os.path.join(path_in, 'metadata.simple.tsv')} {path_out}")

        if bool_accession:
            os.system(f"Rscript scatter.plot.all.accession.R {mat_accession} "
                      f"{os.path.join(path_in, 'metadata.simple.tsv')} "
                      f"{path_out}")

    return


def merge_suborgan(path_in, path_out, meta_suborgan, num_process):
    df_meta_ref = pd.read_csv(meta_suborgan, sep='\t')
    # combine life and organ
    df_meta_ref['Biosample life_organ'] = df_meta_ref.apply(
        lambda x: x['Biosample life stage'] + '_' + x['Biosample organ'],
        axis=1
    )
    df_ref_filter = df_meta_ref.dropna()

    suborgans = list(
        set([organ for organ in df_ref_filter['Biosample suborgan'].tolist()])
    )
    life_organs = list(set(df_ref_filter['Biosample life_organ'].tolist()))

    list_df = []
    for life_organ in life_organs:
        for suborgan in suborgans:
            sub_ref_meta = df_meta_ref.loc[
                (df_meta_ref['Biosample suborgan'] == suborgan) &
                (df_meta_ref['Biosample life_organ'] == life_organ), :]
            if sub_ref_meta.shape[0] == 0:
                continue
            path_organ = os.path.join(path_out, life_organ.replace(' ', '_'))
            path_sub_organ = os.path.join(path_organ,
                                          suborgan.replace(' ', '_'))
            dict_in = \
                dict(sub_ref_meta=sub_ref_meta, life_organ=suborgan,
                     sub_meta=sub_ref_meta, num_process=num_process,
                     path_out=path_sub_organ,
                     path_in=os.path.join(path_in,
                                          life_organ.replace(' ', '_')),
                     bool_accession=False, bool_plot=False)
            sub_df = sub_merge(dict_in)
            list_df.append(sub_df)

    df_overlap = pd.concat(list_df, sort=False)
    df_overlap.to_csv(os.path.join(path_out, 'overlap.suborgan.txt'), sep='\t')
    print("Integration of suborgan ---- completed!")

    return


if __name__ == '__main__':
    time_start = time()
    # parameters
    num_cpu = 40
    # path_root = '/local/zy/PEI'
    path_root = '/lustre/tianlab/zhangyu/PEI'
    path_origin = path_root + '/origin_data'
    path_mid = path_root + '/mid_data_correct'
    # get bed file annotating protein-coding genes
    gtf_file_hg19 = path_origin + '/ENCODE/gencode.v19.annotation.gtf'
    protein_file_hg19 = path_origin + '/gene/genes.protein.gencode.v19.bed'
    promoter_file_hg19 = \
        path_origin + '/gene/promoters.up2k.protein.gencode.v19.bed'
    promoter_file_hg19_merge = \
        path_origin + '/gene/promoters.up2k.protein.gencode.v19.merge.bed'
    exon_file_hg19 = path_origin + '/gene/exon.protein.gencode.v19.bed'
    # generate_gene_file(gtf_file_hg19, protein_file_hg19, promoter_file_hg19,
    #                    promoter_file_hg19_merge, exon_file_hg19)

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

    # read reference organ
    ref_organ = path_origin + '/ENCODE/metadata/organ_ref.txt'
    with open(ref_organ, 'r') as r_ref:
        set_organs = set([organ.strip() for organ in r_ref])
    # organ complement
    file_complement = \
        path_origin + '/ENCODE/metadata/complement_organ.txt'
    df_complement = pd.read_csv(file_complement, sep='\t')
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
    df_meta_dhs = modify_meta(df_meta_dhs, set_organs, df_complement)
    meta_dhs = path_mid + '/tissue/ENCODE/DNase-seq/metadata.simple.tsv'
    df_meta_dhs.to_csv(meta_dhs, sep='\t', index=None)
    print("DHS metadata ---- completed")

    # hg38 to hg19
    path_hg38tohg19 = path_mid + '/tissue/ENCODE/DNase-seq/GRCh38tohg19'
    hg38tohg19(path_dhs, path_hg38tohg19, meta_dhs, num_cpu)
    print("Format conversion: hg38 -> hg19 ---- completed")

    # integrate files from same experiment
    path_exp_dhs = \
        path_mid + '/tissue/ENCODE/DNase-seq/GRCh38tohg19_experiment'
    merge_experiment(path_hg38tohg19, path_exp_dhs, 0.5, num_cpu)
    print("Integration of files from same experiment ---- completed")

    # build DHS reference
    path_dhs_hg38tohg19 = path_mid + '/tissue/DHS/GRCh38tohg19/'
    unique_bed_files(path_exp_dhs, path_dhs_hg38tohg19, 0.5, num_cpu)
    print("Integration of files from same term ---- completed")

    # standardization
    path_dhs_stan = path_mid + '/tissue/DHS/GRCh38tohg19_standard'
    standardize_bed(path_dhs_hg38tohg19, path_dhs_stan, 'DHS', num_cpu)
    print('Standardization of DHS completed!')

    # merge and cluster
    path_dhs_cluster = path_mid + '/tissue/DHS/GRCh38tohg19_cluster'
    merge_organ_cluster(path_dhs_stan, path_dhs_cluster, num_cpu, False, False)
    print('Cluster and merge of DHS completed!')

    # merge sub-organ
    meta_suborgan_dhs = \
        path_origin + '/meta_file/meta.reference.tsv'
    merge_suborgan(path_dhs_stan, path_dhs_cluster, meta_suborgan_dhs, num_cpu)

    # preparation of bed files of histone and TF
    # H3K4me3
    path_h3k4me3 = \
        path_origin + '/ENCODE/histone_ChIP-seq/H3K4me3'
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
    df_meta_h3k4me3 = modify_meta(df_meta_h3k4me3, set_organs, df_complement)
    meta_h3k4me3 = path_mid + '/tissue/ENCODE/histone_ChIP-seq/' \
                              'metadata.simple.H3K4me3.tsv'
    df_meta_h3k4me3.to_csv(meta_h3k4me3, sep='\t', index=None)
    print("H3K4me3 metadata ---- completed")

    # hg38 to hg19
    path_hg38tohg19 = \
        path_mid + '/tissue/ENCODE/histone_ChIP-seq/H3K4me3'
    hg38tohg19(path_h3k4me3, path_hg38tohg19, meta_h3k4me3, num_cpu)
    print("Format conversion: hg38 -> hg19 ---- completed")

    # standardization
    path_h3k4me3_stan = \
        path_mid + '/tissue/ENCODE/histone_ChIP-seq/H3K4me3_standard'
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
    df_meta_h3k27ac = modify_meta(df_meta_h3k27ac, set_organs, df_complement)
    meta_h3k27ac = path_mid + '/tissue/ENCODE/histone_ChIP-seq/' \
                              'metadata.simple.H3K27ac.tsv'
    df_meta_h3k27ac.to_csv(meta_h3k27ac, sep='\t', index=None)
    print("H3K27ac metadata ---- completed")

    # hg38 to hg19
    path_hg38tohg19 = path_mid + '/tissue/ENCODE/histone_ChIP-seq/H3K27ac'
    hg38tohg19(path_h3k27ac, path_hg38tohg19, meta_h3k27ac, num_cpu)
    print("Format conversion: hg38 -> hg19 ---- completed")

    # standardization
    path_h3k27ac_stan = \
        path_mid + '/tissue/ENCODE/histone_ChIP-seq/H3K27ac_standard'
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
    df_meta_ctcf = modify_meta(df_meta_ctcf, set_organs, df_complement)
    meta_ctcf = path_mid + '/tissue/ENCODE/' \
                           'TF_ChIP-seq/metadata.simple.CTCF.tsv'
    df_meta_ctcf.to_csv(meta_ctcf, sep='\t', index=None)
    print("CTCF metadata ---- completed")

    # hg38 to hg19
    path_hg38tohg19 = path_mid + '/tissue/ENCODE/TF_ChIP-seq/CTCF'
    hg38tohg19(path_ctcf, path_hg38tohg19, meta_ctcf, num_cpu)
    print("Format conversion: hg38 -> hg19 ---- completed")

    # standardization
    path_ctcf_stan = path_mid + '/tissue/ENCODE/TF_ChIP-seq/CTCF_standard'
    standardize_bed(path_hg38tohg19, path_ctcf_stan, 'CTCF', num_cpu)
    print('Standardization of CTCF completed!')

    time_end = time()
    print(time_end - time_start)
