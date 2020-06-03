#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: unify_interactions.py
# @time: 2020/4/20 10:19

from time import time
import os
import pandas as pd
import numpy as np
from multiprocessing import Pool
file_chain = '/lustre/tianlab/tools/files_liftOver/hg38ToHg19.over.chain.gz'
liftover = '/lustre/tianlab/tools/liftOver'


def encode_chiapet(file_in, file_out, cutoff=2):
    fmt = "{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t" \
          "{loop_id}\t{score}\n"
    with open(file_in, 'r') as r_in:
        with open(file_out, 'w') as w_out:
            for i, line in enumerate(r_in):
                list_line = line.strip().split('\t')
                pair = list_line[3].split(',')
                score = pair[-1]
                if float(score) < cutoff:
                    continue
                regions = pair[0].split('-')
                region1 = regions[0]
                chrom1 = region1.split(':')[0]
                region2 = regions[1]
                chrom2 = region2.split(':')[0]
                if chrom1 != chrom2:
                    continue
                start1 = region1.split(':')[1].split('..')[0]
                end1 = region1.split(':')[1].split('..')[1]
                start2 = region2.split(':')[1].split('..')[0]
                end2 = region2.split(':')[1].split('..')[1]
                loop_id = str(i)
                w_out.write(fmt.format(**locals()))

    return


def transform_insitu(file_in, file_out):
    with open(file_out, 'w') as w_out:
        fmt = "{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t" \
              "{loop_id}\t{score}\n"
        with open(file_in, 'r') as r_in:
            for i, line in enumerate(r_in):
                list_line = line.strip().split('\t')
                if list_line[0] == 'chr1':
                    continue
                # score = list_line[12]
                score = '.'
                loop_id = str(i)
                chrom1 = 'chr' + list_line[0]
                start1 = (list_line[1])
                end1 = list_line[2]
                chrom2 = 'chr' + list_line[3]
                start2 = list_line[4]
                end2 = list_line[5]
                if chrom1 == chrom2:
                    w_out.write(fmt.format(**locals()))

    return


def tang_cell_2015(file_in, file_out, cutoff=2):
    with open(file_out, 'w') as w_out:
        fmt = "{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t" \
              "{loop_id}\t{score}\n"
        with open(file_in, 'r') as r_in:
            for i, line in enumerate(r_in):
                list_line = line.strip().split('\t')
                score = list_line[6]
                if float(score) < cutoff:
                    continue
                loop_id = str(i)
                chrom1 = list_line[0]
                chrom2 = list_line[3]
                if chrom1 != chrom2:
                    continue
                start1 = list_line[1]
                end1 = list_line[2]
                start2 = list_line[4]
                end2 = list_line[5]
                w_out.write(fmt.format(**locals()))

    return


def heidari_genomeres_2014(file_in, file_out):
    with open(file_out, 'w') as w_out:
        fmt = "{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t" \
              "{loop_id}\t{score}\n"
        with open(file_in, 'r') as r_in:
            for i, line in enumerate(r_in):
                list_line = line.strip().split('\t')
                # score = f"{list_line[9]};{list_line[11]}"
                score = list_line[11]
                loop_id = str(i)
                chrom1 = list_line[0]
                chrom2 = list_line[3]
                if chrom1 != chrom2:
                    continue
                start1 = list_line[1]
                end1 = list_line[2]
                start2 = list_line[4]
                end2 = list_line[5]
                w_out.write(fmt.format(**locals()))

    return


def mumbach(file_in, file_out):
    with open(file_out, 'w') as w_out:
        fmt = "{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t" \
              "{loop_id}\t{score}\n"
        with open(file_in, 'r') as r_in:
            for i, line in enumerate(r_in):
                list_line = line.strip().split('\t')
                if list_line[0][:3] == 'chr':
                    continue
                score = '.'
                loop_id = str(i)
                chrom1 = 'chr' + list_line[0]
                chrom2 = 'chr' + list_line[3]
                if chrom1 != chrom2:
                    continue
                start1 = list_line[1]
                end1 = list_line[2]
                start2 = list_line[4]
                end2 = list_line[5]
                w_out.write(fmt.format(**locals()))

    return


def transform_microc(file_in, file_out):
    file_tmp = file_out + '.tmp'
    with open(file_tmp, 'w') as w_out:
        fmt = "{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t" \
              "{loop_id}\t{score}\n"
        with open(file_in, 'r') as r_in:
            for i, line in enumerate(r_in):
                list_line = line.strip().split('\t')
                if list_line[0] == 'chrom1':
                    continue
                # score = list_line[12]
                score = '.'
                loop_id = str(i)
                chrom1 = list_line[0]
                start1 = (list_line[1])
                end1 = list_line[2]
                chrom2 = list_line[3]
                start2 = list_line[4]
                end2 = list_line[5]
                if chrom1 == chrom2:
                    w_out.write(fmt.format(**locals()))

    file_region1 = file_out + '.region1'
    file_region2 = file_out + '.region2'
    file_region1_hg19 = file_out + '.hg19.region1'
    file_region2_hg19 = file_out + '.hg19.region2'
    file_region1_ummap = file_out + '.region1.umap'
    file_region2_ummap = file_out + '.region2.umap'
    os.system(f"cut -f 1,2,3,7 {file_tmp} > {file_region1}")
    os.system(f"cut -f 4,5,6,7 {file_tmp} > {file_region2}")
    os.system(f"{liftover} {file_region1} {file_chain} "
              f"{file_region1_hg19} {file_region1_ummap}")
    os.system(f"{liftover} {file_region2} {file_chain} "
              f"{file_region2_hg19} {file_region2_ummap}")

    df_region1 = pd.read_csv(file_region1_hg19, sep='\t', header=None)
    df_region2 = pd.read_csv(file_region2_hg19, sep='\t', header=None)
    df_pair = pd.merge(df_region1, df_region2, how='inner', on=3)
    df_pair[4] = df_pair[3]
    df_pair = df_pair.drop(3, axis=1)
    df_pair[5] = np.full((df_pair.shape[0], 1), np.nan)
    df_pair.to_csv(file_out, sep='\t', index=None, header=None, na_rep='.')

    os.remove(file_tmp)
    os.remove(file_region1)
    os.remove(file_region2)
    os.remove(file_region1_hg19)
    os.remove(file_region2_hg19)
    os.remove(file_region1_ummap)
    os.remove(file_region2_ummap)

    return


def uniform_file(dict_in):
    term = dict_in['Biosample term name']
    method = dict_in['Method']
    source = dict_in['Source']
    filename = dict_in['Filename']

    file_in = os.path.join(path_chrom_inter, f"{method}/{source}/{filename}")
    path_term = os.path.join(path_label, term)
    if not os.path.exists(path_term):
        try:
            os.mkdir(path_term)
        except FileExistsError:
            pass
    file_out = os.path.join(
        path_term, f"{method}__{source}__{filename[:-4]}.uniform")

    if source in {'Fullwood_Nature_2009', 'Li_Cell_2012'}:
        encode_chiapet(file_in, file_out, cutoff=cutoff_pet)
    elif source == 'Rao_Cell_2014':
        transform_insitu(file_in, file_out)
    elif source == 'Tang_Cell_2015':
        tang_cell_2015(file_in, file_out, cutoff=cutoff_pet)
    elif source == 'Heidari_GenomeRes_2014':
        heidari_genomeres_2014(file_in, file_out)
    elif source in {'Mumbach_NM_2016', 'Mumbach_NG_2017'}:
        mumbach(file_in, file_out)
    elif source == 'Krietenstein_MolecularCell_2020':
        transform_microc(file_in, file_out)

    return


if __name__ == '__main__':
    time_start = time()
    cutoff_pet = 2
    path_root = '/lustre/tianlab/zhangyu/PEI'
    path_origin = path_root + '/origin_data'
    path_mid = path_root + '/mid_data_correct'

    path_chrom_inter = path_origin + '/Chromatin_interactions/'
    path_label = path_mid + '/training_label/label_interactions'

    flie_meta = os.path.join(path_label, 'meta_label.txt')
    df_meta = pd.read_csv(flie_meta, sep='\t')

    pool = Pool(processes=40)
    pool.map(uniform_file, df_meta.to_dict('records'))
    pool.close()

    time_end = time()
    print(time_end - time_start)
