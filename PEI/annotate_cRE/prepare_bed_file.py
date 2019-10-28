#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: prepare_bed_file.py
# @time: 10/23/19 10:24 PM

from time import time
import pandas as pd
import os
from collections import defaultdict
from multiprocessing import Pool
from functools import partial


def hg19_filter(meta_in):
    df_meta = pd.read_csv(meta_in, sep='\t')
    df_meta_hg19 = df_meta.loc[
        df_meta['Assembly'] == 'hg19',
        ['File accession', 'Experiment accession', 'Biosample term id',
         'Biosample term name', 'Biosample type', 'Biosample treatments',
         'Biosample genetic modifications methods']]
    df_meta_hg19 = df_meta_hg19.loc[
        df_meta['File Status'] == 'released', :]
    df_meta_hg19.index = df_meta_hg19['File accession']

    return df_meta_hg19


def build_dict_attr(path_meta):
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
    df_access = df_meta['File accession']
    df_attr = df_access.apply(lambda x: ','.join(dict_attr[x]))
    df_attr.name = column_name
    df_out = pd.concat([df_meta, df_attr], axis=1, sort=False)

    return df_out


def merge_bed(path_bed, dict_in):
    term_name = dict_in['term_name'].replace(' ', '_').replace('/', ';')
    path_out = dict_in['path']
    cat_out = os.path.join(path_out, f"{term_name}.cat.bed")
    sort_out = os.path.join(path_out, f"{term_name}.sort.bed")
    bed_out = os.path.join(path_out, f"{term_name}.bed")
    cat_in = ' '.join([os.path.join(path_bed, acce_id + '.bed')
                       for acce_id in dict_in['accession_ids']])
    os.system(f"cat {cat_in} > {cat_out}")
    os.system(f"sortBed -i {cat_out} > {sort_out}")
    os.system(f"bedtools merge -i {sort_out} -c 7 -o sum > {bed_out}")

    return


def unique_bed_files(path_bed, metafile, path_out):
    df_meta = pd.read_csv(metafile, sep='\t')
    # new meta file
    new_meta = df_meta.loc[:, ['Biosample term id', 'Biosample term name',
                               'Biosample type', 'Biosample life stage',
                               'Biosample organ', 'Biosample cell']]
    new_meta = new_meta.drop_duplicates()
    new_meta.to_csv(os.path.join(path_out, 'metadata.dhs.tsv'), sep='\t',
                    index=None)
    list_input = []
    types = set(new_meta['Biosample type'].tolist())
    for biotype in types:
        type_meta = new_meta.loc[new_meta['Biosample type'] == biotype, :]
        path_type = os.path.join(path_out, biotype.replace(' ', '_'))
        if not os.path.exists(path_type):
            os.makedirs(path_type)
        life_stages = set(type_meta['Biosample life stage'].tolist())
        for life_stage in life_stages:
            life_meta = type_meta.loc[
                        type_meta['Biosample life stage'] == life_stage, :]
            path_life_stage = \
                os.path.join(path_type, life_stage.replace(' ', '_'))
            if not os.path.exists(path_life_stage):
                os.makedirs(path_life_stage)
            terms = set(life_meta['Biosample term name'].tolist())
            for term in terms:
                filter_meta = \
                    df_meta.loc[
                     (df_meta['Biosample type'] == biotype) &
                     (df_meta['Biosample life stage'] == life_stage) &
                     (df_meta['Biosample term name'] == term), :]
                accession_ids = filter_meta['File accession'].tolist()
                list_input.append(dict(path=path_life_stage,
                                       term_name=term,
                                       accession_ids=accession_ids))

    pool = Pool(processes=30)
    func_merge = partial(merge_bed, path_bed)
    pool.map(func_merge, list_input)
    pool.close()

    return


if __name__ == '__main__':
    time_start = time()
    # get bed file annotating protein-coding genes
    # gtf_file = \
    #     '/home/zy/driver_mutation/data/ENCODE/gencode.v19.annotation.gtf'
    # protein_file = \
    #     '/home/zy/driver_mutation/data/gene/genes.protein.gencode.v19.bed'
    # promoter_file = \
    #     '/home/zy/driver_mutation/data/gene/' \
    #     'promoters.up2k.protein.gencode.v19.bed'
    # with open(protein_file, 'w') as w_gene:
    #     with open(promoter_file, 'w') as w_pro:
    #         fmt_gene = "{chrom}\t{start}\t{end}\t{symbol}\t.\t{strand}\n"
    #         fmt_promoter = "{chrom}\t{start}\t{end}\t{symbol}\t.\t{strand}\n"
    #         with open(gtf_file, 'r') as r_gtf:
    #             for line_gene in r_gtf:
    #                 if line_gene[0] == '#':
    #                     continue
    #                 list_line_gene = line_gene.strip().split('\t')
    #                 if list_line_gene[2] != 'gene':
    #                     continue
    #                 list_attr = list_line_gene[8].strip().split('; ')
    #                 gene_type = list_attr[2][11:-1]
    #                 if list_attr[2][-15:-1] != "protein_coding":
    #                     continue
    #                 gene_name = list_attr[4][11:-1]
    #                 strand = list_line_gene[6]
    #                 dict_gene = dict(chrom=list_line_gene[0],
    #                                  start=list_line_gene[3],
    #                                  end=list_line_gene[4], symbol=gene_name,
    #                                  strand=strand)
    #                 w_gene.write(fmt_gene.format(**dict_gene))
    #                 if strand == '+':
    #                     pro_start = str(int(list_line_gene[3]) - 2000)
    #                     pro_end = list_line_gene[3]
    #                 elif strand == '-':
    #                     pro_start = list_line_gene[4]
    #                     pro_end = str(int(list_line_gene[4]) + 2000)
    #                 else:
    #                     print('Error')
    #                     break
    #                 dict_promoter = dict(chrom=list_line_gene[0],
    #                                      start=pro_start,
    #                                      end=pro_end, symbol=gene_name,
    #                                      strand=strand)
    #                 w_pro.write(fmt_promoter.format(**dict_promoter))

    # build life stage dictionary
    path_lifestage = '/home/zy/driver_mutation/data/ENCODE/metadata/life_stage'
    dict_lifestage = build_dict_attr(path_lifestage)

    # build organ dictionary
    path_organ = '/home/zy/driver_mutation/data/ENCODE/metadata/organ'
    dict_organ = build_dict_attr(path_organ)

    # build organ dictionary
    path_cell = '/home/zy/driver_mutation/data/ENCODE/metadata/cell'
    dict_cell = build_dict_attr(path_cell)

    # filter metafile
    path_dhs = \
        '/home/zy/driver_mutation/data/ENCODE/DNase-seq/hg19/bed_narrowpeak'
    ori_meta = os.path.join(path_dhs, 'metadata.tsv')
    # meta_adult = os.path.join(path_dhs, 'metadata_adult.tsv')
    # meta_embryo = os.path.join(path_dhs, 'metadata_embryo.tsv')
    # meta_unknown = os.path.join(path_dhs, 'metadata_unknown.tsv')
    # meta_newborn = os.path.join(path_dhs, 'metadata_newborn.tsv')
    # meta_child = os.path.join(path_dhs, 'metadata_child.tsv')
    # df_adult = hg19_filter(meta_adult)
    # df_adult['Biosample life stage'] = 'adult'
    # df_embryo = hg19_filter(meta_embryo)
    # df_embryo['Biosample life stage'] = 'embryonic'
    # df_unknown = hg19_filter(meta_unknown)
    # df_unknown['Biosample life stage'] = 'unknown'
    # df_newborn = hg19_filter(meta_newborn)
    # df_newborn['Biosample life stage'] = 'newborn'
    # df_child = hg19_filter(meta_child)
    # df_child['Biosample life stage'] = 'child'
    # df_merge = pd.concat(
    #     [df_adult, df_embryo, df_unknown, df_newborn, df_child])
    # meta_hg19 = os.path.join(path_dhs, 'metadata.hg19.tsv')
    # df_merge.to_csv(meta_hg19, sep='\t', index=None)

    df_meta_dhs = hg19_filter(ori_meta)
    df_meta_dhs = add_attr(df_meta_dhs, dict_lifestage, 'Biosample life stage')
    df_meta_dhs = add_attr(df_meta_dhs, dict_organ, 'Biosample organ')
    df_meta_dhs = add_attr(df_meta_dhs, dict_cell, 'Biosample cell')
    meta_hg19_dhs = os.path.join(path_dhs, 'metadata.hg19.tsv')
    df_meta_dhs.to_csv(meta_hg19_dhs, sep='\t', index=None)

    # unique DHS
    path_bed_dhs = \
        '/home/zy/driver_mutation/data/ENCODE/DNase-seq/hg19/bed_narrowpeak'
    path_dhs_hg19 = '/home/zy/driver_mutation/data/DHS/hg19/'
    unique_bed_files(path_bed_dhs, meta_hg19_dhs, path_dhs_hg19)

    # unique H3K27ac

    time_end = time()
    print(time_end - time_start)
