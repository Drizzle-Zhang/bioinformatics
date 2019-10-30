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
    term_name = dict_in['term_name']
    path_out = dict_in['path']
    cat_out = os.path.join(path_out, f"{term_name}.cat.bed")
    sort_out = os.path.join(path_out, f"{term_name}.sort.bed")
    bed_out = os.path.join(path_out, f"{term_name}.bed")
    cat_in = ' '.join([os.path.join(path_bed, acce_id + '.bed')
                       for acce_id in dict_in['accession_ids']])
    os.system(f"cat {cat_in} > {cat_out}")
    os.system(f"bedtools sort -i {cat_out} > {sort_out}")
    os.system(f"bedtools merge -i {sort_out} "
              f"-c 5,6,7,8,9,10 "
              f"-o collapse,collapse,collapse,collapse,collapse,"
              f"collapse > {bed_out}")
    os.remove(cat_out)
    os.remove(sort_out)

    return


def ref_dhs(path_in, path_ref):
    df_meta = pd.read_csv(os.path.join(path_in, 'metadata.hg19.tsv'), sep='\t')
    cancer_state = df_meta['Biosample cell'].apply(
        lambda x: True if isinstance(x, float) else
        'cancer cell' not in x.strip().split(',')
    )
    organ_state = df_meta['Biosample organ'].apply(
        lambda x: isinstance(x, str)
    )
    df_meta_normal = df_meta.loc[organ_state & cancer_state, :]

    os.system(f"rm -rf {path_ref}")

    list_input = [dict(path=path_ref,
                       term_name='all_organs',
                       accession_ids=
                       df_meta_normal['File accession'].tolist())]
    organs = []
    for line in df_meta_normal['Biosample organ'].tolist():
        organs.extend(line.strip().split(','))
    organs = set(organs)
    for organ in organs:
        organ_meta = df_meta_normal.loc[
                     df_meta_normal['Biosample organ'].apply(
                         lambda x: organ in x.strip().split(',')), :]
        organ_path = \
            os.path.join(path_ref, organ.replace(' ', '_'))
        if not os.path.exists(organ_path):
            os.makedirs(organ_path)
        list_input.append(dict(path=organ_path,
                               term_name=organ.replace(' ', '_'),
                               accession_ids=
                               organ_meta['File accession'].tolist()))
        life_stages = []
        for line in organ_meta['Biosample life stage'].tolist():
            life_stages.extend(line.strip().split(','))
        life_stages = set(life_stages)
        for life_stage in life_stages:
            filter_meta = \
                organ_meta.loc[
                 (organ_meta['Biosample life stage'].apply(
                     lambda x: life_stage in x.strip().split(','))) &
                 (organ_meta['Biosample organ'].apply(
                         lambda x: organ in x.strip().split(','))), :]
            accession_ids = filter_meta['File accession'].tolist()
            list_input.append(dict(path=organ_path,
                                   term_name=life_stage.replace(' ', '_'),
                                   accession_ids=accession_ids))

    pool = Pool(processes=40)
    func_merge = partial(merge_bed, path_in)
    pool.map(func_merge, list_input)
    pool.close()

    return


def unique_bed_files_histone(path_in, path_out):
    df_meta = pd.read_csv(os.path.join(path_in, 'metadata.hg19.tsv'), sep='\t')
    cancer_state = df_meta['Biosample cell'].apply(
        lambda x: True if isinstance(x, float) else
        'cancer cell' not in x.strip().split(',')
    )
    organ_state = df_meta['Biosample organ'].apply(
        lambda x: isinstance(x, str)
    )
    df_meta_normal = df_meta.loc[organ_state & cancer_state, :]

    list_input = [dict(path=path_out,
                       term_name='all_organs',
                       accession_ids=
                       df_meta_normal['File accession'].tolist())]
    organs = []
    for line in df_meta_normal['Biosample organ'].tolist():
        organs.extend(line.strip().split(','))
    organs = set(organs)
    for organ in organs:
        organ_meta = df_meta_normal.loc[
                     df_meta_normal['Biosample organ'].apply(
                         lambda x: organ in x.strip().split(',')), :]
        organ_path = \
            os.path.join(path_out, organ.replace(' ', '_'))
        if not os.path.exists(organ_path):
            os.makedirs(organ_path)
        list_input.append(dict(path=path_out,
                               term_name=organ.replace(' ', '_'),
                               accession_ids=
                               organ_meta['File accession'].tolist()))
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
            list_input.append(dict(path=path_life_stage,
                                   term_name=life_stage.replace(' ', '_'),
                                   accession_ids=
                                   life_meta['File accession'].tolist()))
            terms = set(life_meta['Biosample file name'].tolist())
            for term in terms:
                filter_meta = \
                    life_meta.loc[(df_meta['Biosample file name'] == term), :]
                accession_ids = filter_meta['File accession'].tolist()
                list_input.append(dict(path=path_life_stage,
                                       term_name=
                                       term.replace(' ', '_').replace(
                                           '/', '+').replace("'", '--'),
                                       accession_ids=accession_ids))

    pool = Pool(processes=40)
    func_merge = partial(merge_bed, path_in)
    pool.map(func_merge, list_input)
    pool.close()

    return


if __name__ == '__main__':
    time_start = time()
    # get bed file annotating protein-coding genes
    gtf_file = \
        '/home/zy/driver_mutation/data/ENCODE/gencode.v19.annotation.gtf'
    protein_file = \
        '/home/zy/driver_mutation/data/gene/genes.protein.gencode.v19.bed'
    promoter_file = \
        '/home/zy/driver_mutation/data/gene/' \
        'promoters.up2k.protein.gencode.v19.bed'
    """
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
                    if list_attr[2][-15:-1] != "protein_coding":
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
                        pro_end = list_line_gene[3]
                    elif strand == '-':
                        pro_start = list_line_gene[4]
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
                    w_pro.write(fmt_promoter.format(**dict_promoter))"""

    # build life stage dictionary
    path_lifestage = '/home/zy/driver_mutation/data/ENCODE/metadata/life_stage'
    dict_lifestage = build_dict_attr(path_lifestage)

    # build organ dictionary
    path_organ = '/home/zy/driver_mutation/data/ENCODE/metadata/organ'
    dict_organ = build_dict_attr(path_organ)

    # build organ dictionary
    path_cell = '/home/zy/driver_mutation/data/ENCODE/metadata/cell'
    dict_cell = build_dict_attr(path_cell)

    # DHS
    # metafile
    path_dhs = \
        '/home/zy/driver_mutation/data/ENCODE/DNase-seq/hg19/bed_narrowpeak'
    ori_meta_dhs = os.path.join(path_dhs, 'metadata.tsv')
    df_meta_dhs = hg19_filter(ori_meta_dhs)
    df_meta_dhs = add_attr(df_meta_dhs, dict_lifestage, 'Biosample life stage')
    df_meta_dhs = add_attr(df_meta_dhs, dict_organ, 'Biosample organ')
    df_meta_dhs = add_attr(df_meta_dhs, dict_cell, 'Biosample cell')
    meta_hg19_dhs = os.path.join(path_dhs, 'metadata.hg19.tsv')
    df_meta_dhs.to_csv(meta_hg19_dhs, sep='\t', index=None)

    # build DHS reference
    path_dhs_hg19 = '/home/zy/driver_mutation/data/DHS/hg19/'
    ref_dhs(path_dhs, path_dhs_hg19)

    # H3K27ac
    path_h3k27ac = \
        '/home/zy/driver_mutation/data/ENCODE/histone_ChIP-seq/' \
        'hg19/bed_narrowpeak'
    ori_meta_h3k27ac = os.path.join(path_h3k27ac, 'metadata.tsv')
    df_meta_h3k27ac = hg19_filter(ori_meta_h3k27ac)
    df_meta_h3k27ac = \
        add_attr(df_meta_h3k27ac, dict_lifestage, 'Biosample life stage')
    df_meta_h3k27ac = add_attr(df_meta_h3k27ac, dict_organ, 'Biosample organ')
    df_meta_h3k27ac = add_attr(df_meta_h3k27ac, dict_cell, 'Biosample cell')
    meta_hg19_h3k27ac = os.path.join(path_h3k27ac, 'metadata.hg19.tsv')
    df_meta_h3k27ac.to_csv(meta_hg19_h3k27ac, sep='\t', index=None)

    # unique H3K27ac
    path_h3k27ac_hg19 = \
        '/home/zy/driver_mutation/data/ENCODE/histone_ChIP-seq/hg19/H3K27ac'
    unique_bed_files_histone(path_h3k27ac, path_h3k27ac_hg19)

    time_end = time()
    print(time_end - time_start)
