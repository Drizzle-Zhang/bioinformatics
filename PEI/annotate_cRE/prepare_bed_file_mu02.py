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
                    w_pro.write(fmt_promoter.format(**dict_promoter))

    return


def release_filter(meta_in):
    df_meta = pd.read_csv(meta_in, sep='\t')
    df_meta_released = df_meta.loc[
        df_meta['File Status'] == 'released',
        ['File accession', 'Experiment accession', 'Biosample term id',
         'Biosample term name', 'Biosample type', 'Biosample treatments',
         'Biosample genetic modifications methods', 'Assembly']]
    df_meta_released.index = df_meta_released['File accession']

    return df_meta_released


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


def sub_hg38tohg19(path_hg38, path_hg19, dict_in):
    file_hg38 = os.path.join(path_hg38, dict_in['File accession'] + '.bed')
    file_hg19 = os.path.join(path_hg19, dict_in['File accession'] + '.bed')
    file_chain = \
        '/lustre/tianlab/tools/files_liftOver/hg38ToHg19.over.chain.gz'
    file_ummap = os.path.join(
        path_hg19, dict_in['File accession'] + '.bed.unmap')
    if dict_in['Assembly'] == 'hg19':
        os.system(f"cp {file_hg38} {file_hg19}")
    if dict_in['Assembly'] == 'GRCh38':
        file_prefix = file_hg38 + '.prefix'
        file_suffix = file_hg38 + '.suffix'
        file_hg19_prefix = file_hg19 + '.prefix'
        os.system(f"cut -f 1,2,3,4 {file_hg38} > {file_prefix}")
        os.system(f"cut -f 4,5,6,7,8,9,10 {file_hg38} > {file_suffix}")
        os.system(
            f"liftOver {file_prefix} {file_chain} "
            f"{file_hg19_prefix} {file_ummap}")
        dict_peak_score = defaultdict(list)
        with open(file_suffix, 'r') as r_f:
            for line in r_f:
                list_line = line.strip().split('\t')
                dict_peak_score[list_line[0]].append(list_line[1:])
        with open(file_hg19, 'w') as w_f:
            fmt = "{chrom}\t{start}\t{end}\t{peak_id}\t{score}\t{strand}\t" \
                  "{fold_change}\t{p_value}\t{q_value}\t{peak_location}\n"
            with open(file_hg19_prefix, 'r') as r_hg19:
                for line in r_hg19:
                    list_line = line.strip().split('\t')
                    list_suffix = dict_peak_score[list_line[3]][0]
                    dict_hg19 = dict(
                        chrom=list_line[0], start=list_line[1],
                        end=list_line[2], peak_id=list_line[3],
                        score=list_suffix[0], strand=list_suffix[1],
                        fold_change=list_suffix[2], p_value=list_suffix[3],
                        q_value=list_suffix[4], peak_location=list_suffix[5]
                    )
                    w_f.write(fmt.format(**dict_hg19))
        os.remove(file_ummap)
        os.remove(file_prefix)
        os.remove(file_suffix)
        os.remove(file_hg19_prefix)

    return


def hg38tohg19(path_hg38, path_hg19):
    if os.path.exists(path_hg19):
        os.system(f"rm -rf {path_hg19}")
    os.mkdir(path_hg19)

    file_meta = os.path.join(path_hg38, 'metadata.simple.tsv')
    df_meta = pd.read_csv(file_meta, sep='\t')
    os.system(
        f"cp {file_meta} {os.path.join(path_hg19, 'metadata.simple.tsv')}")
    list_dict = df_meta.to_dict('records')

    pool = Pool(processes=60)
    func_hg38tohg19 = partial(sub_hg38tohg19, path_hg38, path_hg19)
    pool.map(func_hg38tohg19, list_dict)
    pool.close()

    return


def merge_bed(path_bed, col_collapse, dict_in):
    term_name = dict_in['term_name']
    path_out = dict_in['path']
    str_collapse = \
        ','.join([val for val in ['collapse']
                  for i in range(len(col_collapse.split(',')))])
    cat_out = os.path.join(path_out, f"{term_name}.cat.bed")
    sort_out = os.path.join(path_out, f"{term_name}.sort.bed")
    bed_out = os.path.join(path_out, f"{term_name}.bed")
    cat_in = ' '.join([os.path.join(path_bed, acce_id + '.bed')
                       for acce_id in dict_in['accession_ids']])
    os.system(f"cat {cat_in} > {cat_out}")
    os.system(f"bedtools sort -i {cat_out} > {sort_out}")
    # os.system(f"sort -k 1,1 -k2,2n {cat_out} > {sort_out}")
    os.system(f"bedtools merge -i {sort_out} "
              f"-c {col_collapse} -o {str_collapse} > {bed_out}")
    os.remove(cat_out)
    os.remove(sort_out)

    return


def ref_dhs(path_in, path_ref):
    df_meta = pd.read_csv(
        os.path.join(path_in, 'metadata.simple.tsv'), sep='\t')
    cancer_state = df_meta['Biosample cell'].apply(
        lambda x: True if isinstance(x, float) else
        'cancer cell' not in x.strip().split(',')
    )
    organ_state = df_meta['Biosample organ'].apply(
        lambda x: isinstance(x, str)
    )
    df_meta_normal = df_meta.loc[organ_state & cancer_state, :]

    if os.path.exists(path_ref):
        os.system(f"rm -rf {path_ref}")
    os.mkdir(path_ref)
    meta_out = df_meta_normal.loc[
               :, ['Biosample term name', 'Biosample life stage',
                   'Biosample organ']]
    meta_out = meta_out.drop_duplicates()
    meta_out.to_csv(os.path.join(path_ref, 'metadata.tsv'),
                    sep='\t', index=None)

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

    pool = Pool(processes=50)
    func_merge = partial(merge_bed, path_in, '5,6,7,8,9,10')
    pool.map(func_merge, list_input)
    pool.close()

    return


def unique_bed_files_histone(path_in, path_out):
    df_meta = pd.read_csv(
        os.path.join(path_in, 'metadata.simple.tsv'), sep='\t')
    cancer_state = df_meta['Biosample cell'].apply(
        lambda x: True if isinstance(x, float) else
        'cancer cell' not in x.strip().split(',')
    )
    organ_state = df_meta['Biosample organ'].apply(
        lambda x: isinstance(x, str)
    )
    df_meta_normal = df_meta.loc[organ_state & cancer_state, :]

    if os.path.exists(path_out):
        os.system(f"rm -rf {path_out}")
    os.mkdir(path_out)

    list_input = []
    # list_input = [dict(path=path_out,
    #                    term_name='all_organs',
    #                    accession_ids=
    #                    df_meta_normal['File accession'].tolist())]
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
        # list_input.append(dict(path=organ_path,
        #                        term_name=organ.replace(' ', '_'),
        #                        accession_ids=
        #                        organ_meta['File accession'].tolist()))
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
            # list_input.append(dict(path=path_life_stage,
            #                        term_name=life_stage.replace(' ', '_'),
            #                        accession_ids=
            #                        life_meta['File accession'].tolist()))
            terms = set(life_meta['Biosample term name'].tolist())
            for term in terms:
                filter_meta = \
                    life_meta.loc[life_meta['Biosample term name'] == term, :]
                accession_ids = filter_meta['File accession'].tolist()
                path_term = \
                    os.path.join(path_life_stage, term.replace(
                        ' ', '_').replace('/', '+').replace("'", '--'))
                if not os.path.exists(path_term):
                    os.makedirs(path_term)
                list_input.append(dict(path=path_term,
                                       term_name=
                                       term.replace(' ', '_').replace(
                                           '/', '+').replace("'", '--'),
                                       accession_ids=accession_ids))

    pool = Pool(processes=50)
    func_merge = partial(merge_bed, path_in, '5,6,7,8,9,10')
    pool.map(func_merge, list_input)
    pool.close()

    return


if __name__ == '__main__':
    time_start = time()
    # get bed file annotating protein-coding genes
    gtf_file_hg19 = \
        '/lustre/tianlab/zhangyu/driver_mutation/data/ENCODE/' \
        'gencode.v19.annotation.gtf'
    protein_file_hg19 = \
        '/lustre/tianlab/zhangyu/driver_mutation/data/gene/' \
        'genes.protein.gencode.v19.bed'
    promoter_file_hg19 = \
        '/lustre/tianlab/zhangyu/driver_mutation/data/gene/' \
        'promoters.up2k.protein.gencode.v19.bed'
    # generate_gene_file(gtf_file_hg19, protein_file_hg19, promoter_file_hg19)

    # build life stage dictionary
    path_lifestage = '/lustre/tianlab/zhangyu/driver_mutation/data/' \
                     'ENCODE/metadata/life_stage'
    dict_lifestage = build_dict_attr(path_lifestage)

    # build organ dictionary
    path_organ = '/lustre/tianlab/zhangyu/driver_mutation/data/ENCODE/' \
                 'metadata/organ'
    dict_organ = build_dict_attr(path_organ)

    # build organ dictionary
    path_cell = '/lustre/tianlab/zhangyu/driver_mutation/data/ENCODE/' \
                'metadata/cell'
    dict_cell = build_dict_attr(path_cell)

    # DHS
    # metafile
    path_dhs = \
        '/lustre/tianlab/zhangyu/driver_mutation/data/ENCODE/DNase-seq/all'
    ori_meta_dhs = os.path.join(path_dhs, 'metadata.tsv')
    df_meta_dhs = release_filter(ori_meta_dhs)
    df_meta_dhs = add_attr(df_meta_dhs, dict_lifestage, 'Biosample life stage')
    df_meta_dhs = add_attr(df_meta_dhs, dict_organ, 'Biosample organ')
    df_meta_dhs = add_attr(df_meta_dhs, dict_cell, 'Biosample cell')
    meta_dhs = os.path.join(path_dhs, 'metadata.simple.tsv')
    df_meta_dhs.to_csv(meta_dhs, sep='\t', index=None)

    # hg38 to hg19
    path_hg38tohg19 = \
        '/lustre/tianlab/zhangyu/driver_mutation/data/ENCODE/' \
        'DNase-seq/GRCh38tohg19'
    hg38tohg19(path_dhs, path_hg38tohg19)

    # build DHS reference
    path_dhs_hg38tohg19 = \
        '/lustre/tianlab/zhangyu/driver_mutation/data/DHS/GRCh38tohg19/'
    ref_dhs(path_hg38tohg19, path_dhs_hg38tohg19)

    # H3K27ac
    path_h3k27ac = \
        '/lustre/tianlab/zhangyu/driver_mutation/data/ENCODE/' \
        'histone_ChIP-seq/H3K27ac'
    ori_meta_h3k27ac = os.path.join(path_h3k27ac, 'metadata.tsv')
    df_meta_h3k27ac = release_filter(ori_meta_h3k27ac)
    df_meta_h3k27ac = \
        add_attr(df_meta_h3k27ac, dict_lifestage, 'Biosample life stage')
    df_meta_h3k27ac = add_attr(df_meta_h3k27ac, dict_organ, 'Biosample organ')
    df_meta_h3k27ac = add_attr(df_meta_h3k27ac, dict_cell, 'Biosample cell')
    meta_h3k27ac = os.path.join(path_h3k27ac, 'metadata.simple.tsv')
    df_meta_h3k27ac.to_csv(meta_h3k27ac, sep='\t', index=None)

    # hg38 to hg19
    path_hg38tohg19 = \
        '/lustre/tianlab/zhangyu/driver_mutation/data/ENCODE/' \
        'histone_ChIP-seq/GRCh38tohg19/H3K27ac'
    hg38tohg19(path_h3k27ac, path_hg38tohg19)

    # unique H3K27ac
    path_h3k27ac_hg38tohg19 = \
        '/lustre/tianlab/zhangyu/driver_mutation/data/ENCODE/' \
        'histone_ChIP-seq/GRCh38tohg19/H3K27ac_merge'
    unique_bed_files_histone(path_hg38tohg19, path_h3k27ac_hg38tohg19)
    """
    """
    # H3K4me3
    path_h3k4me3 = \
        '/lustre/tianlab/zhangyu/driver_mutation/data/ENCODE/' \
        'histone_ChIP-seq/H3K4me3'
    ori_meta_h3k4me3 = os.path.join(path_h3k4me3, 'metadata.tsv')
    df_meta_h3k4me3 = release_filter(ori_meta_h3k4me3)
    df_meta_h3k4me3 = \
        add_attr(df_meta_h3k4me3, dict_lifestage, 'Biosample life stage')
    df_meta_h3k4me3 = add_attr(df_meta_h3k4me3, dict_organ, 'Biosample organ')
    df_meta_h3k4me3 = add_attr(df_meta_h3k4me3, dict_cell, 'Biosample cell')
    meta_h3k4me3 = os.path.join(path_h3k4me3, 'metadata.simple.tsv')
    df_meta_h3k4me3.to_csv(meta_h3k4me3, sep='\t', index=None)

    # hg38 to hg19
    path_hg38tohg19 = \
        '/lustre/tianlab/zhangyu/driver_mutation/data/ENCODE/' \
        'histone_ChIP-seq/GRCh38tohg19/H3K4me3'
    hg38tohg19(path_h3k4me3, path_hg38tohg19)

    # unique H3K27ac
    path_h3k4me3_hg38tohg19 = \
        '/lustre/tianlab/zhangyu/driver_mutation/data/ENCODE/' \
        'histone_ChIP-seq/GRCh38tohg19/H3K4me3_merge'
    unique_bed_files_histone(path_hg38tohg19, path_h3k4me3_hg38tohg19)

    time_end = time()
    print(time_end - time_start)
