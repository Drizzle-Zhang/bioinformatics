#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: anls.py
# @time: 11/24/19 5:14 PM

from time import time
import pandas as pd
import numpy as np
import os


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
                    gene_type = list_attr[1][11:-1]
                    if gene_type != "protein_coding":
                        continue
                    gene_name = list_attr[3][11:-1]
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


def filter_meta(path_in, list_exp):
    file_meta = os.path.join(path_in, 'metadata.tsv')
    df_meta = pd.read_csv(file_meta, sep='\t')
    bool_exp = df_meta['Experiment accession'].apply(lambda x: x in list_exp)
    df_simple = df_meta.loc[
        (df_meta['File Status'] == 'released') & bool_exp,
        ['File accession', 'Experiment accession', 'Output type']]
    df_simple.to_csv(
        os.path.join(path_in, 'metadata.simple.tsv'), sep='\t'
    )

    return


def annotate_cre(path_root):
    df_ref = pd.read_csv(os.path.join(path_root, 'list_ref.txt'), sep='\t')
    for i in range(df_ref.shape[0]):
        path_out = os.path.join(
            path_root, (df_ref.loc[i, 'Biosample name']).replace(
                ' ', '_').replace('/', '-').replace('(', '').replace(')', '')
        )
        if not os.path.exists(path_out):
            os.mkdir(path_out)
        # DHS
        path_dhs = os.path.join(path_root, 'DNase-seq')
        df_dhs = pd.read_csv(os.path.join(path_dhs, 'metadata.tsv'), sep='\t')
        exp_dhs = df_ref.loc[i, 'Exp_DHS']
        list_access = (df_dhs.loc[
            (df_dhs['File Status'] == 'released') &
            (df_dhs['Experiment accession'] == exp_dhs),
            'File accession']).tolist()
        cat_out = os.path.join(path_out, f"DHS.bed.cat")
        sort_out = os.path.join(path_out, f"DHS.bed.sort")
        merge_out = os.path.join(path_out, f"DHS.bed")
        cat_in = ' '.join([os.path.join(path_dhs, acce_id + '.bed')
                           for acce_id in list_access])
        os.system(f"cat {cat_in} > {cat_out}")
        os.system(f"bedtools sort -i {cat_out} > {sort_out}")
        os.system(f"bedtools merge -i {sort_out} > {merge_out}")
        os.remove(cat_out)
        os.remove(sort_out)

        # annotate promoter locations
        path_promoter_in = os.path.join(
            path_root, 'promoters.2k.protein.gencode.vM7.bed'
        )
        bedtools_out = os.path.join(path_out, "bedtools.intersect.out")
        path_promoter = os.path.join(path_out, 'promoter.bed')
        os.system(
            f"bedtools intersect -a {merge_out} -b {path_promoter_in} -loj "
            f"> {bedtools_out}")
        os.system(f"cut -f 1,2,3,7 {bedtools_out} > {path_promoter}")
        os.remove(bedtools_out)

        # H3K4me3
        path_h3k4me3 = os.path.join(path_root, 'H3K4me3')
        df_h3k4me3 = pd.read_csv(
            os.path.join(path_h3k4me3, 'metadata.tsv'), sep='\t'
        )
        exp_h3k4me3 = df_ref.loc[i, 'Exp_H3K4me3']
        file_h3k4me3 = os.path.join(
            path_h3k4me3, (df_h3k4me3.loc[
                (df_h3k4me3['File Status'] == 'released') &
                (df_h3k4me3['Experiment accession'] == exp_h3k4me3) &
                (df_h3k4me3['Output type'] == 'replicated peaks'),
                'File accession']).tolist()[0] + '.bed'
        )
        os.system(
            f"bedtools intersect -a {path_promoter} -b {file_h3k4me3} -wao "
            f"> {bedtools_out}"
        )
        path_h3k4me3_out = os.path.join(path_out, 'H3K4me3.bed')
        os.system(
            f"cut -f 1,2,3,4,5,6,7,12,15 {bedtools_out} > {path_h3k4me3_out}"
        )
        os.remove(bedtools_out)

        # H3K27ac
        path_h3k27ac = os.path.join(path_root, 'H3K27ac')
        df_h3k27ac = pd.read_csv(
            os.path.join(path_h3k27ac, 'metadata.tsv'), sep='\t'
        )
        exp_h3k27ac = df_ref.loc[i, 'Exp_H3K27ac']
        file_h3k27ac = os.path.join(
            path_h3k27ac, (df_h3k27ac.loc[
                (df_h3k27ac['File Status'] == 'released') &
                (df_h3k27ac['Experiment accession'] == exp_h3k27ac) &
                (df_h3k27ac['Output type'] == 'replicated peaks'),
                'File accession']).tolist()[0] + '.bed'
        )
        os.system(
            f"bedtools intersect -a {path_h3k4me3_out} -b {file_h3k27ac} -wao "
            f"> {bedtools_out}"
        )
        path_h3k27ac_out = os.path.join(path_out, 'H3K27ac.bed')
        os.system(
            f"cut -f 1,2,3,4,5,6,7,8,9,10,11,12,17,20 {bedtools_out} "
            f"> {path_h3k27ac_out}")
        os.remove(bedtools_out)

        # annotate cRE
        file_out = os.path.join(path_out, 'cRE.bed')
        file_out_tmp = os.path.join(path_out, "cRE.tmp")
        file_merge = os.path.join(path_out, "cRE.merge")
        with open(path_h3k27ac_out, 'r') as r_f:
            with open(file_out_tmp, 'w') as w_f:
                fmt_dhs = "{chrom}\t{start}\t{end}\t{dhs_id}\t{score}\t.\t" \
                          "{cre}\t{overlap}\n"
                list_promoter = []
                old_promoter = ''
                for line in r_f:
                    list_line = line.strip().split('\t')
                    chrom = list_line[0]
                    start = list_line[1]
                    end = list_line[2]
                    dhs_id = f"DHS<-{chrom}:{start}-{end}"
                    overlap = int(list_line[13])
                    promoter = list_line[3]
                    h3k4me3 = list_line[4]
                    h3k27ac = list_line[9]
                    if (promoter != '.') & (h3k4me3 != '.'):
                        if (list_promoter != []) & (promoter != old_promoter):
                            df_pro = pd.DataFrame(list_promoter)
                            gene_name = (df_pro['promoter'].tolist()[0]).split(
                                '<-')[0]
                            score = np.sum(
                                df_pro['score'] * df_pro['overlap'])/np.sum(
                                df_pro['overlap'])
                            dict_pro = \
                                dict(chrom=df_pro.loc[0, 'chrom'],
                                     start=np.min(df_pro['start']),
                                     end=np.max(df_pro['end']),
                                     dhs_id=gene_name,
                                     score=score, cre='Promoter',
                                     overlap=np.sum(df_pro['overlap']))
                            w_f.write(fmt_dhs.format(**dict_pro))
                            list_promoter = []
                        list_promoter.append(
                            dict(chrom=chrom, start=start, end=end,
                                 score=float(list_line[7]), promoter=promoter,
                                 overlap=int(list_line[8])))
                        old_promoter = promoter
                        continue
                    elif h3k27ac != '.':
                        if list_promoter:
                            df_pro = pd.DataFrame(list_promoter)
                            gene_name = (df_pro['promoter'].tolist()[0]).split(
                                '<-')[0]
                            score = np.sum(
                                df_pro['score'] * df_pro['overlap'])/np.sum(
                                df_pro['overlap'])
                            dict_pro = \
                                dict(chrom=df_pro.loc[0, 'chrom'],
                                     start=np.min(df_pro['start']),
                                     end=np.max(df_pro['end']),
                                     dhs_id=gene_name,
                                     score=score, cre='Promoter',
                                     overlap=np.sum(df_pro['overlap']))
                            w_f.write(fmt_dhs.format(**dict_pro))
                        cre = 'Enhancer'
                        score = list_line[12]
                        list_promoter = []
                        old_promoter = ''
                    else:
                        if list_promoter:
                            df_pro = pd.DataFrame(list_promoter)
                            gene_name = (df_pro['promoter'].tolist()[0]).split(
                                '<-')[0]
                            score = np.sum(
                                df_pro['score'] * df_pro['overlap'])/np.sum(
                                df_pro['overlap'])
                            dict_pro = \
                                dict(chrom=df_pro.loc[0, 'chrom'],
                                     start=np.min(df_pro['start']),
                                     end=np.max(df_pro['end']),
                                     dhs_id=gene_name,
                                     score=score, cre='Promoter',
                                     overlap=np.sum(df_pro['overlap']))
                            w_f.write(fmt_dhs.format(**dict_pro))
                        list_promoter = []
                        old_promoter = ''
                        continue
                    w_f.write(fmt_dhs.format(**locals()))

        os.system(
            f"bedtools merge -i {file_out_tmp} -c 4,5,7,8 "
            f"-o collapse,collapse,collapse,collapse > {file_merge}"
        )
        with open(file_merge, 'r') as r_f:
            with open(file_out, 'w') as w_f:
                fmt_dhs = "{chrom}\t{start}\t{end}\t{dhs_id}\t{score}\t.\t" \
                          "{cre}\n"
                for line in r_f:
                    list_line = line.strip().split('\t')
                    chrom = list_line[0]
                    start = list_line[1]
                    end = list_line[2]
                    if list_line[3][:3] == 'DHS':
                        dhs_id = list(set(list_line[3].strip().split(',')))[0]
                        array_score = np.array(
                            [float(num) for num in
                             list_line[4].strip().split(',')]
                        )
                        array_overlap = np.array(
                            [int(num) for num in
                             list_line[6].strip().split(',')]
                        )
                        score = np.sum(
                            array_score * array_overlap)/np.sum(array_overlap)
                    else:
                        dhs_id = list_line[3]
                        score = list(set(list_line[4].strip().split(',')))[0]
                    cre = list(set(list_line[5].strip().split(',')))[0]
                    w_f.write(fmt_dhs.format(**locals()))

    return


if __name__ == '__main__':
    time_start = time()
    # get bed file annotating protein-coding genes
    gtf_file_mm10 = \
        '/home/zy/data_processing/anls/gencode.vM7.annotation.gtf'
    protein_file_mm10 = \
        '/home/zy/data_processing/anls/genes.protein.gencode.vM7.bed'
    promoter_file_mm10 = \
        '/home/zy/data_processing/anls/promoters.2k.protein.gencode.vM7.bed'
    generate_gene_file(gtf_file_mm10, protein_file_mm10, promoter_file_mm10)

    path_anls = '/home/zy/data_processing/anls'
    annotate_cre(path_anls)

    time_end = time()
    print(time_end - time_start)
