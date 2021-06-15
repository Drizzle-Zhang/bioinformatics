#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: raQTL.py
# @time: 1/14/20 4:26 PM

from time import time
from compare_hic_V1 import transform_eqtl
import pandas as pd
import numpy as np
import os


def unify_raqtl(file_in, file_out):
    file_tmp = file_out + '.tmp'
    with open(file_tmp, 'w') as w_f:
        fmt = "{chrom}\t{start}\t{end}\t{variant_id}\n"
        with open(file_in, 'r') as r_in:
            for line in r_in:
                list_line = line.strip().split('\t')
                if list_line[0] == 'chr':
                    continue
                chrom = list_line[0]
                start = list_line[2]
                end = str(int(list_line[2]) + 1)
                variant_id = list_line[1]
                w_f.write(fmt.format(**locals()))

    os.system(f"bedtools sort -i {file_tmp} > {file_out}")
    os.remove(file_tmp)

    return


def drop_dup(x):
    if x.shape[0] == 1:
        return x
    else:
        max_overlap = np.max(x.iloc[:, -1])
        row_out = x.loc[x.iloc[:, -1] == max_overlap, :]
        return row_out


def generate_test_file(file_in_eqtl, file_in_snp, file_in_raqtl,
                       file_out, file_inner):
    file_intersect = file_out + '.intersect'
    os.system(f"bedtools intersect -a {file_in_snp} -b {file_in_eqtl} -wao "
              f"> {file_intersect}")
    df_intersect = pd.read_csv(file_intersect, sep='\t', header=None)
    df_intersect = df_intersect.groupby(3).apply(drop_dup)
    df_intersect.columns = \
        ['chrom_snp', 'start_snp', 'end_snp', 'snp_id',
         'chrom_eqtl', 'start_eqtl', 'end_eqtl', 'snp_eqtl', 'overlap']
    df_raqtl = pd.read_csv(file_in_raqtl, sep='\t')
    set_raqtl = set(df_raqtl['SNP_ID'].to_list())
    list_raqtl = df_intersect['snp_id'].apply(
        lambda x: 1 if x in set_raqtl else 0)
    df_intersect['raQTL'] = list_raqtl

    df_out = df_intersect.loc[:, ['snp_id', 'snp_eqtl', 'raQTL']]
    df_out.to_csv(file_out, sep='\t', index=None)
    df_out_inner = df_out.loc[
                   (df_out['snp_eqtl'] != '.') & (df_out['raQTL'] != 0), :]
    df_out_inner.to_csv(file_inner, sep='\t', index=None)

    return


if __name__ == '__main__':
    time_start = time()
    # cRE file
    file_cre_liver_origin = '/local/zy/PEI/data/DHS/cRE_annotation/liver/' \
                            'adult_liver/cRE.txt'
    file_cre_liver = '/local/zy/PEI/raQTL/liver/cRE.txt'
    protein_exon = '/local/zy/PEI/data/gene/exon.up2k.protein.gencode.v19.bed'
    os.system(f"bedtools intersect -a {file_cre_liver_origin} "
              f"-b {protein_exon} -v > {file_cre_liver}")

    # eQTL
    # eGenes
    file_egenes = '/local/zy/PEI/raQTL/liver/eQTL/Liver.v7.egenes.txt'
    protein_gene = '/local/zy/PEI/data/gene/genes.protein.gencode.v19.bed'
    df_protein_genes = pd.read_csv(protein_gene, sep='\t', header=None)
    df_protein_genes.columns = ['chrom', 'start', 'end',
                                'gene_name', 'gene_id', 'strand']
    df_egenes = pd.read_csv(file_egenes, sep='\t')
    df_egenes = df_egenes.loc[
        df_egenes['qval'] <= 0.05,
        ['gene_id', 'gene_name', 'gene_start', 'gene_end']
    ]
    df_egenes_pos = pd.merge(
        df_protein_genes, df_egenes, on=['gene_id', 'gene_name']
    )
    df_egenes_pos = df_egenes_pos.loc[
                    :, ['chrom', 'gene_start', 'gene_end', 'gene_name',
                        'gene_id', 'strand']
                    ]

    file_eqtl = '/local/zy/PEI/raQTL/liver/eQTL/' \
                'Liver.v7.signif_variant_gene_pairs.txt'
    file_uniform = '/local/zy/PEI/raQTL/liver/eQTL/uniform_1000.txt'
    pairs_eqtl_egene = '/local/zy/PEI/raQTL/liver/eQTL/eQTL_cRE_pairs_1000.txt'
    pairs_dhs_egene = '/local/zy/PEI/raQTL/liver/eQTL/cRE_pairs_1000.txt'
    transform_eqtl(
        file_eqtl, file_uniform, pairs_eqtl_egene, pairs_dhs_egene,
        file_cre_liver, df_egenes_pos)
    df_uniform = pd.read_csv(file_uniform, sep='\t', header=None,
                             usecols=[0, 1, 2, 3])
    df_uniform = df_uniform.drop_duplicates()
    df_uniform.index = df_uniform[3]
    df_eqtl = pd.read_csv(pairs_eqtl_egene, sep='\t', header=None, usecols=[0])
    df_eqtl = ((df_eqtl.drop_duplicates()).iloc[:, 0]).tolist()
    df_eqtl = df_uniform.loc[df_eqtl, :]
    file_filterd_eqtl = '/local/zy/PEI/raQTL/liver/eQTL_liver.txt'
    df_eqtl.to_csv(file_filterd_eqtl, sep='\t', header=None, index=None)

    # raQTL
    raqtl = '/local/zy/PEI/raQTL/liver/hepg2.sign.id.LP190708.txt'
    snps = '/local/zy/PEI/raQTL/liver/SuRE_SNP_table_LP190708.txt'
    uniform_snps = '/local/zy/PEI/raQTL/liver/SuRE_SNP.uniform.txt'
    unify_raqtl(snps, uniform_snps)

    # fisher test
    feqtl_raqtl = '/local/zy/PEI/raQTL/liver/filterd_eqtl_raqtl.txt'
    feqtl_raqtl_inner = \
        '/local/zy/PEI/raQTL/liver/filterd_eqtl_raqtl_inner.txt'
    generate_test_file(
        file_filterd_eqtl, uniform_snps, raqtl, feqtl_raqtl, feqtl_raqtl_inner)

    time_end = time()
    print(time_end - time_start)
