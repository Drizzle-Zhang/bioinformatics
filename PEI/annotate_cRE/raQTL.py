#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: raQTL.py
# @time: 1/14/20 4:26 PM

from time import time
from compare_hic_V1 import transform_eqtl
import pandas as pd
import os

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
    file_egenes = '/local/zy/PEI/raQTL/liver/Liver.v7.egenes.txt'
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

    time_end = time()
    print(time_end - time_start)
