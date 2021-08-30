#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: GWAS_analysis.py
# @time: 2021/6/15 12:33

from time import time
import os
import pandas as pd

path_root = '/mdshare/node9/zy/Brain_GWAS/'

# ace DHS
file_acedhs = os.path.join(path_root, 'DHS_gene.bed')
file_alldhs = os.path.join(path_root, 'DHS.pos')

# GWAS sites hg19
file_snp_bed = os.path.join(path_root, 'Brain_SNP_hg19.bed')
file_snp_bed_sort = os.path.join(path_root, 'Brain_SNP_hg19.sort.bed')
os.system(f"bedtools sort -i {file_snp_bed} > {file_snp_bed_sort}")
file_snp_bed_1k = os.path.join(path_root, 'Brain_SNP_hg19_1kb.bed')
file_snp_bed_sort_1k = os.path.join(path_root, 'Brain_SNP_hg19_1kb.sort.bed')
os.system(f"bedtools sort -i {file_snp_bed_1k} > {file_snp_bed_sort_1k}")

# single site
file_acedhs_snp = os.path.join(path_root, 'acedhs_SNP.pos')
file_alldhs_snp = os.path.join(path_root, 'allDHS_SNP.pos')
# acedhs - SNP
os.system(f"bedtools intersect -a {file_acedhs} -b {file_snp_bed_sort} -wo "
          f"> {file_acedhs_snp}")
# allDHS - SNP
os.system(f"bedtools intersect -a {file_alldhs} -b {file_snp_bed_sort} -wo "
          f"> {file_alldhs_snp}")
df_acedhs = pd.read_csv(file_acedhs_snp, sep='\t', header=None)
df_alldhs = pd.read_csv(file_alldhs_snp, sep='\t', header=None)
df_acedhs.shape[0]/3538 # 0.00141
df_alldhs.shape[0]/1867665 # 0.00091

# +- 1kb
file_acedhs_snp_1k = os.path.join(path_root, 'acedhs_SNP_1k.pos')
file_alldhs_snp_1k = os.path.join(path_root, 'allDHS_SNP_1k.pos')
# acedhs - SNP
os.system(f"bedtools intersect -a {file_acedhs} -b {file_snp_bed_sort_1k} -wo "
          f"> {file_acedhs_snp_1k}")
# allDHS - SNP
os.system(f"bedtools intersect -a {file_alldhs} -b {file_snp_bed_sort_1k} -wo "
          f"> {file_alldhs_snp_1k}")
df_acedhs_1k = pd.read_csv(file_acedhs_snp_1k, sep='\t', header=None)
df_alldhs_1k = pd.read_csv(file_alldhs_snp_1k, sep='\t', header=None)
df_acedhs_1k.shape[0]/3538 # 0.00876
df_alldhs_1k.shape[0]/1867665 # 0.00721


# STR
file_str_human = os.path.join(path_root, 'human_expan.uniq.bed')
file_str_nhp = os.path.join(path_root, 'nhp_expan.uniq.bed')
file_str_all = os.path.join(path_root, 'All_ortholog_STR.txt')

# GWAS site hg38
file_snp_bed = os.path.join(path_root, 'Brain_SNP_hg38.bed')
file_snp_bed_sort = os.path.join(path_root, 'Brain_SNP_hg38.sort.bed')
os.system(f"bedtools sort -i {file_snp_bed} > {file_snp_bed_sort}")
file_snp_bed_1k = os.path.join(path_root, 'Brain_SNP_hg38_1kb.bed')
file_snp_bed_sort_1k = os.path.join(path_root, 'Brain_SNP_hg38_1kb.sort.bed')
os.system(f"bedtools sort -i {file_snp_bed_1k} > {file_snp_bed_sort_1k}")

# single STR
file_str_snp = os.path.join(path_root, 'Human_STR_SNP.pos')
file_allstr_snp = os.path.join(path_root, 'All_STR_SNP.pos')
# expand STR - SNP
os.system(f"bedtools intersect -a {file_str_human} -b {file_snp_bed_sort} -wo "
          f"> {file_str_snp}")
# all STR - SNP
os.system(f"bedtools intersect -a {file_str_all} -b {file_snp_bed_sort} -wo "
          f"> {file_allstr_snp}")
df_str_human = pd.read_csv(file_str_snp, sep='\t', header=None)
df_str_all = pd.read_csv(file_allstr_snp, sep='\t', header=None)
df_str_human.shape[0]/1400 # 0
df_str_all.shape[0]/28678 # 6.97e-05

# +- 1kb
file_str_snp_1k = os.path.join(path_root, 'Human_STR_SNP_1k.pos')
file_allstr_snp_1k = os.path.join(path_root, 'All_STR_SNP_1k.pos')
# expand STR - SNP
os.system(f"bedtools intersect -a {file_str_human} -b {file_snp_bed_sort_1k} -wo "
          f"> {file_str_snp_1k}")
# all STR - SNP
os.system(f"bedtools intersect -a {file_str_all} -b {file_snp_bed_sort_1k} -wo "
          f"> {file_allstr_snp_1k}")
df_str_human_1k = pd.read_csv(file_str_snp_1k, sep='\t', header=None)
df_str_all_1k = pd.read_csv(file_allstr_snp_1k, sep='\t', header=None)
df_str_human_1k.shape[0]/1400 # 0.0035
df_str_all_1k.shape[0]/28678 # 0.0051


# genes owning STRs
file_gene = os.path.join(path_root, 'Gene_hg38/genes.protein.gencode.v38.bed')
file_gene_sort = os.path.join(path_root, 'Gene_hg38/genes.protein.sort.bed')
os.system(f"bedtools sort -i {file_gene} > {file_gene_sort}")

# STRGENE
file_str_human = os.path.join(path_root, 'human_expan.uniq.bed')
file_str_all = os.path.join(path_root, 'All_ortholog_STR.txt')
file_gene_str = os.path.join(path_root, 'genes_str.bed')
os.system(f"bedtools intersect -a {file_gene_sort} -b {file_str_human} "
          f"> {file_gene_str}")

# calculate percentage
file_strgene_snp = os.path.join(path_root, 'Human_STRgene_SNP.pos')
file_allgene_snp = os.path.join(path_root, 'ALL_gene_SNP.pos')
# expand STR - SNP
os.system(f"bedtools intersect -a {file_gene_str} -b {file_snp_bed_sort} -wo "
          f"> {file_strgene_snp}")
# all STR - SNP
os.system(f"bedtools intersect -a {file_gene_sort} -b {file_snp_bed_sort} -wo "
          f"> {file_allgene_snp}")
df_strgene = pd.read_csv(file_strgene_snp, sep='\t', header=None) # 0
df_allgene = pd.read_csv(file_allgene_snp, sep='\t', header=None) # 4350

# DHS owning STR
file_dhs = os.path.join(path_root, 'DHS_hg38.txt')

# strDHS
file_str_human = os.path.join(path_root, 'human_expan.uniq.bed')
file_strdhs = os.path.join(path_root, 'DHSs_str.bed')
os.system(f"bedtools intersect -a {file_dhs} -b {file_str_human} "
          f"> {file_strdhs}")

# calculate percentage
file_strdhs_snp = os.path.join(path_root, 'Human_STRDHS_SNP.pos')
file_alldhs_snp = os.path.join(path_root, 'ALL_DHS_SNP.pos')
# expand STR - SNP
os.system(f"bedtools intersect -a {file_strdhs} -b {file_snp_bed_sort} -wo "
          f"> {file_strdhs_snp}")
# all STR - SNP
os.system(f"bedtools intersect -a {file_dhs} -b {file_snp_bed_sort} -wo "
          f"> {file_alldhs_snp}")
file_strdhs_snp = pd.read_csv(file_strdhs_snp, sep='\t', header=None) # 0
file_alldhs_snp = pd.read_csv(file_alldhs_snp, sep='\t', header=None) # 4350


if __name__ == '__main__':
    time_start = time()

    time_end = time()
    print(time_end - time_start)
