#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: gene_disease.py
# @time: 2021/8/17 15:38

import os
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from multiprocessing import Pool
from functools import partial

# merge Hi-C data


# nervous diseases
path_DisGeNET = '/mdshare/node9/zy/Reg_brain/DisGeNET'
file_sel_disease = os.path.join(path_DisGeNET, 'disease_sel.tsv')
sel_disease = pd.read_csv(file_sel_disease, sep='\t', index_col=0)

# aceDHS-genes
path_hic = '/mdshare/node9/zy/Reg_brain/Hi-C'
file_hic = os.path.join(path_hic, 'Promoter-anchored_chromatin_loops.bed')
df_hic = pd.read_csv(file_hic, sep='\t')
df_hic = df_hic.loc[:, ['chromosome', 'TSS_start', 'TSS_end',
                        'Interacting_region_start', 'Interacting_region_end']]
df_hic['chromosome'] = df_hic['chromosome'].apply(lambda x: f"chr{x}")

# gene name
file_EP = os.path.join(path_hic, 'INT-16_HiC_EP_linkages_cross_assembly.csv')
df_EP = pd.read_csv(file_EP)
df_gene = df_EP.loc[:, ['Enhancer_Chromosome_hg19',
                        'Transcription_Start_Site_hg19',
                        'Target_Gene_Name', 'Target_Ensembl_Name']]
df_gene.columns = ['chromosome', 'TSS_start',
                   'Target_Gene_Name', 'Target_Ensembl_Name']

# merge EP
df_region_gene_merge = pd.merge(df_hic, df_gene, on=['chromosome', 'TSS_start'])
df_region_gene = \
    df_region_gene_merge.loc[:, ['chromosome', 'Interacting_region_start',
                                 'Interacting_region_end', 'Target_Gene_Name',
                                 'Target_Ensembl_Name']]
df_region_gene = df_region_gene.dropna()
df_region_gene = df_region_gene.drop_duplicates()
file_EG_bed = os.path.join(path_hic, 'EG_hg19.bed')
df_region_gene.to_csv(file_EG_bed, sep='\t', index=False, header=None)
file_EG_bed_sort = os.path.join(path_hic, 'EG_hg19.sort.bed')
os.system(f"bedtools sort -i {file_EG_bed} > {file_EG_bed_sort}")

# overlap sceDHS gene
path_dhs = '/mdshare/node9/zy/Reg_brain/aceDHS'
file_acedhs = os.path.join(path_dhs, 'DHS_gene.bed')
path_gene = '/mdshare/node9/zy/Reg_brain/gene_bed'
file_aceDHS_gene = os.path.join(path_gene, 'aceDHS_gene.bed')
os.system(f"bedtools intersect -a {file_acedhs} -b {file_EG_bed_sort} -wo "
          f"> {file_aceDHS_gene}")
df_aceDHS_gene = pd.read_csv(file_aceDHS_gene, sep='\t', header=None)
set_aceDHS_gene_distal = set(df_aceDHS_gene.iloc[:, 9])
num_aceDHS_gene_distal = len(set_aceDHS_gene_distal)

# disease-genes
file_disease_gene = os.path.join(path_DisGeNET, 'curated_gene_disease_associations.tsv')
df_disease_gene = pd.read_csv(file_disease_gene, sep='\t')
file_all_gene = os.path.join(path_gene, 'genes.protein.gencode.v19.bed')
df_all_gene = pd.read_csv(file_all_gene, sep='\t', header=None)
all_disease_gene = set(df_disease_gene['geneSymbol'].to_list())
all_gene = set(df_all_gene.iloc[:, 3].to_list())
len(all_gene.intersection(all_disease_gene))
set_aceDHS_gene_distal = all_gene.intersection(set_aceDHS_gene_distal)
num_ace = len(set_aceDHS_gene_distal)  # 597
num_all_gene = len(all_gene)


# function of fisher test
def fisher_test(df_disease_gene, disease):
    disease_gene = df_disease_gene.loc[
        df_disease_gene['diseaseId'] == disease, 'geneSymbol']
    disease_gene = set(disease_gene.to_list()).intersection(all_gene)
    overlap_disease_ace = set(disease_gene).intersection(set_aceDHS_gene_distal)
    num_disease = len(disease_gene)
    num_overlap = len(overlap_disease_ace)
    table = np.array([[num_overlap, num_ace - num_overlap],
                      [num_disease - num_overlap,
                       num_all_gene - num_ace - num_disease + num_overlap]])
    oddsr, p = fisher_exact(table, alternative='two-sided')
    # print({'disease': disease, 'pvalue': p})

    return {'disease': disease, 'pvalue': p}


all_gene_diseases = set(df_disease_gene['diseaseId'].to_list())
input_genes = set(sel_disease.index).intersection(all_gene_diseases)
pool = Pool(30)
func_calc = partial(fisher_test, df_disease_gene)
result = pool.map(func_calc, input_genes)
pool.close()
df_res = pd.DataFrame(result)
df_res['disease_name'] = sel_disease.loc[df_res['disease'], 'name'].to_list()
df_res = df_res.sort_values(by='pvalue')
df_res.to_csv(os.path.join(path_gene, 'results_disease_gene.txt'),
              sep='\t', index=False)
