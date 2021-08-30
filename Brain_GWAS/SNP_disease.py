#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: SNP_disease.py
# @time: 2021/8/12 16:24

import os
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from multiprocessing import Pool
from functools import partial
from subprocess import check_output

path_DisGeNET = '/mdshare/node9/zy/Reg_brain/DisGeNET'

# select disease
file_disease = os.path.join(path_DisGeNET, 'disease_mappings_to_attributes.tsv')
df_disease = pd.read_csv(file_disease, sep='\t', index_col=0)
sel_class = {'C16', 'C10'}
df_no_na = df_disease.loc[
    df_disease['diseaseClassMSH'].apply(lambda x: isinstance(x, str)), ]
df_brain = df_no_na.loc[
    df_no_na['diseaseClassMSH'].apply(
        lambda x: len(set(x.split(';')).intersection(sel_class)) == 2), ]
sel_id = df_brain.index
sel_disease = df_brain.loc[sel_id, 'name']
file_sel_disease = os.path.join(path_DisGeNET, 'disease_sel.tsv')
sel_disease.to_csv(file_sel_disease, sep='\t')

# SNP: build fisher tables
path_disease_SNP = '/mdshare/node9/zy/Reg_brain/disease_SNP/'
path_mid = '/mdshare/node9/zy/Reg_brain/SNP_bed'
path_dhs = '/mdshare/node9/zy/Reg_brain/aceDHS'
path_SNP_out = '/mdshare/node9/zy/Reg_brain/SNP_aceDHS'
file_acedhs = os.path.join(path_dhs, 'DHS_gene.bed')
file_alldhs = os.path.join(path_dhs, 'DHS.pos')
num_dhs = 1867665
num_ace_dhs = 3538


def build_table(path_disease_SNP, path_mid, disease):
    file_snp_bed = os.path.join(path_disease_SNP, f"{disease}.bed")
    sub_path = os.path.join(path_mid, disease)
    if not os.path.exists(sub_path):
        os.mkdir(sub_path)
    file_snp_bed_sort = os.path.join(sub_path, 'Brain_SNP_hg19.sort.bed')
    os.system(f"bedtools sort -i {file_snp_bed} > {file_snp_bed_sort}")
    file_acedhs_snp = os.path.join(sub_path, 'acedhs_SNP.pos')
    file_alldhs_snp = os.path.join(sub_path, 'allDHS_SNP.pos')
    # acedhs - SNP
    os.system(f"bedtools intersect -a {file_acedhs} -b {file_snp_bed_sort} -wo "
              f"> {file_acedhs_snp}")
    # allDHS - SNP
    os.system(f"bedtools intersect -a {file_alldhs} -b {file_snp_bed_sort} -wo "
              f"> {file_alldhs_snp}")
    num_snp_acedhs = int(str(check_output(f"wc -l {file_acedhs_snp}",
                                 shell=True).strip()).split(' ')[0][2:])
    if num_snp_acedhs == 0:
        return
    else:
        df_acedhs = pd.read_csv(file_acedhs_snp, sep='\t', header=None)
        num_snp_acedhs = len(df_acedhs.iloc[:, 3].unique())
        df_alldhs = pd.read_csv(file_alldhs_snp, sep='\t', header=None)
        num_snp_dhs = len(df_alldhs.iloc[:, 3].unique())
        table = np.array([[num_snp_acedhs, num_snp_dhs - num_snp_acedhs],
                          [num_ace_dhs - num_snp_acedhs,
                           num_dhs - num_snp_dhs - num_ace_dhs + num_snp_acedhs]])
        oddsr, p = fisher_exact(table, alternative='two-sided')
        print({'disease': disease, 'pvalue': p})

        return {'disease': disease, 'pvalue': p}


pool = Pool(30)
func_calc = partial(build_table, path_disease_SNP, path_mid)
result = pool.map(func_calc, sel_disease.index)
pool.close()
result_dicts = [one for one in result if one is not None]
df_res = pd.DataFrame(result_dicts)
df_res['disease_name'] = sel_disease.loc[df_res['disease']].to_list()
df_res.to_csv(os.path.join(path_SNP_out, 'results_snp.txt'),
              sep='\t', index=False)

# SNP +- 1kb
path_disease_SNP_1kb = '/mdshare/node9/zy/Reg_brain/disease_SNP_1kb/'
path_mid_1kb = '/mdshare/node9/zy/Reg_brain/SNP_bed_1kb'
pool = Pool(30)
func_calc = partial(build_table, path_disease_SNP_1kb, path_mid_1kb)
result = pool.map(func_calc, sel_disease.index)
pool.close()
result_dicts = [one for one in result if one is not None]
df_res = pd.DataFrame(result_dicts)
df_res['disease_name'] = sel_disease.loc[df_res['disease']].to_list()
df_res.to_csv(os.path.join(path_SNP_out, 'results_snp_1kb.txt'),
              sep='\t', index=False)


