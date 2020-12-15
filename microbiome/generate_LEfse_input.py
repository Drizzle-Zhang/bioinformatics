#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: generate_LEfse_input.py
# @time: 7/14/20 9:02 AM

from time import time
import pandas as pd
import numpy as np

time_start = time()
# file_mat_otu = \
#     '/home/drizzle_zhang/microbiome/result/2.OTUs/OTUs_stat/OTUs_tax_even.csv'
# df_otu = pd.read_csv(file_mat_otu, sep=',', usecols=['OTU_ID', 'taxonomy'])
file_mat_otu = \
    '/home/drizzle_zhang/microbiome/result/2.OTUs/OTUs_stat/OTUs_tax_even.txt'
df_otu = pd.read_csv(file_mat_otu, sep='\t', usecols=['OTU_ID', 'taxonomy'])

# classification
df_otu['kindom'] = df_otu['taxonomy'].apply(
    lambda x: x.split('; ')[0] if len(x.split('; ')) > 0 else 'Unassigned')
df_otu['phylum'] = df_otu['taxonomy'].apply(
    lambda x: x.split('; ')[1] if len(x.split('; ')) > 1 else 'Unassigned')
df_otu['class'] = df_otu['taxonomy'].apply(
    lambda x: x.split('; ')[2] if len(x.split('; ')) > 2 else 'Unassigned')
df_otu['order'] = df_otu['taxonomy'].apply(
    lambda x: x.split('; ')[3] if len(x.split('; ')) > 3 else 'Unassigned')
df_otu['family'] = df_otu['taxonomy'].apply(
    lambda x: x.split('; ')[4] if len(x.split('; ')) > 4 else 'Unassigned')
df_otu['genus'] = df_otu['taxonomy'].apply(
    lambda x: x.split('; ')[5] if len(x.split('; ')) > 5 else 'Unassigned')
df_otu['species'] = df_otu['taxonomy'].apply(
    lambda x: x.split('; ')[6] if len(x.split('; ')) > 6 else 'Unassigned')

file_taxonomy = \
    '/home/drizzle_zhang/microbiome/result/2.OTUs/OTUs_stat/taxonomy.txt'
df_otu.to_csv(file_taxonomy, sep='\t', index=None)

# OTU matrix
file_mat_otu = \
    '/home/drizzle_zhang/microbiome/result/2.OTUs/OTUs_stat/OTUs_tax_even.csv'
mat_otu = pd.read_csv(file_mat_otu, sep=',', index_col=0)
mat_otu = mat_otu.drop('taxonomy', axis=1)
sum_otus = np.sum(mat_otu)

# generate input file
list_df = []
for kindom in df_otu['kindom'].unique().tolist():
    if kindom == 'Unclassified':
        continue
    df_kindom = df_otu.loc[df_otu['kindom'] == kindom, :]
    mat_kindom = mat_otu.loc[df_kindom['OTU_ID'].tolist(), :]
    abundance_kindom = np.sum(mat_kindom / sum_otus)
    sub_mat = abundance_kindom
    sub_mat['Class'] = kindom
    list_df.append(sub_mat)
    for phylum in df_kindom['phylum'].unique().tolist():
        if phylum == 'Unassigned':
            continue
        df_phylum = df_kindom.loc[df_kindom['phylum'] == phylum, :]
        mat_phylum = mat_otu.loc[df_phylum['OTU_ID'].tolist(), :]
        abundance_phylum = np.sum(mat_phylum / sum_otus)
        sub_mat = abundance_phylum
        sub_mat['Class'] = f"{kindom}|{phylum}"
        list_df.append(sub_mat)
        for sub_class in df_phylum['class'].unique().tolist():
            if sub_class == 'Unassigned':
                continue
            df_class = df_phylum.loc[df_phylum['class'] == sub_class, :]
            mat_class = mat_otu.loc[df_class['OTU_ID'].tolist(), :]
            abundance_class = np.sum(mat_class / sum_otus)
            sub_mat = abundance_class
            sub_mat['Class'] = f"{kindom}|{phylum}|{sub_class}"
            list_df.append(sub_mat)
            for order in df_class['order'].unique().tolist():
                if order == 'Unassigned':
                    continue
                df_order = df_class.loc[df_class['order'] == order, :]
                mat_order = mat_otu.loc[df_order['OTU_ID'].tolist(), :]
                abundance_order = np.sum(mat_order / sum_otus)
                sub_mat = abundance_order
                sub_mat['Class'] = f"{kindom}|{phylum}|{sub_class}|{order}"
                list_df.append(sub_mat)
                for family in df_order['family'].unique().tolist():
                    if family == 'Unassigned':
                        continue
                    df_family = df_order.loc[df_order['family'] == family, :]
                    mat_family = mat_otu.loc[df_family['OTU_ID'].tolist(), :]
                    abundance_family = np.sum(mat_family / sum_otus)
                    sub_mat = abundance_family
                    sub_mat['Class'] = \
                        f"{kindom}|{phylum}|{sub_class}|{order}|{family}"
                    list_df.append(sub_mat)
                    for genus in df_family['genus'].unique().tolist():
                        if genus == 'Unassigned':
                            continue
                        df_genus = \
                            df_family.loc[df_family['genus'] == genus, :]
                        mat_genus = mat_otu.loc[df_genus['OTU_ID'].tolist(), :]
                        abundance_genus = np.sum(mat_genus / sum_otus)
                        sub_mat = abundance_genus
                        sub_mat['Class'] = \
                            f"{kindom}|{phylum}|{sub_class}|{order}|" \
                            f"{family}|{genus}"
                        list_df.append(sub_mat)
                        for species in df_genus['species'].unique().tolist():
                            if species == 'Unassigned':
                                continue
                            df_species = \
                                df_genus.loc[df_genus['genus'] == genus, :]
                            mat_species = mat_otu.loc[
                                          df_species['OTU_ID'].tolist(), :]
                            abundance_species = np.sum(mat_species / sum_otus)
                            sub_mat = abundance_species
                            sub_mat['Class'] = \
                                f"{kindom}|{phylum}|{sub_class}|{order}|" \
                                f"{family}|{genus}|{species}"
                            list_df.append(sub_mat)

df_lefse = pd.concat(list_df, sort=True, axis=1)
df_lefse = df_lefse.T
col_names = ['Class'] + [f'S{num}' for num in range(1, df_lefse.shape[1])]
df_lefse = df_lefse.loc[:, col_names]
file_lefse = '/home/drizzle_zhang/microbiome/result/2.OTUs/OTUs_stat/LEfse.txt'
df_lefse.to_csv(file_lefse, sep='\t', index=None)

time_end = time()
print(time_end - time_start)
