#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: Campbell.py
# @time: 12/15/20 10:04 AM

from time import time
import pandas as pd

file_ori = '/home/zy/scRef/sc_data/Campbell/Campbell_exp_sc_mat_cluster_original.txt'
file_mod = '/home/zy/scRef/sc_data/Campbell_exp_sc_mat_cluster_original.txt'

df_label = pd.read_csv(file_ori, sep='\t')
df_label.columns = ['Cell_ID', 'Cluster']

df_label = df_label.replace({'Cluster': 'Astrocyte'}, value='Astrocytes')
df_label = df_label.replace({'Cluster': 'EndothelialCells'}, value='Endothelial cells')
df_label = df_label.replace({'Cluster': 'Fibroblast'}, value='VLMCs')
df_label = df_label.replace({'Cluster': 'microglia'}, value='PVMs & Microglia')
df_label = df_label.replace({'Cluster': 'MuralCells'}, value='Mural cells')
df_label = df_label.replace({'Cluster': 'Oligodend'}, value='Oligodendrocytes')
df_label = df_label.replace({'Cluster': 'ParsTuber'}, value='Pars tuberalis')
df_label = df_label.replace({'Cluster': 'Tanycyte'}, value='Tanycytes')

df_label.to_csv(file_mod, sep='\t', index=None)

if __name__ == '__main__':
    time_start = time()

    time_end = time()
    print(time_end - time_start)
