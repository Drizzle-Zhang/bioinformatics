#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: Baron.py
# @time: 11/11/20 3:18 PM

from time import time
import pandas as pd
import os

path_sc = '/home/zy/scRef/sc_data/Baron/GSE84133_RAW'
path_mouse = '/home/zy/scRef/sc_data/BaronM'
path_human = '/home/zy/scRef/sc_data/BaronH'
files_mouse = ['GSM2230761_mouse1_umifm_counts.csv',
               'GSM2230762_mouse2_umifm_counts.csv']
files_human = ['GSM2230757_human1_umifm_counts.csv',
               'GSM2230758_human2_umifm_counts.csv',
               'GSM2230759_human3_umifm_counts.csv',
               'GSM2230760_human4_umifm_counts.csv']

# mouse
list_exp = []
list_meta = []
for idx, file in enumerate(files_mouse):
    df_sub = pd.read_csv(os.path.join(path_sc, file), sep=',', index_col=0)
    sub_meta = df_sub['assigned_cluster']
    list_meta.append(sub_meta)
    sub_exp = df_sub.iloc[:, 2:].T
    list_exp.append(sub_exp)

df_exp = pd.concat(list_exp, axis=1, sort=False)
df_meta = pd.concat(list_meta, axis=0)
file_mouse_exp = os.path.join(path_mouse, 'cell_exp.txt')
file_mouse_meta = os.path.join(path_mouse, 'cell_meta.txt')
df_exp.to_csv(file_mouse_exp, sep='\t')
df_meta.to_csv(file_mouse_meta, sep='\t', header=True)

# human
list_exp = []
list_meta = []
for idx, file in enumerate(files_human):
    df_sub = pd.read_csv(os.path.join(path_sc, file), sep=',', index_col=0)
    sub_meta = df_sub['assigned_cluster']
    list_meta.append(sub_meta)
    sub_exp = df_sub.iloc[:, 2:].T
    list_exp.append(sub_exp)

df_exp = pd.concat(list_exp, axis=1, sort=False)
df_meta = pd.concat(list_meta, axis=0)
file_human_exp = os.path.join(path_human, 'cell_exp.txt')
file_human_meta = os.path.join(path_human, 'cell_meta.txt')
df_exp.to_csv(file_human_exp, sep='\t')
df_meta.to_csv(file_human_meta, sep='\t', header=True)
