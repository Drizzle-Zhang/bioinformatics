#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: run_LEfse.py
# @time: 7/14/20 4:18 PM

from time import time
import pandas as pd
import os


def run_lefse_one_group(df_meta, df_mat, path_out, sel_dose, sub_time):
    # sub_time = 'B'
    # sel_dose = {0, 3}

    sel_meta = df_meta.loc[
               (df_meta['Time'] == sub_time) &
               (df_meta['Dose'].apply(lambda x: x in sel_dose)), :]
    sel_sample = sel_meta['Sample'].tolist()
    sel_mat = df_mat.loc[:, ['Class'] + sel_sample]

    dict_group = dict(Class='Group')
    for col in sel_mat.columns[1:]:
        dict_group[col] = \
            sel_meta.loc[sel_meta['Sample'] == col, 'Group'].iloc[0]

    mat_out = pd.concat([pd.DataFrame(dict_group, index=[0]), sel_mat])

    file_input = os.path.join(path_out, f"LEfse_{sub_time}.txt")
    mat_out.to_csv(file_input, sep='\t', index=None, header=None)

    # run LEfse
    file_in = os.path.join(path_out, f"LEfse_{sub_time}.in")
    file_res = os.path.join(path_out, f"LEfse_{sub_time}.res")
    path_lefse = '/home/drizzle_zhang/microbiome/galaxy_lefse'
    os.system(f"{os.path.join(path_lefse, 'format_input.py')} {file_input} "
              f"{file_in} -c 1 -o 1000000")
    os.system(
        f"{os.path.join(path_lefse, 'run_lefse.py')} {file_in} {file_res}")

    # plot
    plot_barplot = os.path.join(path_out, f"barplot_{sub_time}.png")
    os.system(f"{os.path.join(path_lefse, 'plot_res.py')} "
              f"{file_res} {plot_barplot} --dpi 600")
    plot_circle = os.path.join(path_out, f"circleplot_{sub_time}.png")
    os.system(f"{os.path.join(path_lefse, 'plot_cladogram.py')} "
              f"{file_res} {plot_circle} --format png --dpi 1000 "
              f"--right_space_prop 0.2")

    return


time_start = time()

file_lefse = '/home/drizzle_zhang/microbiome/result/2.OTUs/OTUs_stat/LEfse.txt'
meta_file = '/home/drizzle_zhang/microbiome/result/meta_sample.out.txt'

df_meta = pd.read_csv(meta_file, sep='\t')
df_mat = pd.read_csv(file_lefse, sep='\t')

sel_dose = {0, 2}
path_out = '/home/drizzle_zhang/microbiome/result/LEfse/plot_0vs2'
list_time = df_meta['Time'].unique().tolist()

for sub_time in list_time:
    run_lefse_one_group(df_meta, df_mat, path_out, sel_dose, sub_time)

time_end = time()
print(time_end - time_start)
