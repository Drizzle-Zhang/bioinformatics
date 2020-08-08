#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: combine_Lefse_results.py
# @time: 7/16/20 5:56 PM

from time import time
import pandas as pd
import os


def combine_res():
    path_files = '/home/drizzle_zhang/microbiome/result/LEfse/plot_0vs3'
    list_df_res = []
    for sub_time in list_time:
        file_res = os.path.join(path_files, f"LEfse_{sub_time}.res")
        df_res = pd.read_csv(file_res, sep='\t', usecols=[0, 2, 3, 4],
                             index_col=0, header=None)
        df_res.columns = [f'Group_{sub_time}',
                          f'LDA_score_{sub_time}', f'P_value_{sub_time}']
        list_df_res.append(df_res)

    df_combine = pd.concat(list_df_res, axis=1, sort=True)
    file_combine = os.path.join(path_files, 'LEfse_combine.res')
    df_combine.to_csv(file_combine, sep='\t')

    df_combine_not_na = df_combine.loc[
        df_combine.apply(lambda x: x.count(), axis=1) != 13, :]
    df_combine_not_na.to_csv(file_combine, sep='\t')
    return


if __name__ == '__main__':
    time_start = time()
    meta_file = '/home/drizzle_zhang/microbiome/result/meta_sample.out.txt'

    df_meta = pd.read_csv(meta_file, sep='\t')
    list_time = df_meta['Time'].unique().tolist()

    time_end = time()
    print(time_end - time_start)
