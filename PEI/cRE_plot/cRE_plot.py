#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: cRE_plot.py
# @time: 11/6/19 10:46 PM

from time import time
from annotate_cRE.annotate_DHS_mu02 import generate_file_list
import os
from multiprocessing import Pool
from functools import partial
import pandas as pd
from subprocess import check_output


def sub_count(dict_in):
    path_in = dict_in['path_in']
    file = dict_in['file']
    if file[-4:] != '.bed':
        return
    file_in = os.path.join(path_in, file)
    total = int(str(check_output(f"wc -l {file_in}",
                                 shell=True).strip()).split(' ')[0][2:])
    promoter = int(check_output(f"grep -w Promoter {file_in} | wc -l ",
                                shell=True).strip())
    enhancer = int(check_output(f"grep -w Enhancer {file_in} | wc -l ",
                                shell=True).strip())
    other = int(check_output(f"grep -w Other {file_in} | wc -l ",
                             shell=True).strip())
    dict_in['total'] = total
    dict_in['promoter'] = promoter
    dict_in['enhancer'] = enhancer
    dict_in['other'] = other

    return dict_in


def count_cre(path_in, path_out, num_process):
    list_cre = generate_file_list(path_in, path_out)
    os.system(f"rm -rf {path_out}/*")
    pool = Pool(processes=num_process)
    list_out = pool.map(sub_count, list_cre)
    pool.close()

    list_out = [one_dict for one_dict in list_out if one_dict is not None]
    df_count = pd.DataFrame(list_out,
                            columns=['organ', 'life_stage', 'term', 'total',
                                     'promoter', 'enhancer', 'other'])
    df_count.to_csv(
        os.path.join(path_out, 'cRE_count.txt'), sep='\t', index=None)

    df_count_organ = df_count.loc[
        (df_count['life_stage'] == '.') & (df_count['term'] == '.'),
        ['organ', 'total', 'promoter', 'enhancer', 'other']]
    df_count_organ.index = range(df_count_organ.shape[0])
    list_plot_organ = []
    for i in range(df_count_organ.shape[0]):
        organ = df_count_organ.loc[i, 'organ']
        if organ == '.':
            organ = 'All'
        list_plot_organ.append(
            {"Organ": organ, "Count": df_count_organ.loc[i, 'promoter'],
             "Type": 'Promoter'}
        )
        list_plot_organ.append(
            {"Organ": organ, "Count": df_count_organ.loc[i, 'enhancer'],
             "Type": 'Enhancer'}
        )
        list_plot_organ.append(
            {"Organ": organ, "Count": df_count_organ.loc[i, 'other'],
             "Type": 'Other'}
        )
    df_plot_organ = pd.DataFrame(list_plot_organ)
    df_plot_organ.to_csv(
        os.path.join(path_out, 'organ_count.txt'), sep='\t', index=None
    )

    return


if __name__ == '__main__':
    time_start = time()

    time_end = time()
    print(time_end - time_start)