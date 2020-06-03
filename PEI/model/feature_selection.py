#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: feature_selection.py
# @time: 6/1/20 6:15 PM

from time import time
import pandas as pd
import numpy as np
import os


def combine_corr_label(file_corr, file_label, file_out):
    df_corr = pd.read_csv(file_corr, sep='\t')
    df_label = pd.read_csv(file_label, sep='\t')
    df_label['label'] = np.full(df_label.shape[0], 1)

    df_combine = pd.merge(
        df_corr, df_label, how='left',
        on=['gene', 'dhs_id', 'type_cre', 'ref_dhs_id'])
    df_combine = df_combine.fillna(0)
    df_combine.to_csv(file_out, sep='\t', index=None)

    return


if __name__ == '__main__':
    time_start = time()
    # path_root = '/local/zy/PEI'
    path_root = '/lustre/tianlab/zhangyu/PEI'
    path_origin = path_root + '/origin_data'
    # path_mid = path_root + '/mid_data'
    path_mid = path_root + '/mid_data_correct'

    file_corr_GM12878 = \
        path_mid + '/cell_line/model_input/GM12878/correlation.txt'
    file_label_GM12878 = \
        path_mid + '/training_label/label_interactions/GM12878/GM12878.txt'
    file_out_GM12878 = \
        path_mid + '/cell_line/model_input/GM12878/selection.txt'
    combine_corr_label(file_corr_GM12878, file_label_GM12878, file_out_GM12878)

    time_end = time()
    print(time_end - time_start)
