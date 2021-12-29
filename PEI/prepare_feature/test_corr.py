#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: test_corr.py
# @time: 6/23/20 10:55 AM

from time import time
import pandas as pd
from scipy.stats import spearmanr, pearsonr, kendalltau
import numpy as np


def test(df_mat_pro, df_mat_dhs):
    vec_gene = df_mat_pro.loc['MYO7A', ]
    vec_dhs = df_mat_dhs.loc['DHS<-chr11:76902426-76903417', :]
    a = pd.concat([vec_gene, vec_dhs], axis=1)
    rho, _ = spearmanr(vec_gene, vec_dhs)

    # cel lines
    variables = vec_gene.index
    cells = [cell for cell in variables if len(cell.split('|')) == 1]
    vec_gene = vec_gene.loc[cells]
    vec_dhs = vec_dhs.loc[cells]

    # negative
    vec_dhs = df_mat_dhs.loc['DHS<-chr11:75945830-75948629', cells]
    spearmanr(vec_gene, vec_dhs)

    return


if __name__ == '__main__':
    time_start = time()

    time_end = time()
    print(time_end - time_start)