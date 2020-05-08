#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: generate_feature.py
# @time: 2020/5/6 18:29

from time import time
import pandas as pd
import os
import numpy as np
from multiprocessing import Pool, Process
from functools import partial


def generate_pairs():
    # generate P-E pairs for each cell lines or tissues

    return


if __name__ == '__main__':
    time_start = time()
    num_cpu = 40

    path_dhs_cell = \
        '/local/zy/PEI/mid_data/cell_line/DHS/GRCh38tohg19_standard'
    path_dhs_tissue_cluster = \
        '/local/zy/PEI/mid_data/tissue/DHS/GRCh38tohg19_cluster'
    path_dhs_tissue_stan = \
        '/local/zy/PEI/mid_data/tissue/DHS/GRCh38tohg19_standard'

    path_cre_cell = '/local/zy/PEI/mid_data/cell_line/DHS/cRE_annotation'
    path_cre_tissue = '/local/zy/PEI/mid_data/tissue/DHS/cRE_annotation'

    path_feature_cell = '/local/zy/PEI/mid_data/cell_line/model_input'
    path_feature_tissue = '/local/zy/PEI/mid_data/tissue/model_input'

    time_end = time()
    print(time_end - time_start)
