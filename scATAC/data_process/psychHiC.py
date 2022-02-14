# _*_ coding: utf-8 _*_
# @author: Drizzle_Zhang
# @file: psychHiC.py
# @time: 2022/1/11 13:59

import os
import pandas as pd
from collections import defaultdict
from shutil import copyfile


# hg19 to hg38
file_chain = '/root/tools/files_liftOver/hg19ToHg38.over.chain.gz'
liftover = '/root/tools/liftOver'

path_hic = '/root/scATAC/pcHi-C/psychHiC/'
path_hg38 = '/root/scATAC/pcHi-C/psychHiC/hg38'
file_hg19 = os.path.join(path_hic, 'EG_hg19.sort.bed')
file_hg38 = os.path.join(path_hg38, 'hg38.bed')
file_hg38_format = file_hg38 + '.format'
file_ummap = os.path.join(path_hg38, 'unmap.bed')
os.system(f"{liftover} {file_hg19} {file_chain} "
          f"{file_hg38_format} {file_ummap}")
df_bed = pd.read_csv(file_hg38_format, sep='\t', header=None)
length = df_bed.iloc[:, 2] - df_bed.iloc[:, 1]
df_bed = df_bed.loc[length < 20000, :]
df_bed = df_bed.drop_duplicates()
df_bed = df_bed.iloc[:, [0, 1, 2, 3]]
df_bed.to_csv(file_hg38, sep='\t', index=None, header=None)

# PP
file_gene_hg38 = '/root/scATAC/Gene_anno/Gene_hg38/promoters.up2k.protein.gencode.v38.bed'
file_PP_intersect = os.path.join(path_hg38, 'PP.intersect')
os.system(f"bedtools intersect -a {file_hg38} -b {file_gene_hg38} -wao > {file_PP_intersect}")
df_PP_intersect = pd.read_csv(file_PP_intersect, sep='\t', header=None)
df_PP_intersect[11] = df_PP_intersect[7].apply(lambda x: x.split('<-')[0])
df_pp = df_PP_intersect.loc[df_PP_intersect[11] != '.', [3, 11]]
df_pp = df_pp.loc[df_pp.apply(lambda x: x.iloc[0] != x.iloc[1], axis=1), :]
df_pp = df_pp.drop_duplicates()
file_PP = os.path.join(path_hg38, 'PP.txt')
df_pp.to_csv(file_PP, sep='\t', index=False, header=False)

# PO
