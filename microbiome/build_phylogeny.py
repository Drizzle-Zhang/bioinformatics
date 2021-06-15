 #!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: build_phylogeny.py
# @time: 8/8/20 12:53 PM

from time import time
import os
import numpy as np
import pandas as pd


def build_phylogeny():
    file_otus = '/home/drizzle_zhang/microbiome/result/2.OTUs/OTUs_stat/OTUs_tax_even.csv'
    file_seq = '/home/drizzle_zhang/microbiome/result/2.OTUs/OTUs_stat/rep_seqs_tax.csv'
    num_top = 500
    file_fa = f"/home/drizzle_zhang/microbiome/result/2.OTUs/rep_seqs_top{num_top}.fasta"
    file_tree = f"/home/drizzle_zhang/microbiome/result/2.OTUs/rep_phylo_top{num_top}.tre"

    df_otus = pd.read_csv(file_otus, sep=',', index_col=0)
    df_otus = df_otus.iloc[:, :-1]
    df_otus_top = \
        df_otus.iloc[np.argsort(-np.mean(df_otus, axis=1), )[:num_top], :]

    df_seq = pd.read_csv(file_seq, sep=',', index_col=0)
    df_seq_top = df_seq.loc[df_otus_top.index, :]
    df_seq_top['OTU_ID'] = df_seq_top.index
    with open(file_fa, 'w') as w_fa:
        for sub_dict in df_seq_top.to_dict('records'):
            line_description = \
                f">{sub_dict['OTU_ID']}|{sub_dict['OTU_taxonomy']}\n"
            line_seq = f"{sub_dict['represent_seq']}\n"
            w_fa.write(line_description)
            w_fa.write(line_seq)

    fasttree = '/home/drizzle_zhang/tools/FastTree'
    os.system(f"{fasttree} -gtr -nt {file_fa} > {file_tree}")

    return


if __name__ == '__main__':
    time_start = time()

    time_end = time()
    print(time_end - time_start)
