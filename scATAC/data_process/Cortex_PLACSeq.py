# _*_ coding: utf-8 _*_
# @author: Drizzle_Zhang
# @file: Cortex_PLACSeq.py
# @time: 2021/11/30 11:01

from time import time
import os
import pandas as pd
from collections import defaultdict
from shutil import copyfile


file_gene_hg19 = '/root/scATAC/Gene_anno/Gene_hg19/promoters.up2k.protein.gencode.v19.merge.bed'
path_plac = '/root/scATAC/pcHi-C/Cortex_PLACSeq'

list_cell = ['Microglia', 'Neuron', 'Oligo']
for cell in list_cell:
    file_ori = os.path.join(path_plac, f"{cell}.txt")
    path_cell = os.path.join(path_plac, f"hg19/{cell}/")
    if not os.path.exists(path_cell):
        os.mkdir(path_cell)
    path_process = os.path.join(path_cell, 'processed_files')
    if not os.path.exists(path_process):
        os.mkdir(path_process)
    file_1 = os.path.join(path_process, 'interactome1.bed')
    file_2 = os.path.join(path_process, 'interactome2.bed')
    with open(file_ori, 'r') as r_ori:
        w_1 = open(file_1, 'w')
        fmt_1 = "{chrom1}\t{start1}\t{end1}\t{interactome_id}\n"
        w_2 = open(file_2, 'w')
        fmt_2 = "{chrom2}\t{start2}\t{end2}\t{interactome_id}\n"
        i = 0
        for line in r_ori:
            i = i + 1
            list_line = line.strip().split('\t')
            w_1.write(fmt_1.format(**dict(chrom1=list_line[0], start1=list_line[1],
                                          end1=list_line[2], interactome_id=f"interactome_{i}")))
            w_2.write(fmt_2.format(**dict(chrom2=list_line[3], start2=list_line[4],
                                          end2=list_line[5], interactome_id=f"interactome_{i}")))
        w_1.close()
        w_2.close()
    file_gene_1 = os.path.join(path_process, 'gene_interactome1.bed')
    file_gene_2 = os.path.join(path_process, 'gene_interactome2.bed')
    os.system(f"bedtools intersect -a {file_1} -b {file_gene_hg19} -wao > {file_gene_1}")
    os.system(f"bedtools intersect -a {file_2} -b {file_gene_hg19} -wao > {file_gene_2}")
    df_gene_1 = pd.read_csv(file_gene_1, sep='\t', header=None)
    df_gene_1 = df_gene_1.loc[:, [0, 1, 2, 3, 7]]
    df_gene_1.columns = ['chrom1', 'start1', 'end1', 'interactome_id', 'gene1']
    df_gene_2 = pd.read_csv(file_gene_2, sep='\t', header=None)
    df_gene_2 = df_gene_2.loc[:, [0, 1, 2, 3, 7]]
    df_gene_2.columns = ['chrom2', 'start2', 'end2', 'interactome_id', 'gene2']
    df_gene = pd.merge(df_gene_1, df_gene_2, on='interactome_id')
    file_gene_12 = os.path.join(path_process, 'gene_interactome12.bed')
    df_gene.to_csv(file_gene_12, sep='\t', header=False, index=False)
    file_pp_pre = os.path.join(path_cell, 'PP_pre.txt')
    file_po_pre = os.path.join(path_cell, 'PO_pre.txt')
    with open(file_gene_12, 'r') as r_gene:
        w_pp = open(file_pp_pre, 'w')
        w_po = open(file_po_pre, 'w')
        for line in r_gene:
            list_line = line.strip().split('\t')
            gene1 = list_line[4]
            gene2 = list_line[8]
            if gene1 != '.' and gene2 != '.':
                for sub_gene1 in gene1.split(','):
                    sub_gene1 = sub_gene1.split('<-')[0]
                    for sub_gene2 in gene2.split(','):
                        sub_gene2 = sub_gene2.split('<-')[0]
                        w_pp.write(f"{sub_gene1}\t{sub_gene2}\n")
            elif gene1 != '.' and gene2 == '.':
                for sub_gene1 in gene1.split(','):
                    sub_gene1 = sub_gene1.split('<-')[0]
                    w_po.write(f"{list_line[5]}\t{list_line[6]}\t{list_line[7]}\t"
                               f"{list_line[5]}:{list_line[6]}-{list_line[7]}\t{sub_gene1}\n")
            elif gene1 == '.' and gene2 != '.':
                for sub_gene2 in gene2.split(','):
                    sub_gene2 = sub_gene2.split('<-')[0]
                    w_po.write(f"{list_line[0]}\t{list_line[1]}\t{list_line[2]}\t"
                               f"{list_line[0]}:{list_line[1]}-{list_line[2]}\t{sub_gene2}\n")
        w_pp.close()
        w_po.close()
    file_pp = os.path.join(path_cell, 'PP.txt')
    file_po = os.path.join(path_cell, 'PO.txt')
    os.system(f"sort {file_pp_pre} | uniq > {file_pp}")
    os.system(f"sort {file_po_pre} | uniq > {file_po}")

# merge all cell types
pp_files = ' '.join([os.path.join(path_plac, f"hg19/{cell}/PP.txt") for cell in list_cell])
file_pp_merge = os.path.join(path_plac, f"hg19/PP.txt")
os.system(f"cat {pp_files} | sort | uniq > {file_pp_merge}")
po_files = ' '.join([os.path.join(path_plac, f"hg19/{cell}/PO.txt") for cell in list_cell])
file_po_merge = os.path.join(path_plac, f"hg19/PO.txt")
os.system(f"cat {po_files} | sort | uniq > {file_po_merge}")

# hg19 to hg38
file_chain = '/root/tools/files_liftOver/hg19ToHg38.over.chain.gz'
liftover = '/root/tools/liftOver'

path_hg38 = '/root/scATAC/pcHi-C/Cortex_PLACSeq/hg38'
path_process = os.path.join(path_hg38, 'processed_files')
if not os.path.exists(path_process):
    os.mkdir(path_process)
file_hg19 = os.path.join(path_process, 'hg19.bed')
df_hg19 = pd.read_csv(file_po_merge, sep='\t', header=None)
df_hg19.loc[:, 'interaction_id'] = df_hg19.apply(lambda x: f"{x.iloc[3]}_{x.iloc[4]}", axis=1)
df_hg19.to_csv(file_hg19, sep='\t', index=None, header=None)
file_hg38 = os.path.join(path_process, 'hg38.bed')
file_prefix = file_hg19 + '.prefix'
file_suffix = file_hg19 + '.suffix'
file_hg38_prefix = file_hg38 + '.prefix'
file_hg38_format = file_hg38 + '.format'
file_ummap = os.path.join(path_process, 'unmap.bed')
os.system(f"cut -f 1,2,3,6 {file_hg19} > {file_prefix}")
os.system(f"cut -f 4,5,6 {file_hg19} > {file_suffix}")
os.system(f"{liftover} {file_prefix} {file_chain} "
          f"{file_hg38_prefix} {file_ummap}")
dict_peak_score = defaultdict(list)
with open(file_suffix, 'r') as r_f:
    for line in r_f:
        list_line = line.strip().split('\t')
        dict_peak_score[list_line[2]].append(list_line[0:2])
with open(file_hg38_format, 'w') as w_f:
    fmt = "{chrom}\t{start}\t{end}\t{interaction_id}\t{peak_id}\t{gene}\n"
    with open(file_hg38_prefix, 'r') as r_hg38:
        for line in r_hg38:
            list_line = line.strip().split('\t')
            list_suffix = dict_peak_score[list_line[3]][0]
            dict_hg38 = dict(
                chrom=list_line[0], start=list_line[1], end=list_line[2],
                interaction_id=list_line[3],
                peak_id=f"{list_line[0]}:{list_line[1]}-{list_line[2]}",
                gene=list_suffix[1]
            )
            w_f.write(fmt.format(**dict_hg38))

df_old = pd.read_csv(file_prefix, sep='\t', header=None)
length_old = df_old.iloc[:, 2] - df_old.iloc[:, 1]
df_old['length'] = length_old
df_bed = pd.read_csv(file_hg38_format, sep='\t', header=None)
length = df_bed.iloc[:, 2] - df_bed.iloc[:, 1]
# df_bed['length'] = length
df_bed = df_bed.loc[length < 20000, :]
df_bed = df_bed.drop_duplicates()
df_bed = df_bed.iloc[:, [0, 1, 2, 4, 5]]
df_bed.to_csv(file_hg38, sep='\t', index=None, header=None)
copyfile(file_pp_merge, os.path.join(path_hg38, 'PP.txt'))
copyfile(file_hg38, os.path.join(path_hg38, 'PO.txt'))

path_hg19 = '/root/scATAC/pcHi-C/Cortex_PLACSeq/hg19'
path_hg38 = '/root/scATAC/pcHi-C/Cortex_PLACSeq/hg38'
list_cell = ['Microglia', 'Neuron', 'Oligo']
for cell in list_cell:
    path_cell_hg19 = os.path.join(path_hg19, cell)
    path_cell_hg38 = os.path.join(path_hg38, cell)
    if not os.path.exists(path_cell_hg38):
        os.mkdir(path_cell_hg38)
    path_process = os.path.join(path_cell_hg38, 'processed_files')
    if not os.path.exists(path_process):
        os.mkdir(path_process)
    file_po_hg19 = os.path.join(path_cell_hg19, 'PO.txt')
    file_hg19 = os.path.join(path_process, 'hg19.bed')
    df_hg19 = pd.read_csv(file_po_hg19, sep='\t', header=None)
    df_hg19.loc[:, 'interaction_id'] = df_hg19.apply(lambda x: f"{x.iloc[3]}_{x.iloc[4]}", axis=1)
    df_hg19.to_csv(file_hg19, sep='\t', index=None, header=None)
    file_hg38 = os.path.join(path_process, 'hg38.bed')
    file_prefix = file_hg19 + '.prefix'
    file_suffix = file_hg19 + '.suffix'
    file_hg38_prefix = file_hg38 + '.prefix'
    file_hg38_format = file_hg38 + '.format'
    file_ummap = os.path.join(path_process, 'unmap.bed')
    os.system(f"cut -f 1,2,3,6 {file_hg19} > {file_prefix}")
    os.system(f"cut -f 4,5,6 {file_hg19} > {file_suffix}")
    os.system(f"{liftover} {file_prefix} {file_chain} "
              f"{file_hg38_prefix} {file_ummap}")
    dict_peak_score = defaultdict(list)
    with open(file_suffix, 'r') as r_f:
        for line in r_f:
            list_line = line.strip().split('\t')
            dict_peak_score[list_line[2]].append(list_line[0:2])
    with open(file_hg38_format, 'w') as w_f:
        fmt = "{chrom}\t{start}\t{end}\t{interaction_id}\t{peak_id}\t{gene}\n"
        with open(file_hg38_prefix, 'r') as r_hg38:
            for line in r_hg38:
                list_line = line.strip().split('\t')
                list_suffix = dict_peak_score[list_line[3]][0]
                dict_hg38 = dict(
                    chrom=list_line[0], start=list_line[1], end=list_line[2],
                    interaction_id=list_line[3],
                    peak_id=f"{list_line[0]}:{list_line[1]}-{list_line[2]}",
                    gene=list_suffix[1]
                )
                w_f.write(fmt.format(**dict_hg38))

    df_old = pd.read_csv(file_prefix, sep='\t', header=None)
    length_old = df_old.iloc[:, 2] - df_old.iloc[:, 1]
    df_old['length'] = length_old
    df_bed = pd.read_csv(file_hg38_format, sep='\t', header=None)
    length = df_bed.iloc[:, 2] - df_bed.iloc[:, 1]
    # df_bed['length'] = length
    df_bed = df_bed.loc[length < 20000, :]
    df_bed = df_bed.drop_duplicates()
    df_bed = df_bed.iloc[:, [0, 1, 2, 4, 5]]
    df_bed.to_csv(file_hg38, sep='\t', index=None, header=None)
    file_pp_hg19 = os.path.join(path_cell_hg19, 'PP.txt')
    copyfile(file_pp_hg19, os.path.join(path_cell_hg38, 'PP.txt'))
    copyfile(file_hg38, os.path.join(path_cell_hg38, 'PO.txt'))


if __name__ == '__main__':
    time_start = time()

    time_end = time()
    print(time_end - time_start)
