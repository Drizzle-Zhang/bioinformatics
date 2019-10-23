#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: prepare_bed_file.py
# @time: 10/23/19 10:24 PM

from time import time

if __name__ == '__main__':
    time_start = time()
    # get bed file annotating protein-coding genes
    gtf_file = \
        '/home/zy/driver_mutation/data/ENCODE/gencode.v19.annotation.gtf'
    protein_file = \
        '/home/zy/driver_mutation/data/gene/genes.protein.gencode.v19.bed'
    promoter_file = \
        '/home/zy/driver_mutation/data/gene/' \
        'promoters.up2k.protein.gencode.v19.bed'
    with open(protein_file, 'w') as w_gene:
        with open(promoter_file, 'w') as w_pro:
            fmt_gene = "{chrom}\t{start}\t{end}\t{symbol}\t.\t{strand}\n"
            fmt_promoter = "{chrom}\t{start}\t{end}\t{symbol}\t.\t{strand}\n"
            with open(gtf_file, 'r') as r_gtf:
                for line in r_gtf:
                    if line[0] == '#':
                        continue
                    list_line = line.strip().split('\t')
                    if list_line[2] != 'gene':
                        continue
                    list_attr = list_line[8].strip().split('; ')
                    gene_type = list_attr[2][11:-1]
                    if list_attr[2][-15:-1] != "protein_coding":
                        continue
                    gene_name = list_attr[4][11:-1]
                    strand = list_line[6]
                    dict_gene = dict(chrom=list_line[0], start=list_line[3],
                                     end=list_line[4], symbol=gene_name,
                                     strand=strand)
                    w_gene.write(fmt_gene.format(**dict_gene))
                    if strand == '+':
                        pro_start = str(int(list_line[3]) - 2000)
                        pro_end = list_line[3]
                    elif strand == '-':
                        pro_start = list_line[4]
                        pro_end = str(int(list_line[4]) + 2000)
                    else:
                        print('Error')
                        break
                    dict_promoter = dict(chrom=list_line[0], start=pro_start,
                                         end=pro_end, symbol=gene_name,
                                         strand=strand)
                    w_pro.write(fmt_promoter.format(**dict_promoter))

    time_end = time()
    print(time_end - time_start)
