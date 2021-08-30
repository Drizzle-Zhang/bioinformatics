#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: Gene_file.py
# @time: 2021/6/15 22:52

from time import time
import os


def generate_gene_file(gtf_file, protein_file, promoter_file, promoter_merge,
                       exon_file):
    exon_file_tmp = exon_file + '.tmp'
    with open(protein_file, 'w') as w_gene:
        with open(promoter_file, 'w') as w_pro:
            with open(exon_file_tmp, 'w') as w_exon:
                fmt_gene = \
                    "{chrom}\t{start}\t{end}\t{symbol}\t{ensg_id}\t{strand}\n"
                with open(gtf_file, 'r') as r_gtf:
                    for line_gene in r_gtf:
                        if line_gene[0] == '#':
                            continue
                        list_line_gene = line_gene.strip().split('\t')
                        list_attr = list_line_gene[8].strip().split('; ')
                        gene_name = list_attr[2][11:-1]
                        ensg_id = list_attr[0][9:-1]
                        strand = list_line_gene[6]
                        gene_type = list_attr[4][14:-1]
                        if gene_type != "protein_coding":
                            continue
                        if list_line_gene[2] == 'gene':
                            dict_gene = dict(chrom='chr'+list_line_gene[0],
                                             start=list_line_gene[3],
                                             end=list_line_gene[4],
                                             symbol=gene_name,
                                             ensg_id=ensg_id,
                                             strand=strand)
                            w_gene.write(fmt_gene.format(**dict_gene))
                            if strand == '+':
                                pro_start = str(int(list_line_gene[3]) - 2000)
                                pro_end = str(int(list_line_gene[3]) + 2000)
                            elif strand == '-':
                                pro_start = str(int(list_line_gene[4]) - 2000)
                                pro_end = str(int(list_line_gene[4]) + 2000)
                            else:
                                print('Error')
                                break
                            dict_promoter = dict(chrom='chr'+list_line_gene[0],
                                                 start=pro_start,
                                                 end=pro_end,
                                                 symbol=f"{gene_name}<-"
                                                        f"{list_line_gene[0]}:"
                                                        f"{pro_start}-"
                                                        f"{pro_end}",
                                                 ensg_id=ensg_id,
                                                 strand=strand)
                            w_pro.write(fmt_gene.format(**dict_promoter))
                        elif list_line_gene[2] == 'exon':
                            dict_exon = dict(chrom='chr'+list_line_gene[0],
                                             start=list_line_gene[3],
                                             end=list_line_gene[4],
                                             symbol=gene_name,
                                             ensg_id=ensg_id,
                                             strand=strand)
                            w_exon.write(fmt_gene.format(**dict_exon))

    promoter_sort = promoter_file + '.sort'
    os.system(f"bedtools sort -i {promoter_file} > {promoter_sort}")
    os.system(f"bedtools merge -i {promoter_sort} "
              f"-c 4,5,6 -o collapse,collapse,collapse > {promoter_merge}")
    os.system(f"bedtools sort -i {exon_file_tmp} > {exon_file}")
    os.remove(promoter_sort)
    os.remove(exon_file_tmp)

    return


if __name__ == '__main__':
    time_start = time()
    # get bed file annotating protein-coding genes
    path_origin = '/mdshare/node9/zy/Brain_GWAS/Gene_hg38'
    gtf_file_hg38 = os.path.join(path_origin, 'Homo_sapiens.GRCh38.87.gtf')
    protein_file_hg38 = os.path.join(path_origin, 'genes.protein.gencode.v38.bed')
    promoter_file_hg38 = os.path.join(path_origin, 'promoters.up2k.protein.gencode.v38.bed')
    promoter_file_hg38_merge = os.path.join(path_origin, 'promoters.up2k.protein.gencode.v38.merge.bed')
    exon_file_hg38 = os.path.join(path_origin, 'exon.protein.gencode.v38.bed')
    generate_gene_file(gtf_file_hg38, protein_file_hg38, promoter_file_hg38,
                       promoter_file_hg38_merge, exon_file_hg38)

    time_end = time()
    print(time_end - time_start)
