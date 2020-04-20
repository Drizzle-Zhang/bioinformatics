#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: node10_preparation_interaction.py
# @time: 2020/2/29 19:38

from time import time
import re
import os
import pandas as pd
from multiprocessing import Pool


def transform_thurman(file_in, file_out):
    with open(file_out, 'w') as w_out:
        fmt = "{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t" \
              "{loop_id}\t{score}\n"
        with open(file_in, 'r') as r_po:
            for i, line in enumerate(r_po):
                list_line = line.strip().split('\t')
                score = list_line[7]
                loop_id = 'po' + str(i)
                chrom1 = list_line[0]
                start1 = list_line[1]
                end1 = list_line[2]
                chrom2 = list_line[4]
                start2 = list_line[5]
                end2 = list_line[6]
                if chrom1 == chrom2:
                    w_out.write(fmt.format(**locals()))

    return


def transform_insitu(file_in, file_out):
    with open(file_out, 'w') as w_out:
        fmt = "{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t" \
              "{loop_id}\t{score}\n"
        with open(file_in, 'r') as r_po:
            for i, line in enumerate(r_po):
                list_line = line.strip().split('\t')
                if list_line[0] == 'chr1':
                    continue
                score = list_line[12]
                loop_id = str(i)
                chrom1 = 'chr' + list_line[0]
                start1 = list_line[1]
                end1 = list_line[2]
                chrom2 = 'chr' + list_line[3]
                start2 = list_line[4]
                end2 = list_line[5]
                if chrom1 == chrom2:
                    w_out.write(fmt.format(**locals()))

    return


def transform_ng2019(file_pp, file_po, file_out, cutoff):
    file_tmp = file_out + '.tmp'
    with open(file_tmp, 'w') as w_out:
        fmt = "{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t" \
              "{loop_id}\t{score}\n"
        pattern_chrom = re.compile(r'.+:')
        pattern_start = re.compile(r':.+-')
        pattern_end = re.compile(r'-.+')
        # with open(file_pp, 'r') as r_pp:
        #     for i, line in enumerate(r_pp):
        #         list_line = line.strip().split('\t')
        #         if list_line[0] == 'type':
        #             continue
        #         score = list_line[11]
        #         if float(score) < cutoff:
        #             continue
        #         loop_id = 'pp' + str(i)
        #         chrom1 = pattern_chrom.search(list_line[1]).group()[:-1]
        #         start1 = pattern_start.search(list_line[1]).group()[1:-1]
        #         end1 = pattern_end.search(list_line[1]).group()[1:]
        #         chrom2 = pattern_chrom.search(list_line[2]).group()[:-1]
        #         start2 = pattern_start.search(list_line[2]).group()[1:-1]
        #         end2 = pattern_end.search(list_line[2]).group()[1:]
        #         if chrom1 == chrom2:
        #             w_out.write(fmt.format(**locals()))
        with open(file_po, 'r') as r_po:
            for i, line in enumerate(r_po):
                list_line = line.strip().split('\t')
                if list_line[0] == 'type':
                    continue
                score = list_line[11]
                if float(score) < cutoff:
                    continue
                loop_id = 'po' + str(i)
                chrom1 = pattern_chrom.search(list_line[1]).group()[:-1]
                start1 = pattern_start.search(list_line[1]).group()[1:-1]
                end1 = pattern_end.search(list_line[1]).group()[1:]
                chrom2 = pattern_chrom.search(list_line[2]).group()[:-1]
                start2 = pattern_start.search(list_line[2]).group()[1:-1]
                end2 = pattern_end.search(list_line[2]).group()[1:]
                if chrom1 == chrom2:
                    w_out.write(fmt.format(**locals()))

    os.system(f"sort -k 1,1 -k2,2n {file_tmp} > {file_out}")
    os.remove(file_tmp)

    return


def transform_3div(file_in, file_out, cutoff):
    with open(file_out, 'w') as w_out:
        fmt = "{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t" \
              "{loop_id}\t{score}\n"
        with open(file_in, 'r') as r_po:
            for i, line in enumerate(r_po):
                list_line = line.strip().split('\t')
                if list_line[0] == 'idx':
                    continue
                score = list_line[4]
                if float(score) < cutoff:
                    continue
                loop_id = str(i)
                chrom1 = list_line[15]
                start1 = list_line[16]
                end1 = int(start1) + 5000
                chrom2 = chrom1
                start2 = list_line[17]
                end2 = int(start2) + 5000
                w_out.write(fmt.format(**locals()))

    return


def get_pairs(file_in, file_out):
    count_effect = 0
    pattern_gene = re.compile(r'.+<-')
    file_out_tmp = file_out + '.tmp'
    with open(file_out_tmp, 'w') as w_out:
        fmt = "{gene}\t{dhs_id}\t{cre_type}\n"
        with open(file_in, 'r') as r_f:
            for line in r_f:
                list_line = line.strip().split('\t')
                cres1 = list_line[4].split(',')
                cres2 = list_line[10].split(',')
                set_cres1 = set(list_line[4].split(','))
                set_cres2 = set(list_line[10].split(','))
                if ('.' in set_cres1) | ('.' in set_cres2):
                    continue
                else:
                    count_add = 1
                count_effect = count_effect + count_add
                if count_add == 1:
                    dhs_ids1 = list_line[3].split(',')
                    dhs_ids2 = list_line[9].split(',')
                    promoters1 = set(list_line[5].split(','))
                    promoters2 = set(list_line[11].split(','))
                    if ('Protein-Promoter' in set_cres1) & \
                            ('Enhancer' in set_cres2):
                        genes1 = [pattern_gene.search(val).group()[:-2]
                                  for val in promoters1 if val != '.']
                        for i, dhs_id in enumerate(dhs_ids2):
                            for gene1 in genes1:
                                dict_pro2 = dict(gene=gene1, dhs_id=dhs_id,
                                                 cre_type=cres2[i])
                                w_out.write(fmt.format(**dict_pro2))
                    elif ('Enhancer' in set_cres1) & \
                            ('Protein-Promoter' in set_cres2):
                        genes2 = [pattern_gene.search(val).group()[:-2]
                                  for val in promoters2 if val != '.']
                        for i, dhs_id in enumerate(dhs_ids1):
                            for gene2 in genes2:
                                dict_pro1 = dict(gene=gene2, dhs_id=dhs_id,
                                                 cre_type=cres1[i])
                                w_out.write(fmt.format(**dict_pro1))

    os.system(f"sort {file_out_tmp} | uniq > {file_out}")
    os.system(f"grep -w Enhancer {file_out_tmp} | sort | uniq > {file_out}")
    os.remove(file_out_tmp)

    return


def annotate_hic(dict_in):
    file_cre = dict_in['file_cre']
    file_hic = dict_in['file_hic']
    file_out = dict_in['file_out']
    file_pair = dict_in['file_pair']
    func_tran = dict_in['func_tran']
    file_uniform = dict_in['file_uniform']
    if func_tran is None:
        file_uniform = file_hic
    elif func_tran == 'in-situ Hi-C':
        transform_insitu(file_hic, file_uniform)
    elif func_tran == '3DIV':
        transform_3div(file_hic, file_uniform, cutoff=10)
    elif func_tran == 'pcHi-C_ng2019':
        transform_ng2019(file_hic[0], file_hic[1], file_uniform, cutoff=1.3)

    file_bin1 = file_out + '.bin1'
    file_bin2 = file_out + '.bin2'
    os.system(f"cut -f 1,2,3,7 {file_uniform} > {file_bin1}")
    os.system(f"cut -f 4,5,6,7 {file_uniform} > {file_bin2}")
    file_intersect_bin1 = file_out + '.intersect.bin1'
    file_intersect_bin2 = file_out + '.intersect.bin2'
    os.system(f"bedtools intersect -a {file_bin1} -b {file_cre} -loj | "
              f"cut -f 1,2,3,4,8,9,11 | "
              f"grep -w 'Promoter\\|Enhancer' | "
              f"bedtools sort -i > {file_intersect_bin1}")
    os.system(f"bedtools intersect -a {file_bin2} -b {file_cre} -loj | "
              f"cut -f 1,2,3,4,8,9,11 | "
              f"grep -w 'Promoter\\|Enhancer' | "
              f"bedtools sort -i > {file_intersect_bin2}")

    def _merge_bin(df_in):
        dhs_ids = ','.join(df_in[4].tolist())
        cres = ','.join(df_in[5].tolist())
        promoters = ','.join(df_in[6].tolist())
        dict_out = {0: df_in.iloc[0, 0], 1: df_in.iloc[0, 1],
                    2: df_in.iloc[0, 2], 3: df_in.iloc[0, 3],
                    4: dhs_ids, 5: cres, 6: promoters}
        df_dict_out = pd.Series(dict_out)

        return df_dict_out

    df_intersect_bin1 = pd.read_csv(file_intersect_bin1, sep='\t', header=None)
    df_intersect_bin1 = df_intersect_bin1.rename(columns={3: 'key'})
    df_merge_bin1 = df_intersect_bin1.groupby('key').apply(_merge_bin)
    df_intersect_bin2 = pd.read_csv(file_intersect_bin2, sep='\t', header=None)
    df_intersect_bin2 = df_intersect_bin2.rename(columns={3: 'key'})
    df_merge_bin2 = df_intersect_bin2.groupby('key').apply(_merge_bin)
    df_bin1 = pd.read_csv(file_bin1, sep='\t', header=None)
    df_bin2 = pd.read_csv(file_bin2, sep='\t', header=None)
    df_res_bin1 = pd.merge(df_bin1, df_merge_bin1,
                           on=[0, 1, 2, 3], how='outer')
    df_res_bin2 = pd.merge(df_bin2, df_merge_bin2,
                           on=[0, 1, 2, 3], how='outer')
    df_out = pd.merge(df_res_bin1, df_res_bin2, on=3)
    df_out[7] = df_out[3]
    df_out = df_out.drop(3, axis=1)
    df_out.to_csv(file_out, sep='\t', header=None, index=None, na_rep='.')

    os.remove(file_bin1)
    os.remove(file_bin2)
    os.remove(file_intersect_bin1)
    os.remove(file_intersect_bin2)

    get_pairs(file_out, file_pair)

    return


def prepare_interaction_data(path_cre, path_in, path_out, bio_type):
    meta_interaction = pd.read_csv(file_meta, sep='\t')
    meta_interaction = meta_interaction.loc[
                       meta_interaction['Biosample type'] == bio_type, :]
    list_meta = meta_interaction.to_dict('records')

    list_input = []
    for sub_dict in list_meta:
        if bio_type == 'cell line':
            term = sub_dict['Biosample term name']
            str_term = term.replace(' ', '_').replace(
                '/', '+').replace("'", '--')
            chrom_term = sub_dict['Chromatin_inter_name']
            str_chrom_term = chrom_term.replace(' ', '_').replace(
                '/', '+').replace("'", '--')
            path_root_out = os.path.join(path_out, str_chrom_term)
            if os.path.exists(path_root_out):
                os.system(f"rm -rf {path_root_out}")
            os.mkdir(path_root_out)
            path_ref_cre = os.path.join(
                path_cre, str_term + '/cRE_no_exon.txt')

        file_insitu = sub_dict['situ_Hi-C']
        path_insitu = os.path.join(path_in, file_insitu)
        sub_path1 = os.path.join(path_root_out,
                                 '/'.join(file_insitu.split('/')[0:1]))
        sub_path2 = os.path.join(path_root_out,
                                 '/'.join(file_insitu.split('/')[0:2]))
        if not os.path.exists(sub_path1):
            os.mkdir(sub_path1)
        if not os.path.exists(sub_path2):
            os.mkdir(sub_path2)
        path_insitu_uniform = os.path.join(
            sub_path2, 'interactions.uniform.txt')
        path_insitu_out = os.path.join(
            sub_path2, 'interactions.cRE.txt')
        path_insitu_pair = os.path.join(
            sub_path2, 'pairs.gene.cRE.txt')
        list_input.append(
            dict(file_cre=path_ref_cre, file_hic=path_insitu,
                 file_out=path_insitu_out, file_pair=path_insitu_pair,
                 func_tran='in-situ Hi-C', file_uniform=path_insitu_uniform))

        file_3div = sub_dict['3DIV']
        path_3div = os.path.join(path_in, file_3div)
        sub_path1 = os.path.join(path_root_out,
                                 '/'.join(file_3div.split('/')[0:1]))
        sub_path2 = os.path.join(path_root_out,
                                 '/'.join(file_3div.split('/')[0:2]))
        if not os.path.exists(sub_path1):
            os.mkdir(sub_path1)
        if not os.path.exists(sub_path2):
            os.mkdir(sub_path2)
        path_3div_uniform = os.path.join(
            sub_path2, 'interactions.uniform.txt')
        path_3div_out = os.path.join(
            sub_path2, 'interactions.cRE.txt')
        path_3div_pair = os.path.join(
            sub_path2, 'pairs.gene.cRE.txt')
        list_input.append(
            dict(file_cre=path_ref_cre, file_hic=path_3div,
                 file_out=path_3div_out, file_pair=path_3div_pair,
                 func_tran='3DIV', file_uniform=path_3div_uniform))

        file_ng2019 = sub_dict['pcHi-C_ng2019'].split(';')
        path_ng2019 = [os.path.join(path_in, file) for file in file_ng2019]
        sub_path1 = os.path.join(path_root_out,
                                 '/'.join(file_ng2019[0].split('/')[0:1]))
        sub_path2 = os.path.join(path_root_out,
                                 '/'.join(file_ng2019[0].split('/')[0:2]))
        if not os.path.exists(sub_path1):
            os.mkdir(sub_path1)
        if not os.path.exists(sub_path2):
            os.mkdir(sub_path2)
        path_ng2019_uniform = os.path.join(
            sub_path2, 'interactions.uniform.txt')
        path_ng2019_out = os.path.join(
            sub_path2, 'interactions.cRE.txt')
        path_ng2019_pair = os.path.join(
            sub_path2, 'pairs.gene.cRE.txt')
        list_input.append(
            dict(file_cre=path_ref_cre, file_hic=path_ng2019,
                 file_out=path_ng2019_out, file_pair=path_ng2019_pair,
                 func_tran='pcHi-C_ng2019', file_uniform=path_ng2019_uniform))

        file_thurman = sub_dict['Thurman']
        path_thurman = os.path.join(path_in, file_thurman)
        sub_path1 = os.path.join(path_root_out,
                                 '/'.join(file_thurman.split('/')[0:1]))
        sub_path2 = os.path.join(path_root_out,
                                 '/'.join(file_thurman.split('/')[0:2]))
        if not os.path.exists(sub_path1):
            os.mkdir(sub_path1)
        if not os.path.exists(sub_path2):
            os.mkdir(sub_path2)
        path_thurman_uniform = os.path.join(
            sub_path2, 'interactions.uniform.txt')
        path_thurman_out = os.path.join(
            sub_path2, 'interactions.cRE.txt')
        path_thurman_pair = os.path.join(
            sub_path2, 'pairs.gene.cRE.txt')
        list_input.append(
            dict(file_cre=path_ref_cre, file_hic=path_thurman,
                 file_out=path_thurman_out, file_pair=path_thurman_pair,
                 func_tran=None, file_uniform=path_thurman_uniform))

    pool = Pool(processes=40)
    pool.map(annotate_hic, list_input)
    pool.close()

    return


if __name__ == '__main__':
    time_start = time()
    interaction_thurman = \
        '/local/zy/PEI/origin_data/Chromatin_interactions/' \
        'Prediction/Thurman_2012/Thurman_2012.txt'
    interaction_thurman_uniform = \
        '/local/zy/PEI/origin_data/Chromatin_interactions/' \
        'Prediction/Thurman_2012/Thurman_2012.uniform.txt'
    transform_thurman(interaction_thurman, interaction_thurman_uniform)

    # cell line
    file_meta = '/local/zy/PEI/origin_data/meta_file/meta_chrom_inter.txt'
    path_ref_cellline = '/local/zy/PEI/mid_data/cell_line/DHS/cRE_annotation'
    path_interaction_cellline = \
        '/local/zy/PEI/origin_data/Chromatin_interactions'
    path_out_cellline = '/local/zy/PEI/mid_data/cell_line/gene_cre_pairs'
    prepare_interaction_data(
        path_ref_cellline, path_interaction_cellline, path_out_cellline,
        'cell line')

    time_end = time()
    print(time_end - time_start)
