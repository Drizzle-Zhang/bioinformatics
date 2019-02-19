#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: cufflinks.py
# @time: 2018/10/24 20:14

from time import time
from argparse import ArgumentParser
import os
import subprocess


def cufflinks(bam_path, cufflinks_out, process=3, thread=10):
    # read samples from results of tophat
    samples = os.listdir(bam_path)
    list_srr = []
    for sample in samples:
        if sample[:3] == 'SRR' and os.path.isdir(os.path.join(bam_path, sample)):
            list_srr.append(sample)

    # run cufflinks
    subprocesses = []
    for i in range(len(list_srr)):
        if i % process == 0:
            for sub in subprocesses:
                sub.wait()
            subprocesses = []

        subprocesses.append(subprocess.Popen("cufflinks -p " + str(thread) +
                                             " -o " + os.path.join(cufflinks_out, list_srr[i]) + ' '
                                             + os.path.join(os.path.join(bam_path, list_srr[i]), 'accepted_hits.bam'),
                                             shell=True))

    return


def cuffmerge(cufflinks_path, cuffmerge_path, gtf_file, fasta_file, thread=20):
    # build assembiles.txt
    samples = os.listdir(cufflinks_path)
    with open(os.path.join(cufflinks_path, 'assemblies.txt'), 'w') as w_asm:
        for sample in samples:
            if sample[:3] == 'SRR' and os.path.isdir(os.path.join(cufflinks_path, sample)):
                w_asm.write(os.path.join(os.path.join(cufflinks_path, sample)) + '/transcripts.gtf' + '\n')

    # run cuffmerge
    os.system("cuffmerge -g " + gtf_file + " -s " + fasta_file + " -o " + cuffmerge_path + " -p " + str(thread) + " " +
              os.path.join(cufflinks_path, 'assemblies.txt'))

    return


def cuffdiff(bam_path, cuffmerge_path, cuffdiff_path, fasta_file, group_name, num_group1, thread=20):
    # divide into groups
    samples = os.listdir(bam_path)
    all_sample = []
    for sample in samples:
        if sample[:3] == 'SRR' and os.path.isdir(os.path.join(bam_path, sample)):
            all_sample.append(os.path.join(os.path.join(bam_path, sample), 'accepted_hits.bam'))

    all_sample.sort()
    group1 = all_sample[:num_group1]
    group2 = all_sample[num_group1:]

    # run cuffdiff
    os.system("cuffdiff -o " + cuffdiff_path + " -p " + str(thread) +
              " -L " + group_name + " -u " + os.path.join(cuffmerge_path, 'merged.gtf') + ' ' +
              ','.join(group1) + ' ' + ','.join(group2))

    return


def main_func():
    # take arguments from the commandline
    parser = ArgumentParser(description='An aggregated python script of differential expression analysis '
                                        'using cufflinks')

    parser.add_argument('--bam_path',
                        help='path of bam-files from trimgalore, each sample has a single folder in the path '
                             '(cufflink, cuffdiff)')
    parser.add_argument('--cufflinks_path',
                        help='the path saving the output of cufflinks function, and needs to be made in advance '
                             '(cufflink, cuffmerge)')
    parser.add_argument('--process', nargs=4, type=int,
                        help='management of multiple process and thread, ther are four parameters and they indicate '
                             'the number of sample processed by cufflinks simultaneously, number of threads of '
                             'cufflinks, number of threads of cuffmerge and number of threads of cuffdiff ')
    parser.add_argument('--cuffmerge_path',
                        help='the path saving the output of cuffmerge function, and needs to be made in advance '
                             '(cuffmerge, cuffdiff)')
    parser.add_argument('--cuffdiff_path',
                        help='the path saving the output of cuffdiff function, and needs to be made in advance '
                             '(cuffdiff)')
    parser.add_argument('--gtf_file', help='path of gtf-file involving genome annotation information (cuffmerge)')
    parser.add_argument('--fasta_file', help='path of fasta-file of reference genome (cuffmerge, cuffdiff)')
    parser.add_argument('--group_name', help='groups in the process of differential expression analysis '
                                             '(cuffdiff)')
    parser.add_argument('--num_group1', type=int, help='number of samples in group1 (cuffdiff)')
    parser.add_argument('--procedure', help="procedures which need to run, input a comma-delimited list, for example, "
                                            "'cufflinks,cuffmerge'", required=True)

    args = parser.parse_args()

    procedures = args.procedure
    procedures = procedures.strip().split(',')

    if args.process:
        list_process = args.process
        # cufflinks
        if args.bam_path and args.cufflinks_path and 'cufflinks' in set(procedures):
            cufflinks(args.bam_path, args.cufflinks_path, list_process[0], list_process[1])

        # cuffmerge
        if args.cufflinks_path and args.cuffmerge_path and args.gtf_file and args.fasta_file \
                and 'cuffmerge' in set(procedures):
            cuffmerge(args.cufflinks_path, args.cuffmerge_path, args.gtf_file, args.fasta_file, list_process[2])

        # cuffdiff
        if args.bam_path and args.cuffmerge_path and args.cuffdiff_path and args.fasta_file and args.group_name \
                and args.num_group1 and 'cuffdiff' in set(procedures):
            cuffdiff(args.bam_path, args.cuffmerge_path, args.cuffdiff_path, args.fasta_file, args.group_name,
                     args.num_group1, list_process[3])

    else:
        # cufflinks
        if args.bam_path and args.cufflinks_path and 'cufflinks' in set(procedures):
            cufflinks(args.bam_path, args.cufflinks_path)

        # cuffmerge
        if args.cufflinks_path and args.cuffmerge_path and args.gtf_file and args.fasta_file \
                and 'cuffmerge' in set(procedures):
            cuffmerge(args.cufflinks_path, args.cuffmerge_path, args.gtf_file, args.fasta_file)

        # cuffdiff
        if args.bam_path and args.cuffmerge_path and args.cuffdiff_path and args.fasta_file and args.group_name \
                and args. num_group1 and 'cuffdiff' in set(procedures):
            cuffdiff(args.bam_path, args.cuffmerge_path, args.cuffdiff_path, args.fasta_file, args.group_name,
                     args. num_group1)

    return


if __name__ == '__main__':
    start = time()
    main_func()
    end = time()
    print(end - start)
