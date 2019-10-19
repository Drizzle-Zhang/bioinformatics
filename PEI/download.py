#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: download.py
# @time: 10/17/19 11:44 AM

import time
import urllib.request
from lxml import etree
import socket
import subprocess
import os
import random


def get_links_roadmap(url, suffix):
    # get download links
    web_html = urllib.request.urlopen(url)
    web_html_str = str(web_html.read())
    html = etree.HTML(web_html_str)
    files = html.xpath('//a/@href')
    files_filter = []
    links = []
    for file in files:
        if file[-len(suffix):] == suffix:
            files_filter.append(file)
            links.append(url + file)

    return files_filter, links


def download_data(files, links, path_out, num_process):
    subprocesses = []
    for i, link in enumerate(links):
        if i % num_process == 0:
            for sub_process in subprocesses:
                sub_process.wait()
            subprocesses = []
        subprocesses.append(
            subprocess.Popen(
                f"wget -O {os.path.join(path_out, files[i])} {link}",
                shell=True))
        random_sleep = random.randint(5, 15)
        time.sleep(random_sleep)
        # if i == 6:
        #     break

    for sub_process in subprocesses:
        sub_process.wait()

    return


def de_gz():
    subprocesses = []
    for i, link in enumerate(links):
        if i % num_process == 0:
            for sub_process in subprocesses:
                sub_process.wait()
            subprocesses = []
        subprocesses.append(
            subprocess.Popen(
                f"wget -O {os.path.join(path_out, files[i])} {link}",
                shell=True))

    for sub_process in subprocesses:
        sub_process.wait()

    return


if __name__ == '__main__':
    time_start = time.time()
    # simulated browser
    headers = ("User-Agent",
               "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 ("
               "KHTML, like Gecko) Chrome/70.0.3538.77 Safari/537.36")
    opener = urllib.request.build_opener()
    opener.addheaders = [headers]
    socket.setdefaulttimeout(2000)

    # narrow_url = 'https://egg2.wustl.edu/roadmap/data/byFileType/peaks/' \
    #              'consolidated/narrowPeak/ucsc_compatible/'
    # path_out = '/home/zy/driver_mutation/data/RoadMap/narrow_peak_chip_DHS'
    # narrow_files, narrow_links = get_links_roadmap(narrow_url)
    # download_data(narrow_files, narrow_links, path_out, 10)

    # broad_url = 'https://egg2.wustl.edu/roadmap/data/byFileType/peaks/' \
    #             'consolidated/narrowPeak/ucsc_compatible/'
    # path_broad = '/home/zy/driver_mutation/data/RoadMap/broad_peak_chip_DHS'
    # broad_files, broad_links = get_links_roadmap(broad_url)
    # download_data(broad_files, broad_links, path_broad, 20)

    signal_p_url = 'https://egg2.wustl.edu/roadmap/data/byFileType/signal/' \
                   'consolidated/macs2signal/pval/'
    path_signal_p = '/home/zy/driver_mutation/data/RoadMap/signal_p_chip_DHS'
    signal_p_files, signal_p_links = get_links_roadmap(signal_p_url, '.bigwig')
    download_data(signal_p_files, signal_p_links, path_signal_p, 20)

    # signal_fc_url = 'https://egg2.wustl.edu/roadmap/data/byFileType/signal/' \
    #                 'consolidated/macs2signal/foldChange/'
    # path_signal_fc = '/home/zy/driver_mutation/data/RoadMap/signal_fc_chip_DHS'
    # signal_fc_files, signal_fc_links = \
    #     get_links_roadmap(signal_fc_url, '.bigwig')
    # download_data(signal_fc_files, signal_fc_links, path_signal_fc, 20)

    # rnaseq_url = \
    #     'https://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/'
    # path_rnaseq = '/home/zy/driver_mutation/data/RoadMap/Rnaseq'
    # rnaseq_files, rnaseq_links = get_links_roadmap(rnaseq_url)
    # download_data(rnaseq_files, rnaseq_links, path_rnaseq, 20)

    # methy_wgbs_url = 'https://egg2.wustl.edu/roadmap/data/byDataType/' \
    #                  'dnamethylation/WGBS/FractionalMethylation_bigwig/'
    # path_methy_wgbs = '/home/zy/driver_mutation/data/RoadMap/methy/WGBS'
    # methy_wgbs_files, methy_wgbs_links = get_links_roadmap(methy_wgbs_url)
    # download_data(methy_wgbs_files, methy_wgbs_links, path_methy_wgbs, 20)

    # methy_rrbs_url = 'https://egg2.wustl.edu/roadmap/data/byDataType/' \
    #                  'dnamethylation/RRBS/FractionalMethylation_bigwig/'
    # path_methy_rrbs = '/home/zy/driver_mutation/data/RoadMap/methy/RRBS'
    # methy_rrbs_files, methy_rrbs_links = get_links_roadmap(
    #     methy_rrbs_url, '.bigwig')
    # download_data(methy_rrbs_files, methy_rrbs_links, path_methy_rrbs, 20)

    # methy_mcrf_url = 'https://egg2.wustl.edu/roadmap/data/byDataType/' \
    #                  'dnamethylation/mCRF/FractionalMethylation_bigwig/'
    # path_methy_mcrf = '/home/zy/driver_mutation/data/RoadMap/methy/mCRF'
    # methy_mcrf_files, methy_mcrf_links = get_links_roadmap(
    #     methy_mcrf_url, '.bigwig')
    # download_data(methy_mcrf_files, methy_mcrf_links, path_methy_mcrf, 20)

    time_end = time.time()
    print(time_end - time_start)
