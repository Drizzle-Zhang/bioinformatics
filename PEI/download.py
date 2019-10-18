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


def get_links_roadmap(url):
    # get download links
    web_html = urllib.request.urlopen(url)
    web_html_str = str(web_html.read())
    html = etree.HTML(web_html_str)
    files = html.xpath('//a/@href')
    files_filter = []
    links = []
    for file in files:
        if file[-3:] == '.gz':
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

    # bigwig_url = 'https://egg2.wustl.edu/roadmap/data/byFileType/peaks/' \
    #              'consolidated/narrowPeak/ucsc_compatible/'
    # path_bigwig = '/home/zy/driver_mutation/data/RoadMap/bigwig_peak_chip_DHS'
    # bigwig_files, bigwig_links = get_links_roadmap(bigwig_url)
    # download_data(bigwig_files, bigwig_links, path_bigwig, 20)

    rnaseq_url = \
        'https://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/'
    path_rnaseq = '/home/zy/driver_mutation/data/RoadMap/Rnaseq'
    rnaseq_files, rnaseq_links = get_links_roadmap(rnaseq_url)
    download_data(rnaseq_files, rnaseq_links, path_rnaseq, 20)

    # bigwig_url = 'https://egg2.wustl.edu/roadmap/data/byFileType/peaks/' \
    #              'consolidated/narrowPeak/ucsc_compatible/'
    # path_bigwig = '/home/zy/driver_mutation/data/RoadMap/bigwig_peak_chip_DHS'
    # bigwig_files, bigwig_links = get_links_roadmap(bigwig_url)
    # download_data(bigwig_files, bigwig_links, path_bigwig, 20)

    # methy_url = 'https://egg2.wustl.edu/roadmap/data/byFileType/peaks/' \
    #             'consolidated/narrowPeak/ucsc_compatible/'
    # path_methy = '/home/zy/driver_mutation/data/RoadMap/methy_peak_chip_DHS'
    # methy_files, methy_links = get_links_roadmap(methy_url)
    # download_data(methy_files, methy_links, path_methy, 20)

    time_end = time.time()
    print(time_end - time_start)
