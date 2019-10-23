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
import re


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


def get_links_geo_supp(url):
    web_html = urllib.request.urlopen(url)
    web_html_str = str(web_html.read())
    html = etree.HTML(web_html_str)
    pre_links = html.xpath(
        "//td[@bgcolor='#DEEBDC' or @bgcolor='#EEEEEE']/a/@href")
    ftp_links = [link for link in pre_links if link[:3] == 'ftp']

    return ftp_links


def get_links_3div_and_download(url, path_out, num_process):
    web_html = urllib.request.urlopen(url)
    web_html_str = str(web_html.read())
    html = etree.HTML(web_html_str)
    text = html.xpath("//*/text()")
    list_text = text[0].split('\\n')
    file_pattern = re.compile(r":[0-9][0-9] .+\\")
    ftp_links = []
    for line in list_text:
        if file_pattern.search(line):
            file = file_pattern.search(line).group()[4:-1]
            if file == 'IMR90_fibroblast,_TNF-\\xa5\\xe1_treated':
                # file = 'IMR90_fibroblast,_TNF-%A5%E1_treated'
                continue
            ftp_links.append(url + file + '/')

    out_links = []
    for link in ftp_links:
        web_html = urllib.request.urlopen(link)
        web_html_str = str(web_html.read())
        html = etree.HTML(web_html_str)
        text = html.xpath("//*/text()")
        list_text = text[0].split('\\n')
        for line in list_text:
            if file_pattern.search(line):
                out_links.append(
                    link + file_pattern.search(line).group()[4:-1])

    subprocesses = []
    for i, link in enumerate(out_links):
        folder = os.path.join(path_out, link.split('/')[-2])
        if not os.path.exists(folder):
            os.makedirs(folder)
        if i % num_process == 0:
            for sub_process in subprocesses:
                sub_process.wait()
            subprocesses = []
        if link[-6:] == '10.zip':
            subprocesses.append(
                subprocess.Popen(
                    "wget -P " + folder.replace('(', '\(').replace(')', '\)') +
                    " " + link.replace('(', '\(').replace(')', '\)'),
                    shell=True))
        random_sleep = random.randint(2, 5)
        time.sleep(random_sleep)
        # if i == 6:
        #     break

    for sub_process in subprocesses:
        sub_process.wait()

    return


def download_data_geo(links, path_out, num_process):
    subprocesses = []
    for i, link in enumerate(links):
        if i % num_process == 0:
            for sub_process in subprocesses:
                sub_process.wait()
            subprocesses = []
        subprocesses.append(
            subprocess.Popen(
                f"wget -P {path_out} {link}",
                shell=True))
        random_sleep = random.randint(5, 15)
        time.sleep(random_sleep)
        # if i == 6:
        #     break

    for sub_process in subprocesses:
        sub_process.wait()

    return


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


def decompress(path, num_process=20):
    files = os.listdir(path)
    subprocesses = []
    for i, file in enumerate(files):
        if i % num_process == 0:
            for sub_process in subprocesses:
                sub_process.wait()
            subprocesses = []
        if file[-3:] == '.gz':
            subprocesses.append(
                subprocess.Popen(
                    f"gzip -d {os.path.join(path, file)}",
                    shell=True))
        elif file[-4:] == '.bz2':
            subprocesses.append(
                subprocess.Popen(
                    f"bzip2 -d {os.path.join(path, file)}",
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

    # RoadMap
    # narrow_url = 'https://egg2.wustl.edu/roadmap/data/byFileType/peaks/' \
    #              'consolidated/narrowPeak/ucsc_compatible/'
    # path_narrow = '/home/zy/driver_mutation/data/RoadMap/narrow_peak_chip_DHS'
    # narrow_files, narrow_links = get_links_roadmap(narrow_url)
    # download_data(narrow_files, narrow_links, path_out, 10)
    # decompress(path_narrow, num_process=20)

    # broad_url = 'https://egg2.wustl.edu/roadmap/data/byFileType/peaks/' \
    #             'consolidated/narrowPeak/ucsc_compatible/'
    # path_broad = '/home/zy/driver_mutation/data/RoadMap/broad_peak_chip_DHS'
    # broad_files, broad_links = get_links_roadmap(broad_url)
    # download_data(broad_files, broad_links, path_broad, 20)
    # decompress(path_broad, num_process=40)

    # signal_p_url = 'https://egg2.wustl.edu/roadmap/data/byFileType/signal/' \
    #                'consolidated/macs2signal/pval/'
    # path_signal_p = '/home/zy/driver_mutation/data/RoadMap/signal_p_chip_DHS'
    # signal_p_files, signal_p_links = get_links_roadmap(signal_p_url, '.bigwig')
    # download_data(signal_p_files, signal_p_links, path_signal_p, 20)

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
    # decompress(path_rnaseq, num_process=40)

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

    # nature genetics promoter capture Hi-C
    # url_jung_ng_2019 = \
    #     'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86189'
    # links_jung_ng_2019 = get_links_geo_supp(url_jung_ng_2019)
    # path_jung_ng_2019 = \
    #     '/home/zy/driver_mutation/data/pcHiC/Jung_NG_2019'
    # download_data_geo(links_jung_ng_2019, path_jung_ng_2019, 5)
    # decompress(path_jung_ng_2019, 20)

    # Segway

    # 3DIV

    time_end = time.time()
    print(time_end - time_start)
