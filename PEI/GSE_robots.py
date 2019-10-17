import urllib.request
import urllib.error
import requests
from lxml import etree
import socket
import time
import random
import multiprocessing
import gzip,tarfile,zipfile,xlrd
import subprocess
import logging
import os
import re
import ssl

headers=("User-Agent","Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.77 Safari/537.36")
# proxy_addr="58.240.220.86:53281"
# proxy=urllib.request.ProxyHandler({"https":proxy_addr})
opener=urllib.request.build_opener()
opener.addheaders=[headers]
# ssl._create_default_https_context = ssl._create_unverified_context
# urllib.request.install_opener(opener)
socket.setdefaulttimeout(2000)

date="0505"
project_base_dir_path="/home/disk/scRef/"
# pmid_file_path="pmid_test.txt"
pmid_file_path=os.path.join(project_base_dir_path,"pmid_file_"+date+"_contain_gse_plus.txt")
base_dir_path="spider_pmid_folder"
cell_set_file_path="M_H_cell_type_set.txt"
organ_file_path="tissue_organ.txt"
filtered_with_availabledata_name="filtered_gse_pmid_"+date+"_available_matrix.txt"    #这个是保存了含有可用细胞类型和可用矩阵的pmid
filtered_with_availabledata_matrix="filtered_gse_file_"+date+"_available.txt"    #这个是包含了可用矩阵和细胞类型的路径及种类
cell_type_reference_file_path_wy="/serverDNA/wuyu/CL_names.txt"
cell_type_reference_file_path_yzj="/home/disk/scRef/M_H_celltype.txt"
human1_file_name="refGenome/Homo_sapiens.GRCh37.75.txt"
human2_file_name="refGenome/Homo_sapiens.GRCh38.87.txt"
mouse1_file_name="refGenome/Mus_musculus.GRCm38.87.txt"
mouse2_file_name="refGenome/Mus_musculus.NCBIM37.67.txt"
pattern_cell_set = '_|\.|>|,|/|\(|\)|\[|\]|"'
if not os.path.exists(base_dir_path):
    os.mkdir(base_dir_path)
logging.basicConfig(level=logging.DEBUG,
                format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                datefmt='%a, %d %b %Y %H:%M:%S',
                filename=os.path.join(project_base_dir_path,"log_file.log"),
                filemode='w')
# id_lst=[]
# with open(pmid_file_path) as pmid_list_file:
#     for line in pmid_list_file:
#         inter=line.strip().split()
#         for GSE_index in range(1,len(inter)):
#             pmid=inter[0];GSE_AN=inter[GSE_index]
#             id_lst.append([pmid,GSE_AN])


def DownOneFile(srcUrl, localFile):
    print('%s\n --->>>\n  %s' % (srcUrl, localFile))

    startTime = time.time()
    with requests.get(srcUrl, stream=True) as r:
        contentLength = int(r.headers['content-length'])
        line = 'content-length: %dB/ %.2fKB/ %.2fMB'
        line = line % (contentLength, contentLength / 1024, contentLength / 1024 / 1024)
        print(line)
        downSize = 0
        with open(localFile, 'wb') as f:
            for chunk in r.iter_content(8192):
                if chunk:
                    f.write(chunk)
                downSize += len(chunk)
                line = '%d KB/s - %.2f MB， 共 %.2f MB'
                line = line % (
                downSize / 1024 / (time.time() - startTime), downSize / 1024 / 1024, contentLength / 1024 / 1024)
                print(line, end='\r')
                if downSize >= contentLength:
                    break
        timeCost = time.time() - startTime
        line = 'Total time: %.2f s, average speed: %.2f KB/s'
        line = line % (timeCost, downSize / 1024 / timeCost)
        print(line)

def callbackfunc(blocknum, blocksize, totalsize):
    '''
    @blocknum: data block downloaded
    @blocksize: size of data block
    @totalsize: total size of file
    '''
    percent = 100.0 * blocknum * blocksize / totalsize
    if percent > 100:
        percent = 100
    print("%.2f%%"% percent)

def getANDsave(id):
    AN=id[1]
    print("working on " + AN)
    if AN.startswith("GSE"):
        return (getANDsave_type1(id))
    elif AN.startswith("E-MTAB"):
        return (getANDsave_type2(id))
    else:
        pass

def getANDsave_type1(id):
    pmid=id[0];GSE_AN=id[1]
    pmid_dir_path=os.path.join(base_dir_path,pmid)
    if not os.path.exists(pmid_dir_path):
        os.mkdir(pmid_dir_path)
    total_print_str = ""
    url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="+GSE_AN
    data=urllib.request.urlopen(url)
    data_str = str(data.read())
    html=etree.HTML(data_str)
    # obj=html.xpath('//tr[@valign="top"]/*/a/@href')
    # print_str=pmid+"\t"+GSE_AN+"\t"+"\t".join([dl_link for dl_link in obj if dl_link.startswith("/geo") or dl_link.startswith("ftp://")])
    # print(print_str)
    multi_link=html.xpath('//table[@cellpadding="2" and @cellspacing="0"]//tr[@valign="top"]//td[@colspan="2"]')
    if multi_link:
        for item in multi_link:
            words=item.xpath('text()')
            if words and words[0].startswith('This SuperSeries is composed of the following '):
                href_base=html.xpath('//table[@cellpadding="2" and @cellspacing="0"]//td[@valign="top"]/a')
                for base_i in href_base:
                    base_i_text=base_i.xpath('text()')
                    if base_i_text:
                        if base_i_text[0].startswith("GSE"):
                            # base_i_href=base_i.xpath('@href')
                            sub_gse_link=[pmid,base_i_text[0]]
                            total_print_str+=getANDsave_type1(sub_gse_link)
                            # getANDsave_type1((sub_gse_link))

    obj=html.xpath('//table[@cellpadding="2" and @cellspacing="2" and @width="600"]//tr[@valign="top"]/td[@bgcolor="#DEEBDC" or @bgcolor="#EEEEEE"]')

    index=0
    print_str=""
    for td in obj:
        if index==4:
            total_print_str+=print_str
            print_str=""
            index=0
            downloadSizeFlag = False
            tmp = td.xpath('text()')
            # print_str += "file name\t" + tmp[0] + "\n"
            print_str+=tmp[0]+"|"
            file_name=tmp[0]
        elif index==0:
            downloadSizeFlag=False    #判断此时大小是否可以下载
            # print_str+=GSE_AN+"\n"+"\n"
            tmp = td.xpath('text()')
            # print_str += "file name\t" + tmp[0] + "\n"
            print_str += tmp[0] + "|"
            file_name=tmp[0]
        elif index==1:
            tmp=td.xpath('text()')
            size=tmp[0].lower()
            if size.endswith("gb") and int(size[0])>0:
                downloadSizeFlag=False
            else:
                downloadSizeFlag=True
            # print_str+="size is \t"+tmp[0]+"\n"
        elif index==2:
            if downloadSizeFlag==True:
                tmp=td.xpath('a/@href')
                # print_str+="downlinks"+"\n"
                dl_link4down=""
                for dl_link in tmp:
                    if dl_link.startswith("/geo"):
                        # print_str+="https://www.ncbi.nlm.nih.gov"+dl_link+"\n"
                        dl_link4down="https://www.ncbi.nlm.nih.gov"+dl_link
                        break
                    elif dl_link.startswith("ftp://"):
                        # print_str+=dl_link+"\n"
                        dl_link4down=dl_link
                if not os.path.exists(os.path.join(pmid_dir_path,GSE_AN)):
                    os.mkdir(os.path.join(pmid_dir_path,GSE_AN))
                if file_name.endswith(".tar") or file_name.endswith(".gz") or file_name.endswith(".zip") or file_name.endswith(".bz2"):
                    decompress_file_name=".".join(file_name.split(".")[:-1])
                elif file_name.endswith(".tar.gz"):
                    decompress_file_name=".".join(file_name.split(".")[:-2])
                else:
                    decompress_file_name=file_name
                if not os.path.exists(os.path.join(pmid_dir_path,GSE_AN,file_name)) and not os.path.exists(os.path.join(pmid_dir_path,GSE_AN,decompress_file_name)):
                    DownOneFile(dl_link4down,os.path.join(pmid_dir_path,GSE_AN,file_name))
        elif index==3:
            tmp=td.xpath('text()')
            # print_str+="file type\t"+tmp[0]+"\n"
        index+=1
    total_print_str += print_str
    random_sleep = random.randint(5, 15)
    time.sleep(random_sleep)
    return (total_print_str)

def getANDsave_type2(id):
    pmid = id[0];
    emtab_AN = id[1]
    pmid_dir_path = os.path.join(base_dir_path, pmid)
    if not os.path.exists(pmid_dir_path):
        os.mkdir(pmid_dir_path)
    total_print_str = ""
    url = "https://www.ebi.ac.uk/arrayexpress/experiments/" + emtab_AN+"/files/"
    data = urllib.request.urlopen(url)
    data_str = str(data.read())
    html = etree.HTML(data_str)
    obj = html.xpath('//table[@border="0" and @cellpadding="0" and @cellspacing="0"]//tr')
    print_str = ""
    downloadSizeFlag = False
    for td in obj:
        tmp = td.xpath('td[@class="col_size"]/div/text()')
        if not tmp:
            continue
        size = tmp[0].lower()
        if size.endswith("gb") and int(size[0]) > 0:
            downloadSizeFlag = False
        else:
            downloadSizeFlag = True
        if downloadSizeFlag == True:
            tmp = td.xpath('td[@class="col_name"]//a/@href')
            file_name_tmp=td.xpath('td[@class="col_name"]//a/text()')
            dl_link=tmp[0]
            file_name=file_name_tmp[0]
            if not file_name.endswith(".fastq.gz") and not file_name.endswith(".fq.gz"):
                dl_link4down = "https://www.ebi.ac.uk" + dl_link
                if not os.path.exists(os.path.join(pmid_dir_path, emtab_AN)):
                    os.mkdir(os.path.join(pmid_dir_path, emtab_AN))
                if file_name.endswith(".tar") or file_name.endswith(".gz") or file_name.endswith(".zip") or file_name.endswith(".bz2"):
                    decompress_file_name=".".join(file_name.split(".")[:-1])
                elif file_name.endswith(".tar.gz"):
                    decompress_file_name=".".join(file_name.split(".")[:-2])
                else:
                    decompress_file_name=file_name
                if not os.path.exists(os.path.join(pmid_dir_path,emtab_AN,file_name)) and not os.path.exists(os.path.join(pmid_dir_path,emtab_AN,decompress_file_name)):
                    DownOneFile(dl_link4down, os.path.join(pmid_dir_path, emtab_AN, file_name))
    random_sleep = random.randint(5, 15)
    time.sleep(random_sleep)
    return (total_print_str)

def process_excel_file(type,abso_path_xls):
    try:
        data = xlrd.open_workbook(abso_path_xls)
        sheet_names_list = data.sheet_names()
        return_name_list=[]    #产生的所有sheet txt文件
        for sheet_name in sheet_names_list:
            table = data.sheet_by_name(sheet_name)
            nrows = table.nrows
            if nrows<=1:
                continue
            abso_path_txt = abso_path_xls.replace(type, "_"+sheet_name+".txt")
            return_name_list.append(abso_path_txt)
            with open(abso_path_txt, "w") as out_file:
                for i in range(nrows - 1):
                    out_file.write("\t".join([str(item) for item in table.row_values(i)]) + "\n")
                out_file.write("\t".join([str(item) for item in table.row_values(nrows - 1)]))
        if os.path.exists(abso_path_xls):
            os.remove(abso_path_xls)
        return return_name_list
    except xlrd.biffh.XLRDError:
        logging.debug(abso_path_xls+" has XLRDError")
        return []

def decompress_file(file_path):
    #要求file_path为绝对路径
    if file_path.endswith(".gz"):
        original_path = file_path.replace(".gz", "")
        if not os.path.exists(original_path):
            with gzip.GzipFile(file_path) as f_c:
                with open(original_path,"wb+") as f_d:
                    file_content=f_c.read()
                    f_d.write(file_content)
        if os.path.exists(file_path):
            os.remove(file_path)
        return original_path
    elif file_path.endswith(".tar"):
        original_path = file_path.replace(".tar", "")
        if not os.path.exists(original_path):
            with tarfile.open(file_path) as f_t:
                os.mkdir(original_path)
                f_t.extractall(path=original_path)
        if os.path.exists(file_path):
            os.remove(file_path)
        return original_path
    elif file_path.endswith("zip"):
        original_path = file_path.replace(".zip", "")
        if not os.path.exists(original_path):
            os.mkdir(original_path)
            try:
                with zipfile.ZipFile(file_path) as f_z:
                    f_z.extractall(path=original_path)
            except zipfile.BadZipFile:
                subprocess.Popen("unzip -d "+original_path+" "+file_path,shell=True).wait()
        if os.path.exists(file_path):
            os.remove(file_path)
        return original_path
    elif file_path.endswith("rar"):
        return file_path
    elif file_path.endswith("bz2"):
        original_path=file_path.replace(".bz2","")
        if not os.path.exists(original_path):
            subprocess.Popen("bunzip2 -k "+file_path,shell=True).wait()
        if os.path.exists(file_path):
            os.remove(file_path)
        return original_path

    elif file_path.endswith(".xls"):
        return (process_excel_file(".xls",file_path))

    elif file_path.endswith(".xlsx"):
        return (process_excel_file(".xlsx",file_path))
    else:
        return file_path

def getAllFiles(gse_base_path,all_file_abso_list):
    if os.path.isdir(gse_base_path):
        gse_sub_file_list=os.listdir(gse_base_path)
        for gse_sub_file_name in gse_sub_file_list:
            gse_sub_file_path=os.path.join(gse_base_path, gse_sub_file_name)
            decompressd_file_path=decompress_file(gse_sub_file_path)
            if decompressd_file_path==None:
                print(gse_sub_file_path)
            elif isinstance(decompressd_file_path,list):
                all_file_abso_list.extend(decompressd_file_path)
            elif os.path.isdir(decompressd_file_path):
                getAllFiles(decompressd_file_path,all_file_abso_list)
            elif decompressd_file_path.endswith(".tar") or decompressd_file_path.endswith(".gz") or decompressd_file_path.endswith(
                    ".zip") or decompressd_file_path.endswith(".bz2"):
                getAllFiles(decompressd_file_path,all_file_abso_list)
            else:
                all_file_abso_list.append(decompressd_file_path)
    else:
        gse_sub_file_path = gse_base_path
        decompressd_file_path = decompress_file(gse_sub_file_path)
        if decompressd_file_path == None:
            print(gse_sub_file_path)
        elif isinstance(decompressd_file_path, list):
            all_file_abso_list.extend(decompressd_file_path)
        elif os.path.isdir(decompressd_file_path):
            getAllFiles(decompressd_file_path, all_file_abso_list)
        elif decompressd_file_path.endswith(".tar") or decompressd_file_path.endswith(
                ".gz") or decompressd_file_path.endswith(
                ".zip") or decompressd_file_path.endswith(".bz2"):
            getAllFiles(decompressd_file_path, all_file_abso_list)
        else:
            all_file_abso_list.append(decompressd_file_path)
    return all_file_abso_list

def isNum(string):
    if string[0]=="-":
        string=string[1:]
    if string.isdigit():
        return True
    inter=string.split(".")
    if len(inter)==2:
        return (inter[0].isdigit() and inter[1].isdigit())
    return False
def confirm2colFile(inter):
    if len(inter)!=2:
        return False
    return isNum(inter[1])
def getNum(line,sep):
    num_shuzi_i=0
    num_gene_i=0
    num_shuzi_t=0
    inter=line.strip().split(sep)
    if len(inter)>100:
        inter_100=inter[:100]
        num_shuzi_t += len(inter_100)
        for item in inter_100:
            if item.startswith('"'):
                item_fix = item[1:-1]
            else:
                item_fix = item
            if item_fix and isNum(item_fix):
                num_shuzi_i+=1
            elif item_fix in gene_name_set or "ENSMUSG" in item_fix or "ENSG" in item_fix or \
                        "Rik" in item_fix:
                num_gene_i+=1
    elif len(inter)==2:
        gene=inter[0];num=inter[1]
        if gene.startswith('"'):
            gene_fix = gene[1:-1]
        else:
            gene_fix = gene
        if num and isNum(num):
            num_shuzi_i+=1
        elif gene_fix in gene_name_set or "ENSMUSG" in gene_fix or "ENSG" in gene_fix or \
        "Rik" in gene_fix:
            num_gene_i += 1
    return ([num_shuzi_i,num_gene_i,num_shuzi_t])

def isMatrixFile(decompressed_file_path):
    # 出现的数字的个数,\t分隔
    num_shuzi_c_t=0
    num_shuzi_t_t=0.1
    # 逗号分隔
    num_shuzi_c_c = 0
    num_shuzi_t_c = 0.1
    # 空格分隔
    num_shuzi_c_s = 0
    num_shuzi_t_s = 0.1
    # 出现的基因名称的个数
    num_gene_c_t=0
    num_gene_c_c=0
    num_gene_c_s=0

    cut_i=0
    with open(decompressed_file_path) as decompressed_file:
        while cut_i<=100:
            line=next(decompressed_file)
            # Tab键分隔
            num_shuzi_i,num_gene_i,add_item=getNum(line,"\t")
            num_shuzi_c_t+=num_shuzi_i
            num_gene_c_t+=num_gene_i
            num_shuzi_t_t += add_item
            # 逗号分隔
            num_shuzi_i, num_gene_i,add_item = getNum(line, ",")
            num_shuzi_c_c+=num_shuzi_i
            num_gene_c_c+=num_gene_i
            num_shuzi_t_c += add_item
            # 空格分隔
            num_shuzi_i,num_gene_i,add_item=getNum(line," ")
            num_shuzi_c_s+=num_shuzi_i
            num_gene_c_s+=num_gene_i
            num_shuzi_t_s+=add_item
            cut_i+=1
    fraction_shuzi = max(num_shuzi_c_t / num_shuzi_t_t,num_shuzi_c_c/num_shuzi_t_c,num_shuzi_c_s/num_shuzi_t_s)
    Num_gene=max(num_gene_c_t,num_gene_c_c,num_gene_c_s)
    if fraction_shuzi >= 0.9 and Num_gene>=90:
        return True
    else:
        return False

def isAvailableMatrix(gse_sub_file_abso_path):
    print_info=""
    # 表达矩阵
    decompressed_file_path = gse_sub_file_abso_path
    #判断是否存在着一个文件夹包含了多种文件的形式
    # 如果是文件怎么判断是表达矩阵
    try:
        flag1= isMatrixFile(decompressed_file_path)
        if flag1:
            print_info+=gse_sub_file_abso_path + " has normal matrix file\n"
            # print(gse_sub_file_abso_path + " has normal matrix file")
            decompressed_file_base_name = os.path.basename(decompressed_file_path)
            last_dir_path="/".join(decompressed_file_path.split("/")[:-1])
            if len(os.listdir(last_dir_path))>100:
                tmp_inter = re.split(pattern_cell_set, decompressed_file_base_name)
                for tmp_item in tmp_inter:
                    if len(tmp_item) > 0:
                        tmp_item = tmp_item.strip()
                        if tmp_item:
                            a, b = neglectCaps(tmp_item)
                            if a in all_cell_type_filter_set or b in all_cell_type_filter_set:
                                print_info += gse_sub_file_abso_path + "\t"+tmp_item+ " has cell type file\n"
                                break
    except StopIteration:
        pass
    except UnicodeDecodeError:
        logging.debug(decompressed_file_path + " has Unicode decode error")
    # Seurat格式_1
    # 是否有三列整数文件
    try:
        cut_i=0
        with open(decompressed_file_path) as gse_sub_file:
            for j in range(100):
                next(gse_sub_file)
            for line in gse_sub_file:
                if cut_i>200:
                    break
                inter = line.strip().split(" ")
                if len(inter) == 3:
                    for sub_num in inter:
                        int(sub_num)
                else:
                    raise ValueError
                cut_i+=1
        print_info+=gse_sub_file_abso_path + " has Seurat matrix file\n"
        available_flag=True
        # print(gse_sub_file_abso_path + " has Seurat matrix file")
    except (ValueError, StopIteration):
        pass
    except UnicodeDecodeError:
        logging.debug(decompressed_file_path+" has UnicodeDecodeError")
    # Seurat格式_2
    try:
        with open(decompressed_file_path) as gse_sub_file:
            for j in range(100):
                next(gse_sub_file)
            num_tmp = 0
            num_total = 0.1
            cut_i=0
            for line in gse_sub_file:
                if cut_i>200:
                    break
                inter = line.strip().split("\t")
                if inter and len(inter) <= 2:
                    if inter[0].startswith('"'):
                        inter_fix = inter[0][1:-1]
                    else:
                        inter_fix = inter[0]
                    if len(inter)==2:
                        if inter[1].startswith('"'):
                            inter_fix2=inter[1][1:-1]
                        else:
                            inter_fix2=inter[1]
                    else:
                        inter_fix2="1"
                    if inter_fix in gene_name_set or "ENSMUSG" in inter_fix or "ENSG" in inter_fix or \
                            "Rik" in inter_fix:
                        if not isNum(inter_fix2):
                            num_tmp += 1
                num_total += 1
                cut_i+=1
            if num_tmp / num_total > 0.8:
                print_info+=gse_sub_file_abso_path + " has Seurat gene file\n"
                # print(gse_sub_file_abso_path + " has Seurat gene file")
    except StopIteration:
        pass
    except UnicodeDecodeError:
        logging.debug(decompressed_file_path+" has UnicodeDecodeError")
    # 细胞类型注释文件
    try:
        threshold=0
        cell_type_set=set()
        list4search=[]
        with open(decompressed_file_path) as gse_sub_file:
            if decompressed_file_path.endswith(".txt"):
                count_row=0
                for line in gse_sub_file:
                    if count_row==500:
                        break
                    count_row+=1
                    list4search.extend(line.strip().split("\t"))
                for item_t in list4search:
                    tmp_inter=re.split(pattern_cell_set,item_t)
                    for tmp_item in tmp_inter:
                        if len(tmp_item)>0:
                            tmp_item = tmp_item.strip()
                            if tmp_item:
                                a,b=neglectCaps(tmp_item)
                                if a in all_cell_type_filter_set or b in all_cell_type_filter_set:
                                    cell_type_set.add(tmp_item)
                                    threshold += 1
                if threshold >= 2 and len(cell_type_set)>=2:
                    print_info += gse_sub_file_abso_path + "\t"+"\t".join(cell_type_set)+" threshold " + str(threshold) + " has cell type file\n"
            elif decompressed_file_path.endswith(".csv"):
                count_row=0
                for line in gse_sub_file:
                    if count_row==500:
                        break
                    count_row+=1
                    list4search.extend(line.strip().split(","))
                for item_c in list4search:
                    tmp_inter=re.split(pattern_cell_set,item_c)
                    for tmp_item in tmp_inter:
                        if len(tmp_item)>0:
                            tmp_item = tmp_item.strip()
                            if tmp_item:
                                a, b = neglectCaps(tmp_item)
                                if a in all_cell_type_filter_set or b in all_cell_type_filter_set:
                                    cell_type_set.add(tmp_item)
                                    threshold += 1
                if threshold >= 2 and len(cell_type_set)>=2:
                    print_info += gse_sub_file_abso_path +"\t"+"\t".join(cell_type_set)+ " threshold " + str(threshold) + " has cell type file\n"
    except StopIteration:
        pass
    except UnicodeDecodeError:
        logging.debug(decompressed_file_path + " has UnicodeDecodeError")
    # RDS格式
    decompressed_file_base_name = os.path.basename(decompressed_file_path)
    if decompressed_file_base_name.endswith(".RDS") or decompressed_file_base_name.endswith(".Rdata") or decompressed_file_base_name.endswith(".rds"):
        print_info+=gse_sub_file_abso_path + "has R data file\n"
    return print_info
def neglectCaps(string):
    a=string[0].upper()
    b=string[0].lower()
    return ([a+string[1:],b+string[1:]])
def isCelltypeInLine(string,threshold):
    cell_num_occur = 0
    for sub_cell_type in all_cell_type_list:
        if sub_cell_type in string:
            cell_num_occur+=1
            if cell_num_occur==threshold:
                return True
    return False

def combine_single_func(pmid_sub_dir):
    gse_list_dir = os.listdir(os.path.join(base_dir_path, pmid_sub_dir))
    file_list=[]
    for gse_sub_dir in gse_list_dir:
        # print("we are processing " + pmid_sub_dir + " " + gse_sub_dir)
        # if gse_sub_dir == "available_mat.log":
        #     os.remove(os.path.join(base_dir_path, pmid_sub_dir, gse_sub_dir))
        #     continue
        # 如下为了递归得到某个gse_sub_dir下所有的文件
        gse_base_path = os.path.join(base_dir_path, pmid_sub_dir, gse_sub_dir)  # 表示gse子文件夹的绝对路径
        gse_sub_all_file_list = []
        file_list.extend(getAllFiles(gse_base_path, gse_sub_all_file_list))
    return(file_list)
def multiprocess4download(id_list):
    n = 4
    pool = multiprocessing.Pool(processes=n)
    pool.map(getANDsave, id_list)
    pool.close()
    pool.join()
def multiprocess4infer(file_list):
    n=30
    pool = multiprocessing.Pool(processes=n)
    available_info=pool.map(isAvailableMatrix,file_list)
    pool.close()
    pool.join()
    return available_info
def multiprocess4combine(pmid_list):
    n=10
    pool=multiprocessing.Pool(processes=n)
    res_list_list=pool.map(combine_single_func,pmid_list)
    pool.close()
    pool.join()
    total_all_files_abso_path_list=[]
    for item in res_list_list:
        if item is not None:
            total_all_files_abso_path_list.extend(item)
    return total_all_files_abso_path_list
if __name__=="__main__":
    #这里可以控制返回的信息
    # id_lst=[["26949938","GSE75382"]]
    #控制是否下载数据
    # multiprocess4download(id_lst)
    #这里是为了产生基因名称的列表
    gene_name_set = set()
    gtf_file_list = [human1_file_name, human2_file_name, mouse1_file_name, mouse2_file_name]
    for gtf_sub_file_name in gtf_file_list:
        with open(os.path.join(project_base_dir_path, gtf_sub_file_name)) as gtf_file:
            for line in gtf_file:
                inter = line.strip().split(":")
                gene_name_set.add(inter[1])
    all_cell_type_list = []
    with open(cell_type_reference_file_path_wy) as wy_cell_file:
        for line in wy_cell_file:
            cell_name = line.strip()
            tmp_inter = re.split(pattern_cell_set, cell_name)
            for tmp_item in tmp_inter:
                if len(tmp_item) > 0:
                    tmp_item = tmp_item.strip()
                    if tmp_item and len(tmp_item) > 1:
                        all_cell_type_list.append(tmp_item)
    with open(cell_type_reference_file_path_yzj) as yzj_cell_file:
        for line in yzj_cell_file:
            cell_name=line.strip()
            tmp_inter = re.split(pattern_cell_set, cell_name)
            for tmp_item in tmp_inter:
                if len(tmp_item) > 0:
                    if tmp_item.endswith("?"):
                        tmp_item=tmp_item[:-1]
                    if tmp_item:
                        tmp_item=tmp_item.strip()
                        if tmp_item and len(tmp_item)>1:
                            all_cell_type_list.append(tmp_item)
    all_cell_type_set=set(all_cell_type_list)    #包含了所有All_ref中的细胞类型集合
    all_cell_type_set=all_cell_type_set.difference(gene_name_set)
    orgen_set=set()
    with open(organ_file_path) as organ_file:
        for line in organ_file:
            string1,string2=neglectCaps(line.strip())
            orgen_set.add(string1);orgen_set.add(string2)
    all_cell_type_set=all_cell_type_set.difference(orgen_set)
    all_cell_type_filter_set=set()
    pattern_filter1=re.compile("^[A-Za-z][0-9]+$")
    pattern_filter2=re.compile("^[0-9].+")
    for cell_type_item in all_cell_type_set:
        if pattern_filter1.search(cell_type_item) or pattern_filter2.search(cell_type_item):
           pass
        else:
            all_cell_type_filter_set.add(cell_type_item)
    with open(cell_set_file_path,"w") as cell_set_file:
        for item in all_cell_type_filter_set:
            cell_set_file.write(item+"\n")
    # print("Unknown" in all_cell_type_set)
    pmid_list_dir=os.listdir(base_dir_path)
    # pmid_list_dir=["30510133"]
    total_all_files_abso_path_list=multiprocess4combine(pmid_list_dir)

    available_info=multiprocess4infer(total_all_files_abso_path_list)

    # print(available_info)

    with open(os.path.join(project_base_dir_path,filtered_with_availabledata_matrix),"w") as available_file1:
            for print_item in available_info:
                available_file1.write(print_item)
    filtered_set_cell_gse = set()
    filtered_set_matrix_gse = set()
    with open(os.path.join(project_base_dir_path,filtered_with_availabledata_matrix)) as available_file1:
        for line in available_file1:
            inter = line.split("/")
            if line.strip().endswith("has cell type file"):
                filtered_set_cell_gse.add(inter[1]+"-"+inter[2])
            else:
                filtered_set_matrix_gse.add(inter[1]+"-"+inter[2])
    # 得到既有cell type又有矩阵的gse号，从而找到pmid
    filtered_set_gse=filtered_set_cell_gse & filtered_set_matrix_gse
    filtered_set_gse_pmid=set([item.split("-")[0] for item in filtered_set_gse])
    with open(os.path.join(project_base_dir_path, filtered_with_availabledata_name), "w") as pmid_matrix:
        for item in filtered_set_gse_pmid:
            pmid_matrix.write(item+"\n")
    print("filtered available pmid file: "+str(len(filtered_set_gse_pmid)))
    print("total pmid: "+str(len(pmid_list_dir)))



