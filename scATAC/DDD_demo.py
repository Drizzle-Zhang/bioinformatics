#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: DDD_demo.py
# @time: 2021/8/23 23:06


import os
import urllib
import torch
import torch.nn as nn
import torch.nn.init as init
import torch.nn.functional as F
import torch.utils.data as data
import torch.optim as optim
import numpy as np
import scipy.sparse as sp
from zipfile import ZipFile
from sklearn.model_selection import train_test_split
import pickle
import pandas as pd
import torch_scatter #注意：torch_scatter 安装时编译需要用到cuda
from collections import Counter


class DDDataset(object):
    # 数据集下载链接
    # url = "https://ls11-www.cs.tu-dortmund.de/people/morris/graphkerneldatasets/DD.zip"
    def __init__(self, url, data_root="data", train_size=0.8):
        self.url = url
        self.data_root = data_root
        # 下载 并解压
        self.maybe_download()
        sparse_adjacency, node_labels, graph_indicator, graph_labels = self.read_data()
        # 把coo格式转换为csr 进行稀疏矩阵运算
        self.sparse_adjacency = sparse_adjacency.tocsr()
        self.node_labels = node_labels
        self.graph_indicator = graph_indicator
        self.graph_labels = graph_labels

        self.train_index, self.test_index = self.split_data(train_size)
        # 得到训练集中所有图对应的类别标签
        self.train_label = graph_labels[self.train_index]
        # 得到测试集中所有图对应的类别标签
        self.test_label = graph_labels[self.test_index]

    def split_data(self, train_size):
        unique_indicator = np.asarray(list(set(self.graph_indicator)))
        # 随机划分训练集和测试集 得到各自对应的图索引   （一个图代表一条数据）
        train_index, test_index = train_test_split(unique_indicator,
                                                   train_size=train_size,
                                                   random_state=1234)
        return train_index, test_index

    def __getitem__(self, index):

        mask = self.graph_indicator == index
        # 得到图索引为index的图对应的所有节点(索引)
        graph_indicator = self.graph_indicator[mask]
        # 每个节点对应的特征标签
        node_labels = self.node_labels[mask]
        # 该图对应的类别标签
        graph_labels = self.graph_labels[index]
        # 该图对应的邻接矩阵
        adjacency = self.sparse_adjacency[mask, :][:, mask]
        return adjacency, node_labels, graph_indicator, graph_labels

    def __len__(self):
        return len(self.graph_labels)

    def read_data(self):
        # 解压后的路径
        data_dir = os.path.join(self.data_root, "DD")
        print("Loading DD_A.txt")
        # 从txt文件中读取邻接表(每一行可以看作一个坐标，即邻接矩阵中非0值的位置)  包含所有图的节点
        adjacency_list = np.genfromtxt(os.path.join(data_dir, "DD_A.txt"),
                                       dtype=np.int64, delimiter=',') - 1
        print("Loading DD_node_labels.txt")
        # 读取节点的特征标签（包含所有图） 每个节点代表一种氨基酸 氨基酸有20多种，所以每个节点会有一个类型标签 表示是哪一种氨基酸
        node_labels = np.genfromtxt(os.path.join(data_dir, "DD_node_labels.txt"),
                                    dtype=np.int64) - 1
        print("Loading DD_graph_indicator.txt")
        # 每个节点属于哪个图
        graph_indicator = np.genfromtxt(os.path.join(data_dir, "DD_graph_indicator.txt"),
                                        dtype=np.int64) - 1
        print("Loading DD_graph_labels.txt")
        # 每个图的标签 （2分类 0，1）
        graph_labels = np.genfromtxt(os.path.join(data_dir, "DD_graph_labels.txt"),
                                     dtype=np.int64) - 1
        # 节点数 （包含所有图的节点）
        num_nodes = len(node_labels)
        # 通过邻接表生成邻接矩阵  （包含所有的图）稀疏存储节省内存（coo格式 只存储非0值的行索引、列索引和非0值）
        # coo格式无法进行稀疏矩阵运算
        sparse_adjacency = sp.coo_matrix((np.ones(len(adjacency_list)),
                                          (adjacency_list[:, 0], adjacency_list[:, 1])),
                                         shape=(num_nodes, num_nodes), dtype=np.float32)
        print("Number of nodes: ", num_nodes)
        return sparse_adjacency, node_labels, graph_indicator, graph_labels

    def maybe_download(self):
        save_path = os.path.join(self.data_root)
        # 本地不存在 则下载
        if not os.path.exists(save_path):
            self.download_data(self.url, save_path)
        # 对数据集压缩包进行解压
        if not os.path.exists(os.path.join(self.data_root, "DD")):
            zipfilename = os.path.join(self.data_root, "DD.zip")
            with ZipFile(zipfilename, "r") as zipobj:
                zipobj.extractall(os.path.join(self.data_root))
                print("Extracting data from {}".format(zipfilename))

    @staticmethod
    def download_data(url, save_path):
        """数据下载工具，当原始数据不存在时将会进行下载"""
        print("Downloading data from {}".format(url))
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        # 下载数据集压缩包 保存在本地
        data = urllib.request.urlopen(url)
        filename = "DD.zip"
        with open(os.path.join(save_path, filename), 'wb') as f:
            f.write(data.read())
        return True


