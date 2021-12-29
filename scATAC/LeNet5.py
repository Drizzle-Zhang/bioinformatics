#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: LeNet5.py
# @time: 2021/9/8 20:01

import torch
from torch.autograd import Variable
import torch.nn.functional as F
import torch.nn as nn
import pdb
import collections
from time import time
import numpy as np
import os
import sys
import tensorflow as tf


sys.path.insert(0, 'lib/')
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = "0"

if torch.cuda.is_available():
    print('cuda available')
    dtypeFloat = torch.cuda.FloatTensor
    dtypeLong = torch.cuda.LongTensor
    torch.cuda.manual_seed(1)
else:
    print('cuda not available')
    dtypeFloat = torch.FloatTensor
    dtypeLong = torch.LongTensor
    torch.manual_seed(1)

mint = tf.keras.datasets.mnist
(train_data, train_labels), (test_data, test_labels) = mint.load_data()
val_data = train_data[55000:]
val_labels = train_labels[55000:]
train_data = train_data[:55000]
train_labels = train_labels[:55000]

if __name__ == '__main__':
    time_start = time()

    time_end = time()
    print(time_end - time_start)
