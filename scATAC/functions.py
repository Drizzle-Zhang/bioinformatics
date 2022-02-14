# _*_ coding: utf-8 _*_
# @author: Drizzle_Zhang
# @file: functions.py
# @time: 2022/1/12 11:54

from time import time
import os
from torch_geometric.data import InMemoryDataset, Data, DataLoader
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GraphConv
from captum.attr import Saliency, IntegratedGradients
from collections import defaultdict
from multiprocessing import Pool
from functools import partial
from scipy.stats import kstest
import episcanpy.api as epi
import scanpy as sc
import pandas as pd
import anndata as ad
import pickle


class GCN(torch.nn.Module):
    def __init__(self, input_channels, output_channels, hidden_channels, num_nodes):
        super(GCN, self).__init__()
        torch.manual_seed(12345)

        self.conv1 = GraphConv(input_channels, hidden_channels//4)
        # self.conv2 = GraphConv(hidden_channels//4, hidden_channels)
        # self.conv3 = GraphConv(hidden_channels, hidden_channels//4)
        # self.conv4 = GraphConv(hidden_channels, hidden_channels)
        self.conv5 = GraphConv(hidden_channels//4, input_channels)

        self.lin1 = nn.Linear(num_nodes, num_nodes//5)
        self.lin2 = nn.Linear(num_nodes//5, num_nodes//25)
        self.lin3 = nn.Linear(num_nodes//25, num_nodes//125)
        self.lin4 = nn.Linear(num_nodes//125, output_channels)

    def forward(self, x, edge_index, batch, edge_weight=None):
        x = self.conv1(x, edge_index, edge_weight).relu()
        # x = self.conv2(x, edge_index, edge_weight).relu()
        # x = self.conv3(x, edge_index, edge_weight).relu()
        # x = self.conv4(x, edge_index, edge_weight).relu()
        x = self.conv5(x, edge_index, edge_weight)
        batch_size = len(torch.unique(batch))
        x = x.view(batch_size, x.shape[0]//batch_size)
        x = self.lin1(x).relu()
        x = F.dropout(x, p=0.5, training=self.training)
        x = self.lin2(x).relu()
        x = F.dropout(x, p=0.5, training=self.training)
        x = self.lin3(x).relu()
        x = F.dropout(x, p=0.5, training=self.training)
        x = self.lin4(x)
        return F.log_softmax(x, dim=1)


def explain(method, device, model, data):

    def model_forward(edge_mask, device, data, model):
        batch = torch.zeros(data.x.shape[0], dtype=int).to(device)
        out = model(data.x, data.edge_index, batch, edge_mask)
        return out

    data = data.to(device)
    target = data.y.to(device)
    input_mask = torch.ones(data.edge_index.shape[1]).requires_grad_(True).to(device)
    if method == 'ig':
        ig = IntegratedGradients(model_forward)
        mask = ig.attribute(input_mask, target=target,
                            additional_forward_args=(data,),
                            internal_batch_size=data.edge_index.shape[1])
    elif method == 'saliency':
        saliency = Saliency(model_forward)
        mask = saliency.attribute(input_mask, target=target,
                                  additional_forward_args=(device, data, model))
    else:
        raise Exception('Unknown explanation method')

    edge_mask = np.abs(mask.cpu().detach().numpy())
    if edge_mask.max() > 0:  # avoid division by zero
        edge_mask = edge_mask / edge_mask.max()

    edge_mask_dict = defaultdict(float)
    for val, u, v in list(zip(edge_mask, *data.edge_index)):
        u, v = u.item(), v.item()
        edge_mask_dict[(u, v)] += val

    return edge_mask_dict
