#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: GraphConv.py
# @time: 2021/9/9 11:58

from time import time
import os
from torch_geometric.data import InMemoryDataset, DataLoader
import networkx as nx
import numpy as np
from torch_geometric.utils import to_networkx
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch_geometric.nn as geo_nn
from captum.attr import Saliency, IntegratedGradients


class ATACGraphDataset(InMemoryDataset):
    def __init__(self, root, data_list=None):
        self.data_list = data_list
        super(ATACGraphDataset, self).__init__(root)
        self.data, self.slices = torch.load(self.processed_paths[0])

    @property
    def raw_file_names(self):
        return ['some_file_1']

    @property
    def processed_file_names(self):
        return ['data.pt']

    def download(self):
        pass

    def process(self):
        # Read data into huge `Data` list.
        data_list = self.data_list
        data, slices = self.collate(data_list)
        torch.save((data, slices), self.processed_paths[0])


class GCN(torch.nn.Module):
    def __init__(self, input_channels, output_channels, hidden_channels):
        super(GCN, self).__init__()
        torch.manual_seed(12345)

        # self.lin0 = nn.Linear(input_channels, hidden_channels//4)
        self.conv1 = geo_nn.GATConv(input_channels, hidden_channels,
                                    heads=4, dropout=0.6)
        # self.conv2 = geo_nn.GATConv(1*hidden_channels, 1*hidden_channels,
        #                             heads=1, dropout=0.7)
        # self.conv3 = geo_nn.GATConv(hidden_channels, hidden_channels)
        # self.conv4 = GraphConv(hidden_channels, hidden_channels)
        self.conv5 = geo_nn.GATConv(4*hidden_channels, output_channels, concat=False,
                                    heads=1, dropout=0.6)

        # self.lin1 = nn.Linear(hidden_channels, hidden_channels)
        # self.lin2 = nn.Sequential(
        #     nn.Linear(4*hidden_channels, hidden_channels),
        #     nn.ReLU(),
        #     nn.Linear(hidden_channels, output_channels)
        # )

    def forward(self, x, edge_index, edge_weight, batch):
        # x = self.lin0(x).relu()
        # edge_weight = None
        x = F.dropout(x, p=0.6, training=self.training)
        x = self.conv1(x=x, edge_index=edge_index)
        x = F.elu(x)
        # x = self.conv2(x, edge_index)
        # x = F.elu(x)
        # x = self.conv3(x, edge_index, edge_weight).relu()
        # x = self.conv4(x, edge_index, edge_weight).relu()
        # x = geo_nn.global_mean_pool(x, batch)
        # x = self.lin1(x).relu()
        x = F.dropout(x, p=0.6, training=self.training)
        # x = self.lin2(x)
        x = self.conv5(x, edge_index)
        x = geo_nn.global_mean_pool(x, batch)

        return F.log_softmax(x, dim=1)


def train(loader):
    model.train()
    for data in loader:  # Iterate in batches over the training dataset.
        data = data.to(device)
        out = model(data.x, data.edge_index, data.edge_attr, data.batch)  # Perform a single forward pass.
        loss = criterion(out, data.y)  # Compute the loss.
        loss.backward()  # Derive gradients.
        optimizer.step()  # Update parameters based on gradients.
        optimizer.zero_grad()  # Clear gradients.


def test(loader):
    model.eval()
    correct = 0
    for data in loader:  # Iterate in batches over the training/test dataset.c
        data = data.to(device)
        out = model(data.x, data.edge_index, data.edge_attr, data.batch)
        pred = out.argmax(dim=1)  # Use the class with highest probability.
        correct += int((pred == data.y).sum())  # Check against ground-truth labels.
    return correct / len(loader.dataset)  # Derive ratio of correct predictions.


def model_forward(edge_mask, data):
    batch = torch.zeros(data.x.shape[0], dtype=int).to(device)
    out = model(data.x, data.edge_index, edge_mask, batch)
    return out


def explain(method, data, target=0):
    input_mask = torch.ones(data.edge_index.shape[1]).requires_grad_(True).to(device)
    if method == 'ig':
        ig = IntegratedGradients(model_forward)
        mask = ig.attribute(input_mask, target=target,
                            additional_forward_args=(data,),
                            internal_batch_size=data.edge_index.shape[1])
    elif method == 'saliency':
        saliency = Saliency(model_forward)
        mask = saliency.attribute(input_mask, target=target,
                                  additional_forward_args=(data,))
    else:
        raise Exception('Unknown explanation method')

    edge_mask = np.abs(mask.cpu().detach().numpy())
    if edge_mask.max() > 0:  # avoid division by zero
        edge_mask = edge_mask / edge_mask.max()
    return edge_mask


def aggregate_edge_directions(edge_mask, data):
    edge_mask_dict = defaultdict(float)
    for val, u, v in list(zip(edge_mask, *data.edge_index)):
        u, v = u.item(), v.item()
        if u > v:
            u, v = v, u
        edge_mask_dict[(u, v)] += val
    return edge_mask_dict


if __name__ == '__main__':
    time_start = time()
    path_data_root = '/root/scATAC/ATAC_data/Forebrain/Mulqueen_human_brain'
    path_graph_input = os.path.join(path_data_root, 'input_graph')
    dataset_atac_graph = ATACGraphDataset(path_graph_input)
    torch.manual_seed(12345)
    dataset = dataset_atac_graph.shuffle()
    train_dataset = dataset[:1700]
    test_dataset = dataset[1700:]

    device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")
    model = GCN(input_channels=dataset.num_node_features,
                output_channels=dataset.num_classes, hidden_channels=16).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
    criterion = torch.nn.CrossEntropyLoss()

    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)

    time_start = time()
    for epoch in range(1, 101):
        train(train_loader)
        train_acc = test(train_loader)
        test_acc = test(test_loader)
        print(f'Epoch: {epoch:03d}, '
              f'Train Acc: {train_acc:.4f}, Test Acc: {test_acc:.4f}')

        optimizer = torch.optim.Adam(model.parameters(), lr=0.005)
    for epoch in range(1, 101):
        train(train_loader)
        train_acc = test(train_loader)
        test_acc = test(test_loader)
        print(f'Epoch: {epoch:03d}, Train Acc: {train_acc:.4f}, Test Acc: {test_acc:.4f}')

    optimizer = torch.optim.Adam(model.parameters(), lr=0.003)
    for epoch in range(1, 101):
        train(train_loader)
        train_acc = test(train_loader)
        test_acc = test(test_loader)
        print(f'Epoch: {epoch:03d}, Train Acc: {train_acc:.4f}, Test Acc: {test_acc:.4f}')

    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    for epoch in range(1, 101):
        train(train_loader)
        train_acc = test(train_loader)
        test_acc = test(test_loader)
        print(f'Epoch: {epoch:03d}, Train Acc: {train_acc:.4f}, Test Acc: {test_acc:.4f}')
    time_end = time()
    print(time_end - time_start)
    model = model.cpu()
    torch.cuda.empty_cache()
