#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: NNConv.py
# @time: 2021/9/9 17:43

from time import time
import os
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.datasets import MNISTSuperpixels
import torch_geometric.transforms as T
from torch_geometric.data import InMemoryDataset, DataLoader
from torch_geometric.utils import normalized_cut
import torch_geometric.nn as geo_nn
from torch_geometric.nn import (NNConv, graclus, max_pool, max_pool_x, global_mean_pool)


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


class SAGPool(nn.Module):
    def __init__(self, input_channels, output_channels, hidden_channels):
        super(SAGPool, self).__init__()
        torch.manual_seed(12345)
        self.lin1 = nn.Linear(input_channels, 8)
        self.gcn1 = geo_nn.GraphConv(8, hidden_channels)
        self.pool1 = geo_nn.SAGPooling(hidden_channels, 0.5)
        self.gcn2 = geo_nn.GraphConv(hidden_channels, hidden_channels)
        self.pool2 = geo_nn.SAGPooling(hidden_channels, 0.5)
        self.gcn3 = geo_nn.GraphConv(hidden_channels, hidden_channels)
        self.pool3 = geo_nn.SAGPooling(hidden_channels, 0.5)

        self.mlp = nn.Sequential(
            # nn.Linear(hidden_channels * 2, hidden_channels),
            # nn.ReLU(),
            nn.Linear(hidden_channels, hidden_channels // 2),
            nn.ReLU(),
            nn.Linear(hidden_channels // 2, output_channels))

    def forward(self, x, edge_index, edge_weight, batch):
        x = self.lin1(x)
        x = x.relu()
        gcn1 = F.relu(self.gcn1(x, edge_index, edge_weight))
        x1, edge_index1, edge_attr1, batch1, _, _ = \
            self.pool1(gcn1, edge_index, edge_weight, batch=batch)
        global_pool1 = geo_nn.global_mean_pool(x1, batch1)
        # global_pool1 = torch.cat(
        #     [geo_nn.global_mean_pool(x1, batch1),
        #      geo_nn.global_max_pool(x1, batch1)],
        #     dim=1)

        gcn2 = F.relu(self.gcn2(x1, edge_index1, edge_attr1))
        x2, edge_index2, edge_attr2, batch2, _, _ = \
            self.pool2(gcn2, edge_index1, edge_attr1, batch=batch1)
        global_pool2 = geo_nn.global_mean_pool(x2, batch2)
        # global_pool2 = torch.cat(
        #     [geo_nn.global_mean_pool(x2, batch2),
        #      geo_nn.global_max_pool(x2, batch2)],
        #     dim=1)

        gcn3 = F.relu(self.gcn3(x2, edge_index2, edge_attr2))
        x3, edge_index3, edge_attr3, batch3, _, _ = \
            self.pool3(gcn3, edge_index2, edge_attr2, batch=batch2)
        global_pool3 = geo_nn.global_mean_pool(x3, batch3)
        # global_pool3 = torch.cat(
        #     [geo_nn.global_mean_pool(x3, batch3),
        #      geo_nn.global_max_pool(x3, batch3)],
        #     dim=1)

        x = global_pool1 + global_pool2 + global_pool3
        x = self.mlp(x)
        return x


def train(loader):
    model.train()

    # if epoch == 50:
    #     for param_group in optimizer.param_groups:
    #         param_group['lr'] = 0.003
    #
    # if epoch == 100:
    #     for param_group in optimizer.param_groups:
    #         param_group['lr'] = 0.001

    for data in loader:
        data = data.to(device)
        out = model(data.x, data.edge_index, data.edge_attr, data.batch)  # Perform a single forward pass.
        loss = criterion(out, data.y)  # Compute the loss.
        optimizer.zero_grad()
        loss.backward()  # Derive gradients.
        optimizer.step()  # Update parameters based on gradients.


def test(loader):
    model.eval()
    correct = 0

    for data in loader:
        data = data.to(device)
        out = model(data.x, data.edge_index, data.edge_attr, data.batch)
        pred = out.argmax(dim=1)  # Use the class with highest probability.
        correct += int((pred == data.y).sum())  # Check against ground-truth labels.
    return correct / len(loader.dataset)


if __name__ == '__main__':
    time_start = time()
    path_data_root = '/root/scATAC/ATAC_data/Forebrain/Mulqueen_human_brain'
    path_graph_input = os.path.join(path_data_root, 'input_graph')
    dataset_atac_graph = ATACGraphDataset(path_graph_input)
    dataset = dataset_atac_graph.shuffle()
    train_dataset = dataset[:1700]
    test_dataset = dataset[1700:]
    device = torch.device("cuda:3" if torch.cuda.is_available() else "cpu")
    torch.manual_seed(123)

    model = SAGPool(input_channels=dataset.num_node_features,
                    output_channels=dataset.num_classes, hidden_channels=32).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
    criterion = torch.nn.CrossEntropyLoss()

    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)
    transform = T.Cartesian(cat=False)

    for epoch in range(1, 201):
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

    optimizer = torch.optim.Adam(model.parameters(), lr=0.0005)
    for epoch in range(1, 101):
        train(train_loader)
        train_acc = test(train_loader)
        test_acc = test(test_loader)
        print(f'Epoch: {epoch:03d}, Train Acc: {train_acc:.4f}, Test Acc: {test_acc:.4f}')

    time_end = time()
    print(time_end - time_start)
