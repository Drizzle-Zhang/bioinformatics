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


class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        nn1 = nn.Sequential(nn.Linear(1, 30), nn.ReLU(),
                            nn.Linear(30, dataset.num_features * 48))
        self.conv1 = NNConv(dataset.num_features, 48, nn1, aggr='mean')

        nn2 = nn.Sequential(nn.Linear(1, 30), nn.ReLU(),
                            nn.Linear(30, 48 * 96))
        self.conv2 = NNConv(48, 96, nn2, aggr='mean')

        self.fc1 = torch.nn.Linear(96, 128)
        self.fc2 = torch.nn.Linear(128, dataset.num_classes)

    def forward(self, data):
        data.x = F.elu(self.conv1(data.x, data.edge_index, data.edge_attr))
        cluster = graclus(data.edge_index,
                          torch.reshape(data.edge_attr,
                                        (data.edge_attr.shape[0],)),
                          data.x.size(0))
        data = max_pool(cluster, data)

        data.x = F.elu(self.conv2(data.x, data.edge_index, data.edge_attr))
        cluster = graclus(data.edge_index,
                          torch.reshape(data.edge_attr,
                                        (data.edge_attr.shape[0],)),
                          data.x.size(0))
        x, batch = max_pool_x(cluster, data.x, data.batch)

        x = global_mean_pool(x, batch)
        x = F.elu(self.fc1(x))
        x = F.dropout(x, training=self.training)
        return F.log_softmax(self.fc2(x), dim=1)


def train(epoch):
    model.train()

    # if epoch == 50:
    #     for param_group in optimizer.param_groups:
    #         param_group['lr'] = 0.003
    #
    # if epoch == 100:
    #     for param_group in optimizer.param_groups:
    #         param_group['lr'] = 0.001

    for data in train_loader:
        data = data.to(device)
        optimizer.zero_grad()
        F.nll_loss(model(data), data.y).backward()
        optimizer.step()


def test(loader):
    model.eval()
    correct = 0

    for data in loader:
        data = data.to(device)
        pred = model(data).max(1)[1]
        correct += pred.eq(data.y).sum().item()
    return correct / len(loader.dataset)


if __name__ == '__main__':
    time_start = time()
    path_data_root = '/root/scATAC/ATAC_data/Forebrain/Mulqueen_human_brain'
    path_graph_input = os.path.join(path_data_root, 'input_graph')
    dataset_atac_graph = ATACGraphDataset(path_graph_input)
    dataset = dataset_atac_graph.shuffle()
    train_dataset = dataset[:1700]
    test_dataset = dataset[1700:]
    device = torch.device("cuda:2" if torch.cuda.is_available() else "cpu")
    torch.manual_seed(123)

    model = Net().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.005)

    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)
    transform = T.Cartesian(cat=False)

    for epoch in range(1, 201):
        train(epoch)
        train_acc = test(train_loader)
        test_acc = test(test_loader)
        print(f'Epoch: {epoch:03d}, Train Acc: {train_acc:.4f}, Test Acc: {test_acc:.4f}')

    optimizer = torch.optim.Adam(model.parameters(), lr=0.003)
    for epoch in range(1, 101):
        train(epoch)
        train_acc = test(train_loader)
        test_acc = test(test_loader)
        print(f'Epoch: {epoch:03d}, Train Acc: {train_acc:.4f}, Test Acc: {test_acc:.4f}')

    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    for epoch in range(1, 101):
        train(epoch)
        train_acc = test(train_loader)
        test_acc = test(test_loader)
        print(f'Epoch: {epoch:03d}, Train Acc: {train_acc:.4f}, Test Acc: {test_acc:.4f}')

    optimizer = torch.optim.Adam(model.parameters(), lr=0.0005)
    for epoch in range(1, 101):
        train(epoch)
        train_acc = test(train_loader)
        test_acc = test(test_loader)
        print(f'Epoch: {epoch:03d}, Train Acc: {train_acc:.4f}, Test Acc: {test_acc:.4f}')

    time_end = time()
    print(time_end - time_start)
