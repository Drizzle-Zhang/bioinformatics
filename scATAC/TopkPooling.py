# _*_ coding: utf-8 _*_
# @author: Drizzle_Zhang
# @file: TopkPooling.py
# @time: 2022/1/5 10:08

import torch
import torch.nn.functional as F
import os
from torch_geometric.data import InMemoryDataset, DataLoader
from torch_geometric.nn import GraphConv, TopKPooling
from torch_geometric.nn import global_mean_pool as gap, global_max_pool as gmp


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


path_data_root = '/root/scATAC/ATAC_data/Forebrain/Mulqueen_human_brain'
path_graph_input = os.path.join(path_data_root, 'input_graph')
dataset_atac_graph = ATACGraphDataset(path_graph_input)
torch.manual_seed(12345)
dataset = dataset_atac_graph.shuffle()
train_dataset = dataset[:1700]
test_dataset = dataset[1700:]


class Net(torch.nn.Module):
    def __init__(self):
        super().__init__()

        self.conv1 = GraphConv(dataset.num_features, 32)
        self.pool1 = TopKPooling(32, ratio=0.3)
        self.conv2 = GraphConv(32, 32)
        self.pool2 = TopKPooling(32, ratio=0.3)
        self.conv3 = GraphConv(32, 32)
        self.pool3 = TopKPooling(32, ratio=0.3)

        self.lin1 = torch.nn.Linear(64, 32)
        # self.lin2 = torch.nn.Linear(128, 64)
        self.lin3 = torch.nn.Linear(32, dataset.num_classes)

    def forward(self, data):
        x, edge_index, edge_weight, batch = data.x, data.edge_index, data.edge_attr, data.batch

        x = F.relu(self.conv1(x, edge_index))
        x, edge_index, edge_weight, batch, _, _ = self.pool1(x, edge_index, edge_weight, batch)
        x1 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)

        x = F.relu(self.conv2(x, edge_index))
        x, edge_index, edge_weight, batch, _, _ = self.pool2(x, edge_index, edge_weight, batch)
        x2 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)

        x = F.relu(self.conv3(x, edge_index))
        x, edge_index, edge_weight, batch, _, _ = self.pool3(x, edge_index, edge_weight, batch)
        x3 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)

        x = x1 + x2 + x3

        x = F.relu(self.lin1(x))
        x = F.dropout(x, p=0.5, training=self.training)
        # x = F.relu(self.lin2(x))
        x = F.log_softmax(self.lin3(x), dim=-1)

        return x


device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")
model = Net().to(device)
train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)

optimizer = torch.optim.Adam(model.parameters(), lr=0.005)


def train(epoch):
    model.train()

    loss_all = 0
    for data in train_loader:
        data = data.to(device)
        optimizer.zero_grad()
        output = model(data)
        loss = F.nll_loss(output, data.y)
        loss.backward()
        loss_all += data.num_graphs * loss.item()
        optimizer.step()
    return loss_all / len(train_dataset)


def test(loader):
    model.eval()

    correct = 0
    for data in loader:
        data = data.to(device)
        pred = model(data).max(dim=1)[1]
        correct += pred.eq(data.y).sum().item()
    return correct / len(loader.dataset)


for epoch in range(1, 201):
    loss = train(epoch)
    train_acc = test(train_loader)
    test_acc = test(test_loader)
    print(f'Epoch: {epoch:03d}, Loss: {loss:.5f}, Train Acc: {train_acc:.5f}, '
          f'Test Acc: {test_acc:.5f}')
