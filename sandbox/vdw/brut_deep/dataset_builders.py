
import torch
from torch.utils.data import Dataset
import pandas as pd


class ScaleCellDataset(Dataset):

    def __init__(self, xa, xb, ta):
        """
        this
        :param nr_experiments:
        """
        x = []
        y = []
        for i in range(len(xa)):
            inp = [xa[i][0][0], xa[i][0][1], xa[i][1][0], xa[i][1][1],
                   xb[i][0][0], xb[i][0][1], xb[i][1][0], xb[i][1][1]]
            # print(inp)
            x.append(inp)
            y.append([ta[i][0][0], ta[i][0][1], ta[i][1][0], ta[i][1][1]])

        self.x_train = torch.tensor(x, dtype=torch.float32)
        # print("xt", self.x_train[0])
        self.y_train = torch.tensor(y, dtype=torch.float32)

    def __len__(self):
        return len(self.x_train)

    def __getitem__(self, idx):
        return self.x_train[idx], self.y_train[idx]
