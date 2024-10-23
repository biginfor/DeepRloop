import torch
from torch import nn
import torch.nn.functional as F
import numpy as np
import pandas as pd
__all__ = ['BassetNet2D', 'BassetNet1D']

device = torch.device("cpu") if not torch.cuda.is_available() \
    else torch.device("cuda:0")

class BassetNet2D(nn.Module):

    def __init__(self, dropout=0.3, num_targets=164, seq_len=600, KS1_len=19, Conv1_nus=300):
        super(BassetNet2D, self).__init__()
        #conv layer1
        self.conv1 = nn.Conv2d(in_channels=4,
                               out_channels=Conv1_nus,
                               kernel_size=(1, KS1_len),
                               stride=1,
                               padding=(0, (KS1_len-1)))
        self.batchnorm2d1 = nn.BatchNorm2d(Conv1_nus)
        self.relu1 = nn.ReLU()
        self.pool1 = nn.MaxPool2d(kernel_size=(1,3), ceil_mode=False)
        #conv layer2
        self.conv2 = nn.Conv2d(in_channels=Conv1_nus,
                               out_channels=200,
                               kernel_size=(1, 11),
                               stride=1,
                               padding=(0,5))
        self.batchnorm2d2 = nn.BatchNorm2d(200)
        self.relu2 = nn.ReLU()
        self.pool2 = nn.MaxPool2d(kernel_size=(1,4),ceil_mode=True)
        #conv layer3
        self.conv3 = nn.Conv2d(in_channels=200,
                               out_channels=200,
                               kernel_size=(1, 7),
                               stride=1,
                               padding=(0,3))
        self.batchnorm2d3 = nn.BatchNorm2d(200)
        self.relu3 = nn.ReLU()
        self.pool3 = nn.MaxPool2d(kernel_size=(1,4),ceil_mode=True)
        #flatten layer1
        self.flatten = nn.Flatten()
        self.fc1_in_dim = np.ceil((((seq_len - 3) / 3 + 1 - 4) / 4 + 1 - 4) / 4 + 1).astype(int)
        self.fc1 = nn.Linear(in_features=200*1*self.fc1_in_dim, out_features=1000)
        self.batchnormlinear1 = nn.BatchNorm1d(1000)
        self.relu4 = nn.ReLU()
        self.dropout1 = nn.Dropout(dropout)
        #flatten layer2
        self.fc2 = nn.Linear(1000,1000)
        self.batchnormlinear2 = nn.BatchNorm1d(1000)
        self.relu5 = nn.ReLU()
        self.dropout2 = nn.Dropout(dropout)
        #flatten layer3
        self.fc3 = nn.Linear(1000,num_targets)
        self.sigmoid = nn.Sigmoid()

        self.init_weight()

    def init_weight(self):
        for m in self.modules():
            if isinstance(m, nn.Conv2d):
                m.weight.data.normal_(0, 1)
            elif isinstance(m, nn.BatchNorm2d):
                m.weight.data.fill_(1)
                m.bias.data.zero_()

    def forward(self, x):
        x6 = self.conv1(x)
        x  = self.batchnorm2d1(x6)
        x1 = self.relu1(x)
        x2 = self.pool1(x1)
        x  = self.conv2(x2)
        x  = self.batchnorm2d2(x)
        x  = self.relu2(x)
        x3 = self.pool2(x)
        x  = self.conv3(x3)
        x  = self.batchnorm2d3(x)
        x = F.relu(x)
        x4 = self.pool3(x)
        x5 = torch.flatten(x4, 1)
        x = self.fc1(x5)
        x = self.batchnormlinear1(x)
        x = self.relu4(x)
        x = self.dropout1(x)
        x = self.fc2(x)
        x = self.batchnormlinear2(x)
        x = self.relu5(x)
        x = self.dropout2(x)
        x = self.fc3(x)
        x = self.sigmoid(x)
        return x, x1, x2, x3, x4, x5 ,x6


class BassetNet1D(nn.Module):

    def __init__(self, dropout=0.3, num_targets=164, seq_len=600, KS1_len=19, Conv1_nus=300):
        print("KS1_len: {}".format(KS1_len))
        super(BassetNet1D, self).__init__()
        #conv layer1
        self.conv1 = nn.Conv1d(in_channels=4,
                               out_channels=Conv1_nus,
                               kernel_size=KS1_len,
                               stride=1,
                               padding=(KS1_len-1)//2
                               )
        self.batchnorm2d1 = nn.BatchNorm1d(Conv1_nus)
        self.relu1 = nn.ReLU()
        self.pool1 = nn.MaxPool1d(kernel_size=3, ceil_mode=True)
        #conv layer2
        self.conv2 = nn.Conv1d(in_channels=Conv1_nus,
                               out_channels=200,
                               kernel_size=11,
                               stride=1,
                               padding=5)
        self.batchnorm2d2 = nn.BatchNorm1d(200)
        self.relu2 = nn.ReLU()
        self.pool2 = nn.MaxPool1d(kernel_size=4, ceil_mode=True)
        #conv layer3
        self.conv3 = nn.Conv1d(in_channels=200,
                               out_channels=200,
                               kernel_size=7,
                               stride=1,
                               padding=3)
        self.batchnorm2d3 = nn.BatchNorm1d(200)
        self.relu3 = nn.ReLU()
        self.pool3 = nn.MaxPool1d(kernel_size=4,ceil_mode=True)
        #flatten layer1
        self.flatten = nn.Flatten()
        self.fc1_in_dim = np.ceil((((seq_len-3)/3+1-4)/4+1-4)/4+1).astype(int)

        self.fc1 = nn.Linear(in_features=200*1*self.fc1_in_dim, out_features=1000)
        self.batchnormlinear1 = nn.BatchNorm1d(1000)
        self.relu4 = nn.ReLU()
        self.dropout1 = nn.Dropout(dropout)
        #flatten layer2
        self.fc2 = nn.Linear(1000,1000)
        self.batchnormlinear2 = nn.BatchNorm1d(1000)
        self.relu5 = nn.ReLU()
        self.dropout2 = nn.Dropout(dropout)
        self.fc3 = nn.Linear(1000, num_targets)

        self.sigmoid = nn.Sigmoid()

        self.init_weight()

    def init_weight(self):
        for m in self.modules():
            if isinstance(m, nn.Conv1d):
                m.weight.data.normal_(0, 1)
            elif isinstance(m, nn.BatchNorm1d):
                m.weight.data.fill_(1)
                m.bias.data.zero_()

    def forward(self, x):
            x6 = self.conv1(x)
            x = self.batchnorm2d1(x6)
            x1 = self.relu1(x)
            x2 = self.pool1(x1)
            x = self.conv2(x2)
            x = self.batchnorm2d2(x)
            x = self.relu2(x)
            x3 = self.pool2(x)
            x = self.conv3(x3)
            x = self.batchnorm2d3(x)
            x = F.relu(x)
            x4 = self.pool3(x)
            x5 = torch.flatten(x4, 1)
            x = self.fc1(x5)
            x = self.batchnormlinear1(x)
            x = self.relu4(x)
            x = self.dropout1(x)
            x = self.fc2(x)
            x = self.batchnormlinear2(x)
            x = self.relu5(x)
            x = self.dropout2(x)
            x7 = self.fc3(x)
            x = self.sigmoid(x7)
            return x, x1, x2, x3, x4, x5 ,x6, x7


if __name__ == "__main__":
    x0 = torch.randn((128, 4, 1, 600))
    x1 = torch.randn((128, 4, 600))
    model0 = BassetNet2D()
    out0 = model0(x0)
    model1 = BassetNet1D()