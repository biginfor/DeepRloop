import h5py
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader

# A:1000,C:0100,G:0010,T:0001
class BassetDatasetH5(Dataset):

    def __init__(self, filename, n_samples=None, dim=2, phase='train'):
        super(BassetDatasetH5, self).__init__()
        assert dim in [1, 2], "模型仅支持2D, 1D"
        self.dim = dim
        self.len_data = n_samples
        self.phase = phase
        self.data_open = h5py.File(filename, 'r')
    
    def __len__(self):
        return self.len_data
    
    def __getitem__(self, index):
        datas = self.data_open[f'{self.phase}_in'][index]
        targs = self.data_open[f'{self.phase}_out'][index]

        datas = np.expand_dims(datas, 0)
        targs = np.expand_dims(targs, 0)

        datas, targs = data_arrange(datas, targs)
        datas = torch.from_numpy(datas)
        targs = torch.from_numpy(targs)
        if self.dim == 1:
            datas = datas.squeeze(2)
        datas = datas[0]
        targs = targs[0]
        return datas, targs


class BassetDataset(Dataset):
    
    def __init__(self, data, label, n_samples=None, dim=2, phase='train'):
        super(BassetDataset, self).__init__()
        assert dim in [1, 2], "模型仅支持2D, 1D"
        self.dim = dim
        self.len_data = n_samples
        self.x_data = torch.from_numpy(data)
        self.phase=phase
        if self.phase == "predict":
           pass
        else:
            self.y_data = torch.from_numpy(label)
    
    def __len__(self):
        if self.len_data is not None:
            return self.len_data
        return self.x_data.size(0)

    def __getitem__(self, idx):
        x = self.x_data[idx]
        if self.dim == 1:
            x = x.squeeze(1)
        if self.phase != "predict":
            y = self.y_data[idx]
            return x, y
        else:
            return x
    

def load_data_from_pkl(filename, phase='train'):
    """ 从pickle中加载数据 """
    import pickle as pkl
    with open(filename, 'rb') as f:
        data_dict = pkl.load(f)
    
    return data_dict[f'{phase}_in'], data_dict[f'{phase}_out'], data_dict['total_samples']


def load_data_from_h5py(filename, phase='test', max_size=-1):
    """ 从h5py中加载数据 """
    data_open = h5py.File(filename, 'r')
    if max_size <= 0:
        datas = data_open[f'{phase}_in'][:]
        targs = data_open[f'{phase}_out'][:]
    else:
        datas = data_open[f'{phase}_in'][:max_size]
        targs = data_open[f'{phase}_out'][:max_size]
    return datas, targs, data_open['target_labels']


def bool2int(matrix):
    return matrix.astype(np.int_)


def data_arrange(seqs, target):
    """ 数据整理 """
    data = bool2int(seqs)
    data = data.flatten()
    # one_hot_encoded1,70000行，600列，4面
    data = data.reshape(seqs.shape[0], seqs.shape[-1], 4)
    data = data.swapaxes(1, 2) # 列和面交换
    data = data[:, :, np.newaxis]
    data = data.astype(np.float32)
    target = target.astype(np.float32)
    return data, target
