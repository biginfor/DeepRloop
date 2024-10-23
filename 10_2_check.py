import os
import h5py
import numpy as np
import argparse
from tqdm import tqdm
from sklearn import metrics
import torch
from models.Basset_net import BassetNet1D, BassetNet2D
from models.datasets import BassetDatasetH5, DataLoader
from models.utils import get_loss
from myutils.misc import show_cfgs


def get_args():

    parser = argparse.ArgumentParser("DNA Testing")

    # path of data
    parser.add_argument('--inputdata', default='datas', type=str)

    # path of model weight
    parser.add_argument('--save_dir', default='outputs', type=str)

    # device, cpu or gpu
    parser.add_argument('--device', default=0)

    # model dim
    parser.add_argument('--model_dim', default=1, type=int, 
                        help='model 2D or 1D, default 1D')

    # best model parameters
    parser.add_argument('--model_path', default='', type=str)

    # batch size
    parser.add_argument('--batch_size', default=1000, type=int)

    # sample size
    parser.add_argument('--sample_size', default=10000, type=int)

    # threads
    parser.add_argument('--n_worker', default=5, type=int)

    # num_targets
    parser.add_argument('--num_targets', default=26, type=int)

    # dropout
    parser.add_argument('--dropout', default=0.3885, type=float)

    #
    parser.add_argument('--seq_len', default=600, type=int)

    #
    parser.add_argument('--KS1_len', default=19, type=int)
    
    #
    parser.add_argument('--Conv1_nus', default=300, type=int)

    args = parser.parse_args()
    return args


@torch.no_grad()
def test_epoch(model, dataloader, loss_fnc, device, dim, num_targets,save_dir=None,Conv1_nus=300,sample_size=10000,seq_len=600):
    model.eval()
    model.to(device)

    losses, current, n_zeros, nloop = 0, 0, 0, 0
    pi = 0
    mi = 0
    init_reprs = True
    preds  = torch.FloatTensor(len(dataloader.dataset), num_targets)
    labels = torch.FloatTensor(len(dataloader.dataset), num_targets)
    tqdm_loader = tqdm(dataloader, ncols=120, desc='test')
    count = 0
    #
    for idx, (X, y) in enumerate(tqdm_loader):
        count += 1
        # forward 
        X, y = X.to(device), y.to(device)
        print("The input data shape is {}".format(X.shape))
        outs = model(X)
        output = outs[0]
        relu1  = outs[1]
        loss = loss_fnc(output, y)

        for i in range(X.shape[0]):
            preds[pi, :] = output[i, :].float()
            labels[pi, :] = y[i, :].float()
            pi += 1

        # initialize reprs
        if init_reprs:
            reprs = torch.FloatTensor(sample_size,
                                      relu1.shape[1],
                                      relu1.shape[2])
        # initialize test seqs
            test_seqs = torch.FloatTensor(sample_size,
                                      4,
                                      seq_len)
        
        # larger tensor
        for i in range(relu1.shape[0]):
            reprs[mi,:,:] = relu1[i,:,:].float()
            test_seqs[mi,:,:] = X[i,:,:]
            mi += 1
        init_reprs = False
        # compute auc
        AUC = torch.Tensor(num_targets)     
        for yi in range(num_targets):
            vec1 = y[:, yi]
            vec2 = output[:, yi]
            if vec1.sum()== 0:
                AUC[yi] = 1.0
            else:
                fpr, tpr, thr = metrics.roc_curve(
                    y_true=vec1.detach().cpu().numpy(),
                    y_score=vec2.detach().cpu().numpy(),
                    pos_label=1)
                AUC[yi] = metrics.auc(fpr, tpr)
         # backward
        cur_auc = torch.mean(AUC)
        losses += loss.item()
        current += cur_auc.item()

        # number of zeros
        if dim == 2:
            relu1 = relu1.squeeze(2)
        for i in range(Conv1_nus):
            if torch.sum(relu1[:, i, :] != 0) == 0:
                n_zeros += 1
        
        nloop += 1

        tqdm_loader.set_postfix(
            {
                "loss": losses / nloop,
                "zeros": n_zeros / nloop,
                "auc": current / nloop
            }
        )

    total_zeros = 0
    lst = []
    for i in range(Conv1_nus):
        nnn = torch.sum(reprs[:, i, :] != 0)
        if nnn != 0:
            total_zeros += 1
        lst.append(nnn.item())
    print("non zeros in {} channels: {}".format(Conv1_nus,total_zeros))
    print("max z in {} channels: {}, min z in {} channels: {}".format(Conv1_nus,int(max(lst)), Conv1_nus,int(min(lst))))

    no_0 = []
    filter_out = reprs.numpy()
    for i in range(Conv1_nus):
        if np.unique(filter_out[:, i, :]).shape[0] > 1:
            no_0.append(1)
        else:
            no_0.append(0)
    print("\nno_0: {}".format(no_0))
    print("no_0 == 1: ", np.count_nonzero(no_0 == 1))

    with h5py.File(os.path.join(save_dir, "reprs.hdf5"), "w") as hdf_out:
        hdf_out.create_dataset('outs', data=reprs)
    hdf_out.close()

    with h5py.File(os.path.join(save_dir, "test_and_weights.hdf5"), 'w') as f_out:
        X_tensor = test_seqs.cpu()
        weights_conv = model.conv1.weight.data.cpu()
        f_out.create_dataset('weights', data=weights_conv)
        f_out.create_dataset('test_in', data=X_tensor)
    f_out.close()


def main(args):
    
    show_cfgs(args)

    # parameters
    inputdata = args.inputdata
    save_dir = args.save_dir
    sample_size =args.sample_size
    batch_size = args.batch_size
    n_worker = args.n_worker
    model_dim = args.model_dim
    num_targets =args.num_targets
    seq_len=args.seq_len
    KS1_len=args.KS1_len
    Conv1_nus=args.Conv1_nus
    # device
    device = torch.device("cpu")  if not torch.cuda.is_available() \
                                else torch.device("cuda:{}".format(args.device))
    
    # loss
    loss_fnc = get_loss()

    # model
    model = BassetNet1D(dropout=args.dropout,
                        num_targets=num_targets,seq_len=seq_len,KS1_len=KS1_len,Conv1_nus=Conv1_nus) if model_dim == 1 else BassetNet2D(
        args.dropout, num_targets=args.num_targets,seq_len=seq_len,KS1_len=KS1_len,Conv1_nus=Conv1_nus)
    state = torch.load(args.model_path, map_location='cpu')
    model.load_state_dict(state)

    # dataloader
    test_loader = DataLoader(
        dataset=BassetDatasetH5(filename=inputdata,
                                n_samples=args.sample_size, dim=model_dim, phase='test'),
        num_workers=n_worker, batch_size=batch_size,
        pin_memory=False, shuffle=False, drop_last=False)
    
    named = os.path.dirname(args.model_path).split('/')[-1]
    save_root = os.path.join(save_dir, named)
    if not os.path.exists(save_root):
        os.makedirs(save_root)
    test_epoch(model, test_loader, loss_fnc, device, model_dim, num_targets, save_root,Conv1_nus,sample_size,seq_len)


if __name__ == "__main__":
    main(get_args())