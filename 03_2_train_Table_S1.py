import os
import h5py
import argparse
from tqdm import tqdm
from sklearn import metrics
import torch
from torch.optim.lr_scheduler import StepLR
from models.Basset_net import BassetNet1D, BassetNet2D
from models.datasets import BassetDatasetH5, DataLoader
from models.utils import get_loss, get_optim
from myutils.misc import show_cfgs, AverageMeter, draw_curve

import optuna

def get_args():
    parser = argparse.ArgumentParser("DNA Training")

    # path of dataset
    parser.add_argument('--input_data', default='./', type=str, help="Please set the name of dataset :all_dataset.h5")

    # path of model weight
    parser.add_argument('--save_dir', default='./weights_kfold_narrow_test_tmp', type=str)

    # device, cpu or gpu
    parser.add_argument('--device', default=0)

    # model dim
    parser.add_argument('--model_dim', default=1, type=int,
                        help='model 2D or 1D, default 1D')

    # well-trained model
    parser.add_argument('--resume', default='', type=str)

    # epoch
    parser.add_argument('--epochs', default=200, type=int)

    # num_targets
    parser.add_argument('--num_targets', default=53, type=int)

    # cpu workers to load data
    parser.add_argument('--n_worker', default=20, type=int)

    # early stopping patience
    parser.add_argument('--e_iters', default=10, type=int)

    # the kernel size of first conv layer
    parser.add_argument('--KS1_len', default=12, type=int)

    # number of optuna trails
    parser.add_argument('--n_trials', default=1, type=int)

    args = parser.parse_args()
    return args


def train_epoch(
    model: torch.nn.Module,
    train_loader,
    optim,
    loss_fnc,
    device,
    scheduler,
    dim,
    records,
    num_targets,
    Conv1_nus,
):
    # prepare
    model.train()
    model.to(device)
    losses, current, n_zeros, nloop = 0, 0, 0, 0

    tqdm_loader = tqdm(train_loader, ncols=100, desc="train",mininterval=10,leave=False)
        
    for idx, (X, y) in enumerate(tqdm_loader):
        optim.zero_grad()
        # forward 
        X, y = X.to(device), y.to(device)
        outs = model(X)
        output = outs[0]
        relu1  = outs[1]
        loss_init = loss_fnc(output, y)
        loss=loss_init
        print("init_loss:{}".format(loss_init))
        print("batch_loss:{}".format(loss))

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
        print("batch_AUC:{}".format(cur_auc))
        loss.backward()
        optim.step()
        losses += loss.item()
        current += cur_auc.item()

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
                "auc": current / nloop,
                "lr": optim.state_dict()['param_groups'][0]['lr']
            }
        )
    records['loss'].update(losses / nloop)
    records['auc'].update(current / nloop)
    records['zero'].update(n_zeros / nloop)
    scheduler.step()
    return model, scheduler, records


@torch.no_grad()
def valid_epoch(model, dataloader, loss_fnc, device, dim, records, num_targets,Conv1_nus):
    model.eval()
    model.to(device)

    losses, current, n_zeros, nloop = 0, 0, 0, 0
    tqdm_loader = tqdm(dataloader, ncols=120, desc='valid')

    for idx, (X, y) in enumerate(tqdm_loader):
        # forward 
        X, y = X.to(device), y.to(device)
        outs = model(X)
        output = outs[0]
        relu1  = outs[1]
        loss = loss_fnc(output, y)

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

        cur_auc = torch.mean(AUC)
        losses += loss.item()
        current += cur_auc.item()


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
    records['loss'].update(losses / nloop)
    records['auc'].update(current / nloop)
    records['zero'].update(n_zeros / nloop)
    
    return records


def train_valid(args, params,save_model):
    show_cfgs(args)

    # parameters
    input_data = args.input_data
    save_dir = args.save_dir
    epochs = args.epochs
    n_worker = args.n_worker
    model_dim = args.model_dim
    num_targets = args.num_targets
    early_stopping_iter = args.e_iters

    device = torch.device("cpu") if not torch.cuda.is_available() \
                                else torch.device("cuda:{}".format(args.device))

    loss_fnc = get_loss()
    with h5py.File(input_data, 'r') as dataset_1:
        train_in = dataset_1['train_in']
        valid_in = dataset_1['valid_in']
        n_samples_train = train_in.shape[0]
        n_samples_valid = valid_in.shape[0]
        seq_len = train_in.shape[-1]
    print("seqs:{}".format(seq_len))
    train_loader = DataLoader(
        dataset=BassetDatasetH5(
            filename=input_data,
            n_samples=n_samples_train, dim=model_dim, phase='train'),
            num_workers=n_worker, batch_size=params["batch_size"],
            pin_memory=False, shuffle=True,
            drop_last=True)
    valid_loader = DataLoader(
        dataset=BassetDatasetH5(
            filename=input_data,
            n_samples=n_samples_valid, dim=model_dim, phase='valid'),
            num_workers=n_worker, batch_size=params["batch_size"],
            pin_memory=False, shuffle=False, drop_last=False)

    train_rec = {
        "loss": AverageMeter(),
        "auc":  AverageMeter(),
        "zero": AverageMeter(),
    }
    valid_rec = {
        "loss": AverageMeter(),
        "auc":  AverageMeter(),
        "zero": AverageMeter(),
    }

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    weight_files = os.listdir(save_dir)
    types = 'm{}d_opt{}_drop{:.2f}_lr{:.5f}_step{}_gama{:.2f}_mom{:.2f}_batch{}_ks1_{}_conv1_{}_norm_{}'.format(str(model_dim),
           params["optimizer_name"],
           params["dropout"],
           params["learning_rate"],
           params["lr_step"],
           params["lr_gamma"],
           params["momentum"],
           params["batch_size"],
           params["KS1_len"],
           params["Conv1_nus"]
           )
    if len(weight_files) <= 0:
        save_dir = os.path.join(save_dir, '{}_run_0'.format(types))
    else:
        num = [int(f.split('_')[-1]) for f in weight_files if f.find(types) >= 0]
        num = [-1] if len(num) <= 0 else num
        save_dir = os.path.join(save_dir, '{}_run_{}'.format(types, max(num) + 1))

    print("optuna begin")
    model = BassetNet1D(dropout=params["dropout"],
                        num_targets=num_targets,seq_len=seq_len, KS1_len=params["KS1_len"],Conv1_nus=params["Conv1_nus"]) if model_dim == 1 else BassetNet2D(
        params["dropout"], num_targets=num_targets,seq_len=seq_len, KS1_len=params["KS1_len"],Conv1_nus=params["Conv1_nus"])
    if args.resume:
        state = torch.load(args.resume, map_location='cpu')
        model.load_state_dict(state)

    optim = get_optim(model, name=params["optimizer_name"],
                      lr=params["learning_rate"],
                      momentum=params["momentum"])
    scheduler = StepLR(optim, step_size=params["lr_step"],
                       gamma=params["lr_gamma"])
    early_stopping_counter = 0
    min_auc = 0
    #
    for epoch in range(epochs):
        print(" ======== Epoch: {} ======== ".format(epoch))

        model, scheduler, train_rec = train_epoch(
            model, train_loader, optim, loss_fnc, device,
            scheduler, model_dim, train_rec, num_targets,Conv1_nus=params["Conv1_nus"]
        )

        train_loss = train_rec['loss'].avg
        train_auc = train_rec['auc'].avg
        train_rec['loss'].clear()
        train_rec['auc'].clear()
        train_rec['zero'].clear()

        valid_rec = valid_epoch(
            model, valid_loader, loss_fnc, device, model_dim, valid_rec,
            num_targets,Conv1_nus=params["Conv1_nus"]
        )

        val_loss = valid_rec['loss'].avg
        val_auc = valid_rec['auc'].avg
        valid_rec['loss'].clear()
        valid_rec['auc'].clear()
        valid_rec['zero'].clear()
        if not os.path.exists("{}".format(save_dir)):
            os.makedirs("{}".format(save_dir))
        print("auc_min: {}".format(min_auc))
        print("auc_val: {}".format(val_auc))

        # loss
        draw_curve(
            x_train=train_rec['loss'].lst,
            x_valid=valid_rec['loss'].lst,
            title="loss curve",
            y_label='',
            save_path=os.path.join("{}".format(save_dir),
                                   'loss.png'),
            legents=['train_loss', 'valid_loss'])

        # auc
        draw_curve(
            x_train=train_rec['auc'].lst,
            x_valid=valid_rec['auc'].lst,
            title="auc curve",
            y_label='',
            save_path=os.path.join("{}".format(save_dir),
                                   'auc.png'),
            legents=['train_auc', 'valid_auc'])

        # zeros
        draw_curve(
            x_train=train_rec['zero'].lst,
            x_valid=valid_rec['zero'].lst,
            title="zero curve",
            y_label='',
            save_path=os.path.join("{}".format(save_dir),
                                   'zeros.png'),
            legents=['train_zeros', 'valid_zeros'])
        # save the best model
        print("current_val_auc:{}\ncurrent_min_auc:{}".format(val_auc,min_auc))
        if val_auc > min_auc:
            min_auc = val_auc
            if save_model:
                save_to = os.path.join("{}".format(save_dir),
                                       'best.pt')
                torch.save(model.state_dict(), save_to)
                print("模型保存至 '{}'.".format(save_to))
        # early stopping
        else:
            early_stopping_counter += 1
        print("early_stopping_iter:{}".format(early_stopping_iter))
        print("early_stopping_counter:{}".format(early_stopping_counter))
        if early_stopping_counter > early_stopping_iter:
            print("auc couldn't be better")
            break
        print("第{}个epoch内的min_auc:{}".format(epoch, min_auc))
    print("current trial best auc: {}".format(min_auc))
    return min_auc

def objective(trial,args):
    params = {
        "optimizer_name": trial.suggest_categorical("optimizer_name",
                                                    ["rmsprop"]),
        "dropout": trial.suggest_float("dropout", 0.3885, 0.3885), #0.25-0.60
        "learning_rate": trial.suggest_float("learning_rate", 2.4637e-05, 2.4637e-05, log=True), #5e-06 - 1e-04
        "lr_step": trial.suggest_int("lr_step", 14, 14), # 5 - 15
        "lr_gamma": trial.suggest_float("lr_gamma", 0.2845, 0.2845), # 0.15 - 0.35
        "momentum": trial.suggest_float("momentum", 0.9711, 0.9711), # 0.90 - 0.99
        "batch_size": trial.suggest_int("batch_size", 512, 512), 
        "KS1_len":trial.suggest_int("KS1_len", args.KS1_len, args.KS1_len), # 2 - 50
        "Conv1_nus":trial.suggest_int(name="Conv1_nus", low=300, high=300, step=50)
    }
    accuracy = train_valid(get_args(), params, save_model=True)
    return accuracy

if __name__ == "__main__":
    if not os.path.exists(get_args().save_dir):
        os.makedirs(get_args().save_dir)
    study = optuna.create_study(study_name="test", direction="maximize",storage='sqlite:////{}/db.sqlite3'.format(get_args().save_dir))
    study.optimize(lambda trial: objective(trial, get_args()), get_args().n_trials)
    print("best trial:")
    best_trial = study.best_trial
    print(best_trial.values)
    print(best_trial.params)
    scr = train_valid(get_args(), best_trial.params, save_model=True)
    print("The best valid auc is {}".format(scr))
