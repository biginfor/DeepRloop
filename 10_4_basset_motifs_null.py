#ref https://github.com/davek44/Basset
import h5py
import argparse
from tqdm import tqdm
import torch
from myutils.misc import show_cfgs
from models.Basset_net import BassetNet1D
from models.datasets import BassetDatasetH5, DataLoader
from models.utils import get_loss


def get_args():

    parser = argparse.ArgumentParser("DNA ConvNet hidden layer visualizations")

    # the well-trained model
    parser.add_argument('--model_file', default='/path/to/model/best.pt', type=str, help="the pre-trained model")
    
    # path of input file
    parser.add_argument('--data_file', default='/path/to/all_dataset.h5', type=str,help="the origin data file")

    # path of output hdf5 file
    parser.add_argument('--out_file', default='/path/to/model_out.h5', type=str,help="out file name ")

    # batch
    parser.add_argument('--batch_size', default=100, type=int)

    # dropout
    parser.add_argument('--dropout', default=0.30, type=float)

    #
    parser.add_argument('--seq_len', default=600, type=int)

    #
    parser.add_argument('--KS1_len', default=19, type=int)
    
    #
    parser.add_argument('--Conv1_nus', default=300, type=int)

    args = parser.parse_args()
    return args
    
    # fix seed
    torch.manual_seed(1)


@torch.no_grad()
def test_epoch(test_seqs,test_targets,model, dataloader, loss_fnc, device, out_file,num_filters):
    model.eval()
    model.to(device)
    # define some params
    num_seqs = test_seqs.shape[0]
    num_targets = test_targets.shape[1]

    filter_means = torch.zeros(num_filters) # 300 means
    filter_stds = torch.zeros(num_filters)
    filter_infl = torch.zeros(num_filters)
    filter_infl_targets = torch.zeros(num_filters, num_targets) # torch.Size([300, 164])
    
    # test params
    batches = 1
    tqdm_loader = tqdm(dataloader, ncols=120, desc='test') #loading bar 

    for idx, (X, y) in enumerate(tqdm_loader):
        # forward 
        X, y = X.to(device), y.to(device)
        outs = model(X)
        relu1  = outs[1]
        fc3 = outs[7]

        preds = outs[0]
        loss = loss_fnc(preds, y)
        final_output = fc3.clone().detach()
        final_means = torch.mean(final_output,dim=0).reshape(num_targets)
        final_stds = torch.std(final_output,dim=0).reshape(num_targets)
        final_var = torch.pow(final_stds,2)
        
        #save original
        actual_rep = relu1.squeeze().clone()
        
        batch_size = actual_rep.shape[0] # 100
        seq_len = actual_rep.shape[2] # 600
        # compute filter means
        filter_means_batch = torch.mean(actual_rep,dim=(0,2)) # [300]
        filter_means =filter_means.cuda() + filter_means_batch.cuda()
        
        # compute filter stds
        filters_out1 = actual_rep.permute(1,0,2).reshape(num_filters,batch_size*seq_len) # swap the dim0  and dim1 [100,300,600]>[300,100,600]>[300,100*600]
        filter_stds = filter_stds.cuda() + filters_out1.std(dim=1).squeeze().cuda()  # cal the std along row [300]
        # nullify and re-measure
        for fi in range(num_filters):
            print(" Filter {}".format(fi))
            zloss = loss
            zfinal_means = final_means
            # if the unit is inactive
            if filter_means_batch[fi] == 0:
                print("Inactivte")
            else:
                # zero the hidden unit
                zero_rep = actual_rep.clone()
                zero_rep[:, fi, :] = filter_means_batch[fi]
                relu1_null = zero_rep.reshape(batch_size, num_filters, seq_len) # relu1_null needs to replace origin layer

                # propagate the change through the network
                def get_output():
                    for name, module in model.named_modules():
                        if name in ['', 'conv1', 'batchnorm2d1']:
                            #print("init layers")
                            pass
                        elif name.startswith("relu1"):
                            input_init = relu1_null
                            #print(input_init.shape)
                        else:
                            #print(name)
                            input_init = module(input_init)
                            if name == 'fc3':
                                    zfinal_output = input_init
                            if name == 'sigmoid':
                                return  input_init, zfinal_output
                zpreds, zfinal_output = get_output()
                zfinal_means = torch.mean(zfinal_output, dim=0).reshape(num_targets)
                # re-compute loss
                print(torch.allclose(preds, zpreds))
                zloss = loss_fnc(zpreds, y)

        # z normalized difference
            batch_infl_targets = final_means-zfinal_means
            batch_infl_targets.div(final_var)
            filter_infl_targets[fi] = filter_infl_targets[fi].cuda() + batch_infl_targets.cuda()
        
            # save unit delta
            filter_infl[fi] = filter_infl[fi] + (zloss - loss)
            print(filter_infl[fi])
        batches += 1
        print(" one batch ")
    batches = batches -1

    #norm by batch number
    filter_means = filter_means / batches
    filter_stds = filter_stds /batches
    filter_infl_targets = filter_infl_targets /batches

    filter_infl = filter_infl / num_seqs
    with h5py.File(out_file, 'w') as hdf_out:
        hdf_out.create_dataset('filter_means', data=filter_means.cpu())
        hdf_out.create_dataset('filter_stds', data=filter_stds.cpu())
        hdf_out.create_dataset('filter_infl', data=filter_infl.cpu())
        hdf_out.create_dataset('filter_infl_targets', data=filter_infl_targets.cpu())
    hdf_out.close()

def main(args):
    show_cfgs(args)
    data_file = args.data_file
    num_filters = args.Conv1_nus
    device = torch.device("cuda:0")
    with h5py.File(data_file, 'r') as data_open:
        test_seqs = data_open['test_in'][:]
        test_targets = data_open['test_out'][:]

    sample_size = args.batch_size
    n_worker = 20

    # loss
    loss_fnc = get_loss()

    # model
    model = BassetNet1D(dropout=args.dropout, num_targets=test_targets.shape[1],seq_len=args.seq_len,KS1_len=args.KS1_len,Conv1_nus=num_filters)
    state = torch.load(args.model_file, map_location='cuda:0')
    model.load_state_dict(state)

    # dataloader
    test_loader = DataLoader(
        dataset=BassetDatasetH5(filename=data_file,
                                n_samples=18000, dim=1, phase='test'),
        num_workers=n_worker, batch_size=sample_size,
        pin_memory=False, shuffle=False, drop_last=False)
    
    out_file = args.out_file
    test_epoch(test_seqs,test_targets,model, test_loader, loss_fnc, device, out_file,num_filters)


if __name__ == "__main__":
    main(get_args())
