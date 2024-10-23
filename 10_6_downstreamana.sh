#!/bin/bash
source /home/fangzj/Workdir/software/anaconda/etc/profile.d/conda.sh
conda activate pytorch
scripts_loc='./'
inputdata=${1:-'/path/to/input.h5'}
model_path=${2:-'/path/to/best_model'}
save_dir=${model_path}'samples_1000/'
motif_script='./motifs_unblocked.py'
cell_act='./data/cell_activity.txt'
#check
cd ${model_path}
python ${scripts_loc}check.py --inputdata ${inputdata} --save_dir ${save_dir} --model_path ./best.pt --sample_size 1000 --batch_size 100 --n_worker 20 --num_targets 51 --dropout 0 --seq_len 600 --KS1_len 12 --Conv1_nus 300
pid1=$!
wait $pid1
#motif
nohup  ${motif_script} ${save_dir}reprs.hdf5 ${save_dir}test_and_weights.hdf5 -a 0.5 -m  ./TFs_RLBPs_withheader_IDs.meme -o ${save_dir}TF_RLBP_comb > ${save_dir}motifs_comb.log 2>&1 &
pid2=$!
wait $pid2
#motif null
python ${scripts_loc}basset_motifs_null.py --model_file ${model_path}best.pt --data_file ${inputdata} --out_file ${model_path}/model_out.h5 --batch_size 128 --dropout 0 --seq_len 600 --KS1_len 12 --Conv1_nus 300
pid3=$!
wait $pid3
#infl 
python ${scripts_loc}basset_motifs_infl.py -d ${model_path}/model_out.h5 -m ${save_dir}/TF_RLBP_comb/table.txt -s 2000 -b 500 -o ${save_dir}TF_RLBP_comb_inflout -t ${cell_act} --subset '' ${model_path}best.pt ${inputdata}
pid4=$!
wait $pid4
