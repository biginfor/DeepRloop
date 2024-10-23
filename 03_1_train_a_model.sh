#!/bin/bash
source /home/fangzj/Workdir/software/anaconda/etc/profile.d/conda.sh
conda activate pytorch
#config
#whether to run a quick check 
quick_check=${1:-no}
genomeloc="/home/fangzj/Workdir/Basset_Unblocked/data/genomes/hg38/"
all=`echo {0..100..5}`

#loop
for i in ${all}
do
    echo ${i}
    clean_data_loc="/home/fangzj/Workdir/Basset_Unblocked/data/cleandata/flank_seqs/all_classes_flank/both_ext/up_split/1e-90/200_600_${i}/"
	if [ ! -d ${clean_data_loc} ]; then
	    mkdir ${clean_data_loc}
	fi
cd  ${clean_data_loc}
        if [ ! -f "bin_200_peak_600_${i}.fa" ]; then
	bedtools getfasta -fi ${genomeloc}chroms/hg38.fa -bed ${clean_data_loc}bin_200_600_${i}.bed -s -fo bin_200_peak_600_${i}.fa 
	pid1=$!
	wait $pid1
	else
	echo "fasta has been created"
	fi
	#
	if [ ! -f "bin_200_peak_600_${i}_uniq_act.h5" ]; then
	/home/fangzj/Workdir/Basset_Unblocked/scripts/rerrangefiles/seq_hdf5.py -t 0.15 -v 0.15 ./bin_200_peak_600_${i}.fa ./bin_200_600_${i}_act.txt ./bin_200_peak_600_${i}_uniq_act.h5
	pid2=$!
	wait $pid2
	else
	echo "h5 has been created"
	fi
        if [[ "$quick_check" == "yes" ]]; then
        quick_check_dir=${clean_data_loc}quick_check
		if [ ! -d ${quick_check_dir} ]; then
		    mkdir ${quick_check_dir}
		fi
		if [  -f ${quick_check_dir}/db.sqlite3 ]; then
		    echo 'db existed! removing!'
		    rm ${quick_check_dir}/db.sqlite3
		fi
        nohup python /home/fangzj/Workdir/Basset_Unblocked/scripts/optuna/train_with_h5.py --input_data ./bin_200_peak_600_${i}_uniq_act.h5 --save_dir ${quick_check_dir} --num_targets 51 --n_worker 20 --epochs 20 --CV '' --e_iters 10   --KS1_len 12 --n_trials 1 > ${quick_check_dir}flank_${i}.log  2>&1 &
        pid3=$!
	wait $pid3
        else
	nohup python /home/fangzj/Workdir/Basset_Unblocked/scripts/optuna/train_with_h5.py --input_data ./bin_200_peak_600_${i}_uniq_act.h5 --save_dir ${clean_data_loc} --num_targets 51 --n_worker 20 --epochs 20 --CV '' --e_iters 10   --KS1_len 12 --n_trials 2 > ${clean_data_loc}flank_${i}.log  2>&1 &
        pid3=$!
	wait $pid3
	fi
done
