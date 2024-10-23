#!/bin/bash

merge_overlap=200
flank_len=`echo {0..100..5}`
ext_len=600
flank_ratio=1e-90
direc_for_ext=('up' 'both' 'dn')
direc_for_split=('up' 'both' 'dn')

for i in "${direc_for_ext[@]}"
do
    for j in "${direc_for_split[@]}"
    do
        output_dir="./data/all_classes_flank/${i}_ext/${j}_split/${flank_ratio}"
        if [ ! -d  ${output_dir} ]; then
            mkdir -p ${output_dir}
        else
            echo "${output_dir} already exists"
        fi
        cd ${output_dir}
            for k in ${flank_len}
            do
            prefix=bin_${merge_overlap}_${ext_len}_${k}
            echo flank_len_is:${k}
                if [ -f ${prefix}.bed ] && [ -f ${prefix}_act.txt ]; then
                    echo "file exists! next!"  
                else
                    python 02_2_peak_merging_splitting.py -c ./data/hg38.chrom.sizes -m ${merge_overlap} -o ${prefix} -s ${ext_len} -f ${k} -y -r ${flank_ratio} -d ${i},${j} ./data/Rloop_E_R_meta.txt
                fi
            done
    done
done



