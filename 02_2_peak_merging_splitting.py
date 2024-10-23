#!/usr/bin/env python
#ref. https://github.com/davek44/Basset
from optparse import OptionParser
import gzip
import os
import subprocess
import sys
from collections import defaultdict
import h5py
import numpy as np
import pandas as pd
import shutil
import math

def main():
    usage = 'usage: %prog [options] <target_beds_file>'
    parser = OptionParser(usage)
    parser.add_option('-a', dest='db_act_file',help='Existing database activity table.')
    parser.add_option('-b', dest='db_bed', help='Existing database BED file.')
    parser.add_option('-c', dest='chrom_lengths_file', default="/home/fangzj/Workdir/Basset_Unblocked/data/genomes/hg38/hg38.chrom.sizes", help='Table of chromosome lengths')
    parser.add_option('-m', dest='merge_overlap', default=200, type='int', help='Overlap length (after extension to feature_size) above which to merge features [Default: %default]')
    parser.add_option('-n', dest='no_db_activity', default=False, action='store_true', help='Do not pass along the activities of the database sequences [Default: %default]')
    parser.add_option('-o', dest='out_prefix', default='features', help='Output file prefix [Default: %default]')
    parser.add_option('-s', dest='feature_size', default=600, type='int', help='Extend features to this size [Default: %default]')
    parser.add_option('-f', dest='flank_size', default=300, type='int', help='Flanking length of peak start and stop sites, feature_size [Default: %default]')
    parser.add_option('-y', dest='ignore_y', default=True, action='store_true', help='Ignore Y chromsosome features [Default: %default]')
    parser.add_option('-r', dest='flank_ratio', default=0.8, type='float', help='The ratio of flank seqs length to full peak length  [Default: %default]')
    parser.add_option('-d', dest='direction', default='both,up', type='string', help='flank to upstream, downstream, or both for ext and split,comma-split  [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide file labeling the targets and providing BED file paths.')
    else:
        target_beds_file = args[0]

    db_targets = []
    db_add = False
    if options.db_bed:
        db_add = True
        if not options.no_db_activity:
            if options.db_act_file is None:
                parser.error('Must provide both activity table or specify -n if you want to add to an existing database')
            else:
                db_act_in = open(options.db_act_file)
                db_targets = db_act_in.readline().strip().split('\t')
                db_act_in.close()


    target_beds = []
    target_dbi = []
    for line in open(target_beds_file):
        a = line.split()
        if len(a) != 2:
            print(a)
            print('Each row of the target BEDS file must contain a label and BED file separated by whitespace',file=sys.stderr)
            exit(1)
        target_dbi.append(len(db_targets))
        db_targets.append(a[0])
        target_beds.append(a[1])

    # read in chromosome lengths
    chrom_lengths = {}
    if options.chrom_lengths_file:
        chrom_lengths = {}
        for line in open(options.chrom_lengths_file):
            a = line.split()
            chrom_lengths[a[0]] = int(a[1])
    else:
        print('Warning: chromosome lengths not provided, so regions near ends may be incorrect.',file=sys.stderr)


    chrom_files = defaultdict(list)
    chrom_outs = {}

    peak_beds = target_beds
    if db_add:
        peak_beds.append(options.db_bed)

    for bi in range(len(peak_beds)):
        if peak_beds[bi][-3:] == '.gz':
            peak_bed_in = gzip.open(peak_beds[bi])
        else:
            peak_bed_in = open(peak_beds[bi])

        for line in peak_bed_in:
            a = line.split('\t')
            a[-1] = a[-1].rstrip()
            a = a[:7]

            chrom = a[0]
            strand = '+'
            if len(a) > 5 and a[5] in '+-':
                strand = a[5]
            chrom_key = (chrom, strand)

            if chrom_key not in chrom_outs:
                chrom_files[chrom_key].append('%s_%s_%s.bed' % (options.out_prefix, chrom, strand))
                chrom_files[chrom_key].append('%s_%s_%s_cut.bed' % (options.out_prefix, chrom, strand))
                chrom_outs[chrom_key] = open(chrom_files[chrom_key][1], 'w')


            if db_add and bi == len(peak_beds)-1:
                if options.no_db_activity:

                    a[6] = '.'
                    print('\t'.join(a[:7]), file = chrom_outs[chrom_key])
                else:
                    print(line, chrom_outs[chrom_key])


            else:
                # specify the target index
                while len(a) < 7:
                    a.append('')
                a[5] = strand
                a[6] = str(target_dbi[bi])

                tmp = a[:]

                direction_lst = options.direction
                direc_for_ext = direction_lst.split(',')[0]
                direc_for_split = direction_lst.split(',')[1]
                subpeaks = seperate(start=int(a[1]),
                                    end=int(a[2]),
                                    flank_len=options.flank_size,
                                    ext_len=options.feature_size,
                                    flank_ratio=options.flank_ratio,
                                    chrom_len=chrom_lengths.get(chrom, None),
                                    direc_for_ext=direc_for_ext,
                                    direc_for_split=direc_for_split)
                for (key, value) in subpeaks.items():
                    tmp.append(tmp[3]+"_"+key)
                    tmp.append(str(value[0]))
                    tmp.append(str(value[1]))
                    print('\t'.join(tmp[:]), file=chrom_outs[chrom_key])
                    tmp = a[:]
        peak_bed_in.close()


    for chrom_key in chrom_outs:
        chrom_outs[chrom_key].close()


    if options.ignore_y:
        for orient in '+-':
            chrom_key = ('chrY',orient)
            if chrom_key in chrom_files:
                os.remove(chrom_files[chrom_key][1])
                del chrom_files[chrom_key]



    for chrom_key in chrom_files:
        chrom, strand = chrom_key
        sort_dir = os.getcwd()+"/sort"
        if not os.path.exists(sort_dir):
            os.mkdir(sort_dir)
        chrom_sbed = '%s_%s_%s_sort.bed' % (sort_dir+"/"+options.out_prefix,chrom,strand)
        sort_cmd = 'sort -k 1,1 -k9,9n %s > %s' % (chrom_files[chrom_key][1], chrom_sbed)
        subprocess.call(sort_cmd, shell=True)
        os.remove(chrom_files[chrom_key][1])
        chrom_files[chrom_key] = chrom_sbed



    final_bed_out = open('%s.bed' % options.out_prefix, 'w')
    for chrom_key in chrom_files:
        chrom, strand = chrom_key

        open_peaks = []

        for line in open(chrom_files[chrom_key]):
            a = line.split('\t')
            a[-1] = a[-1].rstrip()

            peak_start = int(a[8])
            peak_end = int(a[9])
            peak_act = activity_set(a[6])
            peak = Peak(peak_start, peak_end, peak_act)

            if len(open_peaks) == 0:


                open_end = peak.end
                open_peaks = [peak]
            else:

                if open_end - options.merge_overlap <= peak.start:
                    mpeaks = merge_peaks(open_peaks, options.feature_size, options.merge_overlap, chrom_lengths.get(chrom,None))
                    for mpeak in mpeaks:

                        print(mpeak.bed_str(chrom, strand), file=final_bed_out)

                    open_end = peak.end
                    open_peaks = [peak]

                else:
                    # extend open peak
                    open_peaks.append(peak)

                    open_end = max(open_end, peak.end)


        if len(open_peaks) > 0:

            mpeaks = merge_peaks(open_peaks, options.feature_size, options.merge_overlap, chrom_lengths.get(chrom,None))


            for mpeak in mpeaks:

                print(mpeak.bed_str(chrom, strand), file=final_bed_out)

    final_bed_out.close()

    #clean
    for chrom_key in chrom_files:
        os.remove(chrom_files[chrom_key])

    final_act_out = open('%s_act.txt' % options.out_prefix, 'w')

    cols = [''] + db_targets
    print('\t'.join(cols), file=final_act_out)


    for line in open('%s.bed' % options.out_prefix):
        a = line.rstrip().split('\t')

        peak_id = '%s:%s-%s(%s)' % (a[0], a[1], a[2], a[5])


        peak_act = [0]*len(db_targets)
        for ai in a[6].split(','):
            if ai != '.':
                peak_act[int(ai)] = 1

        cols = [peak_id] + peak_act
        print('\t'.join([str(c) for c in cols]), file=final_act_out)

    final_act_out.close()

    final_bed_out_loc = "{}.bed".format(options.out_prefix)
    final_act_out_loc = "{}_act.txt".format(options.out_prefix)
    df_bed = pd.read_table(final_bed_out_loc, header=None, index_col=None)
    df_bed.columns = ['chr', 'start', 'end', 'peak_name', 'score', 'strand', 'act']
    dup_bed = df_bed[df_bed.duplicated(subset=['chr', 'start', 'end'], keep=False)]
    if len(dup_bed) != 0:

        grouped_df_bed = df_bed.groupby(list(df_bed.columns[df_bed.columns!="act"]), as_index=False,sort=False).agg({'act': lambda x: ','.join(x)})
        grouped_df_bed.to_csv("{}_uniq.bed".format(options.out_prefix),sep="\t",index=False,header=False)
        dup_bed.to_csv("{}_dup.bed".format(options.out_prefix), sep="\t",index=False, header=False)

        df_act = pd.read_table(final_act_out_loc, header=0, index_col=None)
        df_act = df_act.rename(columns={df_act.columns[0]: 'peak'})
        grouped_df_act = df_act.groupby('peak', as_index=False, sort=False).sum()
        grouped_df_act = grouped_df_act.map(lambda x: 1 if (isinstance(x, (int, float)) and x > 1) else x)
        grouped_df_act.rename(columns={grouped_df_act.columns[0]:''},inplace=True)

        dup_act = df_act[df_act.duplicated(subset=['peak'], keep=False)]
        grouped_df_act.to_csv("{}_uniq_act.txt".format(options.out_prefix), sep="\t",
                              index=False, header=True)
        dup_act.to_csv("{}_dup_act.txt".format(options.out_prefix), sep="\t",
                       index=False, header=True)

def activity_set(act_cs):
    ''' Return a set of ints from a comma-separated list of int strings.

    Attributes:
        act_cs (str) : comma-separated list of int strings

    Returns:
        set (int) : int's in the original string
    '''
    ai_strs = [ai for ai in act_cs.split(',')]

    if ai_strs[-1] == '':
        ai_strs = ai_strs[:-1]

    if ai_strs[0] == '.':
        aset = set()
    else:
        aset = set([int(ai) for ai in ai_strs])

    return aset

def find_midpoint(start, end):

    mid = (start + end)//2
    return mid

def merge_peaks(peaks, peak_size, merge_overlap, chrom_len):

    max_overlap = merge_overlap
    while len(peaks) > 1 and max_overlap >= merge_overlap:

        max_i = 0
        max_overlap = peaks[0].end - peaks[1].start
        for i in range(1, len(peaks)-1):
            peaks_overlap = peaks[i].end - peaks[i+1].start
            if peaks_overlap > max_overlap:
                max_i = i
                max_overlap = peaks_overlap

        if max_overlap >= merge_overlap:
            peaks[max_i].merge(peaks[max_i+1], peak_size, chrom_len)
            peaks = peaks[:max_i+1] + peaks[max_i+2:]

    return peaks



class Peak:

    def __init__(self, start, end, act):
        self.start = start
        self.end = end
        self.act = act


    def bed_str(self, chrom, strand):

        if len(self.act) == 0:
            act_str = '.'
        else:
            act_str = ','.join([str(ai) for ai in sorted(list(self.act))])
        cols = (chrom, str(self.start), str(self.end), '.', '1', strand, act_str) #全变成1了
        return '\t'.join(cols)

    def merge(self, peak2, ext_len, chrom_len):


        peak_mids = [find_midpoint(self.start,self.end)]
        peak_mids.append(find_midpoint(peak2.start,peak2.end))

        # weight peaks
        peak_weights = [1+len(self.act)]
        peak_weights.append(1+len(peak2.act))


        merge_mid = int(0.5+np.average(peak_mids, weights=peak_weights))


        merge_start = max(0, merge_mid - ext_len//2)
        merge_end = merge_start + ext_len
        if chrom_len and merge_end > chrom_len:
            merge_end = chrom_len
            merge_start = merge_end - ext_len



        merge_act = self.act | peak2.act


        self.start = merge_start
        self.end = merge_end
        self.act = merge_act
def split_region(lst,ext_len,direc,chrom_len):
    cut_off=int(np.ceil(ext_len*0.5))
    region_range = lst[:]
    if direc=='up':
        if len(region_range)-1 <= ext_len:
            region_range=[list(range(region_range[-1]-ext_len, region_range[-1]+1))]
        else:
            region_range = [lst[i-(ext_len+1):i] if i - (ext_len+1) >= 0 else lst[0:i] for i in range(len(lst), 0, -(ext_len+1))][::-1]
            if len(region_range[0])<=cut_off:
                region_range=region_range[1:]
            else:
                region_range=supply_region(lst=region_range[:], ext_len=ext_len,direc=direc,chrom_len=chrom_len)
    elif direc=='dn':
        if len(region_range)-1<=ext_len:
            region_range = [list(range(region_range[0], region_range[0] + ext_len+1))]
        else:
            region_range = [lst[i:i+(ext_len+1)] for i in range(0, len(lst), (ext_len+1))]
            if len(region_range[-1])<=cut_off:
                region_range=region_range[:-1]
            else:
                region_range=supply_region(lst=region_range[:], ext_len=ext_len,direc=direc,chrom_len=chrom_len)
    return region_range

def supply_region(lst,ext_len,direc,chrom_len):
    if direc=='up':
        lst[0] = [x for x in range(lst[0][-1] + 1 - (ext_len+1), lst[0][-1] + 1)]
        if lst[0][0] < 0:
            lst[0] = list(range(0, (ext_len+1)))
    elif direc=='dn':
        lst[-1] = [x for x in range(lst[-1][0], lst[-1][0]+(ext_len+1))]
        if chrom_len and lst[-1][-1] > chrom_len:
            lst[-1]=list(range(chrom_len+1 - (ext_len+1), chrom_len + 1))
    return lst

def seperate(start, end, flank_len, ext_len, chrom_len, flank_ratio,direc_for_ext,direc_for_split):

    if flank_len != 0:
        peak_width = end - start
        cut_off = flank_ratio * peak_width  # 2400
        flank_num = math.ceil(cut_off / flank_len)#1
        if direc_for_ext=="both":
            start = max(0, start - flank_len * flank_num)
            end = end + flank_len * flank_num
        elif direc_for_ext=="up":
            start = max(0, start - flank_len * flank_num)
            end = end
        elif direc_for_ext=="dn":
            start = start
            end = end + flank_len * flank_num
        else:
            raise ValueError("Please input both/up/dn in -d ")
    else:
        start = start
        end = end
    if chrom_len and end > chrom_len:
        end = chrom_len

    peak_full_range = list(range(start, end+1))

    if direc_for_split=='up':
        peak_sep_range = split_region(lst=peak_full_range[:], ext_len=ext_len, direc='up', chrom_len=chrom_len)
    elif direc_for_split=='dn':
        peak_sep_range = split_region(lst=peak_full_range[:], ext_len=ext_len, direc='dn', chrom_len=chrom_len)
    elif direc_for_split=='both':
        mid = find_midpoint(start, end)
        mid_start = mid-ext_len//2
        mid_end = mid_start+ext_len
        if len(peak_full_range)-1 <= ext_len+0.5*ext_len:
            peak_sep_range = [list(range(mid_start, mid_end+1))]
        else:
            mid_up_region = peak_full_range[:peak_full_range.index(mid_start)]
            mid_dn_region = peak_full_range[peak_full_range.index(mid_end)+1:]
            mid_up_split = split_region(lst=mid_up_region[:], ext_len=ext_len, direc='up', chrom_len=chrom_len)
            mid_dn_split = split_region(lst=mid_dn_region[:], ext_len=ext_len, direc='dn', chrom_len=chrom_len)
            peak_sep_range=[]
            for i in mid_up_split:
                peak_sep_range.append(i)
            peak_sep_range.append(list(range(mid_start,mid_end+1)))
            for j in mid_dn_split:
                peak_sep_range.append(j)

    obj_dict = {}
    counter = 0
    for item in peak_sep_range:
        key = f"subpeak_{counter}"
        value = [item[0], item[-1]]
        obj_dict[key]= value
        counter+=1
    return obj_dict




if __name__ == '__main__':

    main()
