#ref https://github.com/davek44/Basset
#!/usr/bin/env python
from optparse import OptionParser
import os, subprocess
import h5py
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import myutils.dna_io as dna_io

#####################################################################
#explore the 1st Conv layer
#####################################################################
#weblogo param
weblogo_opts = '-X YES -Y YES --errorbars YES --fineprint ""'
weblogo_opts += ' -C "#CB2026" A A'
weblogo_opts += ' -C "#34459C" C C'
weblogo_opts += ' -C "#FBB116" G G'
weblogo_opts += ' -C "#0C8040" T T'

def main():
    usage = 'usage: %prog [options] <input_reprs_file> <input_weights_file>' #<test_hdf5_file>'
    parser = OptionParser(usage)
    parser.add_option('-a', dest='act_t', default=0.5, type='float', help='Activation threshold (as proportion of max) to consider for PWM [Default: %default]')
    parser.add_option('-o', dest='out_dir', default='.')
    parser.add_option('-m', dest='meme_db', default='./data/TFs_RLBPs_withheader_IDs.meme')
    parser.add_option('-t', dest='trim_filters', default=True, action='store_true', help='Trim uninformative positions off the filter ends [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('must provide input_file')
    else:
        input_reprs_file = args[0]
        input_weights_file = args[1]


    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)
    out_dir = options.out_dir 
    #################################################################
    # load data
    #################################################################
    with h5py.File(input_weights_file,'r') as h5_2:
        seqs_bool = np.array(h5_2['test_in'])
        filter_weights = np.array(h5_2['weights'])

    num_filters = filter_weights.shape[0]  
    filter_size = filter_weights.shape[2]

    with h5py.File(input_reprs_file,'r') as h5_1:
        filter_outs = np.array(h5_1['outs'])
    # bool 2 letters
    seqs = dna_io.vecs2dna(seqs_bool)
    
    #################################################################
    # individual filter plots
    #################################################################
    filters_ic = []
    meme_out = meme_intro(meme_file="{}/filters_meme.txt".format(out_dir),seqs=seqs)
    meme_out = "{}/filters_meme.txt".format(out_dir)

    for f in range(0,num_filters):
        print('Filter {}'.format(f))

        # write possum motif file
        filter_possum(filter_weights[f,:,:], f"filter{f}", f"{out_dir}/filter{f}_possum.txt", False)

        # plot weblogo of high scoring outputs
        plot_filter_logo(filter_outs[:,f,:], filter_size, seqs, f"{out_dir}/filter{f}_logo", maxpct_t=options.act_t)

        # make a PWM for the filter
        filter_pwm, nsites = make_filter_pwm(f"{out_dir}/filter{f}_logo.fa")
        print(nsites)
        if nsites < 10:
            #no information
            filters_ic.append(0)
        else:
            #compute and save information content
            filters_ic.append(info_content(filter_pwm))

            #add to the meme motif file
            meme_add(meme_out, f, filter_pwm, nsites, True)
            
    #################################################################
    # annotate filters
    #################################################################
    # tom tom 
    tomtom_cmd = 'conda run -n MEME tomtom -dist pearson -thresh 0.5 -oc {}/tomtom {}/filters_meme.txt {}'.format(out_dir, out_dir,options.meme_db)
    subprocess.call(tomtom_cmd, shell=True)
    subprocess.call("sed '/^\s*#/d; /^$/d' {}/tomtom/tomtom.tsv > {}/tomtom/tomtom_del.tsv".format(out_dir,out_dir), shell=True)
    filter_names = name_filters(num_filters, '{}/tomtom/tomtom_del.tsv'.format(out_dir), options.meme_db)
    
    
    #################################################################
    # print a table of information
    #################################################################

    table_out = open('{}/table.txt'.format(out_dir), 'w')
    # print header for later panda reading
    header_cols = ('', 'consensus', 'annotation', 'ic', 'mean', 'std')
    print('{:<3s}  {:<19s}  {:<10s}  {:<5s}  {:<6s}  {:<6s}'.format(*header_cols), file=table_out)

    for f in range(num_filters):
        # collapse to a consensus motif
        consensus = filter_motif(filter_weights[f,:,:])
        #print(consensus)
        # grab annotation
        annotation = '.'
        name_pieces = filter_names[f].split('_')
        if len(name_pieces) > 1:
            annotation = name_pieces[1]

        # plot density of filter output scores
        fmean, fstd = plot_score_density(np.ravel(filter_outs[:,f,:]), '{}/filter{}_dens.pdf'.format(options.out_dir,f))
        row_cols = (f, consensus, annotation, filters_ic[f],fmean, fstd)
        print(row_cols)

        print('{:<-3d}  {:<19s}  {:<10s}  {:<5.2f}  {:<6.4f}  {:<6.4f}'.format(*row_cols), file=table_out)

    table_out.close()
    
#################################################################
# utility functions
#################################################################
def get_motif_proteins(meme_db_file):
    ''' Hash motif_id's to protein names using the MEME DB file '''
    motif_protein = {}
    for line in open(meme_db_file):
        a = line.split()
        if len(a) > 0 and a[0] == 'MOTIF':
            if a[2][0] == '(':
                motif_protein[a[1]] = a[2][1:a[2].find(')')]
            else:
                motif_protein[a[1]] = a[2]
    return motif_protein


def info_content(pwm, transpose=False, bg_gc=0.415):
    ''' Compute PWM information content.
    In the original analysis, I used a bg_gc=0.5. For any
    future analysis, I ought to switch to the true hg38
    value of 0.415.
    '''
    pseudoc = 1e-9

    if transpose:
        pwm = np.transpose(pwm)

    bg_pwm = [1-bg_gc, bg_gc, bg_gc, 1-bg_gc]

    ic = 0
    for i in range(pwm.shape[0]):
        for j in range(4):
            # ic += 0.5 + pwm[i][j]*np.log2(pseudoc+pwm[i][j])
            ic += -bg_pwm[j]*np.log2(bg_pwm[j]) + pwm[i][j]*np.log2(pseudoc+pwm[i][j])

    return ic

def make_filter_pwm(filter_fasta):
    ''' Make a PWM for this filter from its top hits '''

    nts = {'A':0, 'C':1, 'G':2, 'T':3}
    pwm_counts = []
    nsites = 4 # pseudocounts
    for line in open(filter_fasta):
        if line[0] != '>':
            seq = line.rstrip()
            nsites += 1
            if len(pwm_counts) == 0:
                # initialize with the length
                for i in range(len(seq)):
                    pwm_counts.append(np.array([1.0]*4))

            # count
            for i in range(len(seq)):
                try:
                    pwm_counts[i][nts[seq[i]]] += 1
                except KeyError:
                    pwm_counts[i] += np.array([0.25]*4)

    # normalize
    pwm_freqs = []
    for i in range(len(pwm_counts)):
        pwm_freqs.append([pwm_counts[i][j]/float(nsites) for j in range(4)])

    return np.array(pwm_freqs), nsites-4



def meme_add(meme_out, f, filter_pwm, nsites, trim_filters=True):
    ''' Print a filter to the growing MEME file
    Attrs:
        meme_out : open file
        f (int) : filter index #
        filter_pwm (array) : filter PWM array
        nsites (int) : number of filter sites
    '''
    if not trim_filters:
        ic_start = 0 #
        ic_end = filter_pwm.shape[0]-1 #18
    else:
        ic_t = 0.2

        # trim PWM of uninformative prefix
        ic_start = 0
        while ic_start < filter_pwm.shape[0] and info_content(filter_pwm[ic_start:ic_start+1]) < ic_t:
            ic_start += 1

        # trim PWM of uninformative suffix
        ic_end = filter_pwm.shape[0]-1
        while ic_end >= 0 and info_content(filter_pwm[ic_end:ic_end+1]) < ic_t:
            ic_end -= 1        
    if ic_start < ic_end:
        with open(meme_out, 'a') as meme_file:
            print('MOTIF filter%d' % f, file=meme_file)
            print('letter-probability matrix: alength= 4 w= %d nsites= %d E= 0' % (ic_end-ic_start+1, nsites), file=meme_file)
            for i in range(ic_start, ic_end+1):#0 , 19
                print('%.4f %.4f %.4f %.4f' % tuple(filter_pwm[i]), file=meme_file) #trouble
            print('', file=meme_file)



def meme_intro(meme_file, seqs):
    ''' Open MEME motif format file and print intro
    Attrs:
        meme_file (str) : filename
        seqs [str] : list of strings for obtaining background freqs
    Returns:
        mem_out : open MEME file
    '''
    nts = {'A':0, 'C':1, 'G':2, 'T':3}

    # count
    nt_counts = [1]*4
    for i in range(len(seqs)):
        for nt in seqs[i]:
            try:
                nt_counts[nts[nt]] += 1
            except KeyError:
                pass

    # normalize
    nt_sum = float(sum(nt_counts))
    nt_freqs = [nt_counts[i]/nt_sum for i in range(4)]

    # open file for writing
    meme_out = open(meme_file, 'w')

    # print intro material
    print('MEME version 5.0.5', file=meme_out)
    print('', file=meme_out)
    print('ALPHABET= ACGT', file=meme_out)
    print('', file=meme_out)
    print('Background letter frequencies:', file=meme_out)
    print(f'A {nt_freqs[0]:.4f} C {nt_freqs[1]:.4f} G {nt_freqs[2]:.4f} T {nt_freqs[3]:.4f}', file=meme_out)
    print('', file=meme_out)
    #meme_out.close()
    return meme_out



def name_filters(num_filters, tomtom_file, meme_db_file):
    ''' Name the filters using Tomtom matches.
    Attrs:
        num_filters (int) : total number of filters
        tomtom_file (str) : filename of Tomtom output table.
        meme_db_file (str) : filename of MEME db
    Returns:
        filter_names [str] :
    '''
    # name by number
    filter_names = ['f%d'%fi for fi in range(num_filters)]

    # name by protein
    if tomtom_file is not None and meme_db_file is not None:
        motif_protein = get_motif_proteins(meme_db_file)

        # hash motifs and q-value's by filter
        filter_motifs = {}

        tt_in = open(tomtom_file)
        tt_in.readline()
        for line in tt_in:
            a = line.split()
            fi = int(a[0][6:])
            motif_id = a[1]
            qval = float(a[5])

            filter_motifs.setdefault(fi,[]).append((qval,motif_id))

        tt_in.close()

        # assign filter's best match
        for fi in filter_motifs:
            top_motif = sorted(filter_motifs[fi])[0][1]
            filter_names[fi] += '_%s' % motif_protein[top_motif]

    return np.array(filter_names)



def filter_motif(param_matrix):
    nts = 'ACGT'

    motif_list = []
    for v in range(param_matrix.shape[1]):
        max_n = 0
        for n in range(1,4):
            if param_matrix[n,v] > param_matrix[max_n,v]:
                max_n = n

        if param_matrix[max_n,v] > 0:
            motif_list.append(nts[max_n])
        else:
            motif_list.append('N')

    return ''.join(motif_list)


def filter_possum(param_matrix, motif_id, possum_file, trim_filters=False, mult=200):
    # possible trim
    trim_start = 0
    trim_end = param_matrix.shape[1]-1
    trim_t = 0.3
    if trim_filters:
        # trim PWM of uninformative prefix
        while trim_start < param_matrix.shape[1] and np.max(param_matrix[:,trim_start]) - np.min(param_matrix[:,trim_start]) < trim_t:
            trim_start += 1

        # trim PWM of uninformative suffix
        while trim_end >= 0 and np.max(param_matrix[:,trim_end]) - np.min(param_matrix[:,trim_end]) < trim_t:
            trim_end -= 1
 
    if trim_start < trim_end:
        with open(possum_file, 'w') as possum_out:
            print('BEGIN GROUP', file=possum_out)
            print('BEGIN FLOAT', file=possum_out)
            print('ID %s' % motif_id, file=possum_out)
            print('AP DNA', file=possum_out)
            print('LE %d' % (trim_end+1-trim_start), file=possum_out)
            for ci in range(trim_start,trim_end+1):
                print('MA %s' % ' '.join(['%.2f'%(mult*n) for n in param_matrix[:,ci]]), file=possum_out)
            print('END', file=possum_out)
            print('END', file=possum_out)

def plot_filter_logo(filter_outs, filter_size, seqs, out_prefix, raw_t=0, maxpct_t=None):
    if maxpct_t:
        all_outs = np.ravel(filter_outs) # flatten
        all_outs_mean = all_outs.mean()
        all_outs_norm = all_outs - all_outs_mean
        raw_t = maxpct_t * all_outs_norm.max() + all_outs_mean

    # SAME padding
    pad_side = (filter_size - 1) // 2

    # print fasta file of positive outputs
    filter_fasta_out = open('%s.fa' % out_prefix, 'w')
    filter_count = 0
    for i in range(filter_outs.shape[0]):
        for j in range(pad_side, filter_outs.shape[1]-pad_side):
            if filter_outs[i,j] > raw_t:
                js = j - pad_side
                kmer = seqs[i][js:js+filter_size]
                print('>%d_%d' % (i, j), file=filter_fasta_out)
                print(kmer, file=filter_fasta_out)
                filter_count += 1
    filter_fasta_out.close()
    print(filter_count)
    # make weblogo
    if filter_count > 0:
        weblogo_cmd = 'conda run -n MEME `weblogo {} < {}.fa > {}.eps`'.format(weblogo_opts, out_prefix, out_prefix)
        subprocess.call(weblogo_cmd, shell=True)
        print(weblogo_cmd)



def plot_score_density(f_scores, out_pdf):
    sns.set(font_scale=1.3)
    plt.figure()
    sns.displot(f_scores, kde=False)
    plt.xlabel('ReLU output')
    plt.savefig(out_pdf)
    plt.close()

    return f_scores.mean(), f_scores.std()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()