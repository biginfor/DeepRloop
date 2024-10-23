library("dplyr")
#
output_dir<-"./data/curation/"
jaspar_fil<-read.table(paste0(output_dir,"internal_motifs_curation_filtered.tab"),sep = "\t",header = T)
#remove redundancy
TF_max_idx<-jaspar_fil %>%group_by(Dataset_TF_NAME) %>%reframe(max_score = max(Cent_pval), max_score_TF= curation_id[Cent_pval == max(Cent_pval)])
sel_motifs_col<-file(paste0(output_dir,'selected_motifs.jaspar'),'r')
max_motifs_col<-file(paste0(output_dir,'selected_max_motifs.jaspar'),'w')
sel_motifs<-readLines(sel_motifs_col)
close(sel_motifs_col)
for(i in TF_max_idx$max_score_TF){
  idx_start<-grep(i,sel_motifs)
  idx_end<-idx_start+4
  writeLines(sel_motifs[idx_start:idx_end],max_motifs_col)
}
close(max_motifs_col)

#change file format in shell!
#inputjaspar=./data/curation/selected_max_motifs.jaspar
#outputmeme=./data/allRLBPs/denovo_sel_max_RLBPs.meme
#bg_file=./data/bkgd.txt
#conda activate MEME
#jaspar2meme -bundle -bg ${bg_file} ${inputjaspar}> ${outputmeme}
#tail -n +10 ${outputmeme} > ./data/allRLBPs/denovo_sel_max_RLBPs_noheader.meme