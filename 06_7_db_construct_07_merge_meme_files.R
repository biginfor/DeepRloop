#just TFs
header<-"/home/fangzj/Workdir/Basset_Unblocked/data/meme_db/header_for_tf_chip.txt"
JASPAR_TFs_human<-"/home/fangzj/Workdir/Basset_Unblocked/data/meme_db/JASPAR_23_11_28/JASPAR_human_fil.meme"
JASPAR_TFs_mouse<-"/home/fangzj/Workdir/Basset_Unblocked/data/meme_db/JASPAR_23_11_28/JASPAR_mouse_fil.meme"
CISBP_TFs_remained<-"/home/fangzj/Workdir/Basset_Unblocked/data/meme_db/Homo_sapiens_CIS_BP_2023_10_25/remained_CISBP_TFs_withoutheader.meme"
RLBP_denovo_human<-"/home/fangzj/Workdir/Basset_Unblocked/data/encode/Chip_Seq/meme_files/allRLBPs/denovo_sel_max_RLBPs_noheader.meme"
#
TFs_PFMs_output<-"/home/fangzj/Workdir/Basset_Unblocked/data/meme_db/just_TFs/"
#
merge_script<-paste0("cat ",header," ", JASPAR_TFs_human," ",JASPAR_TFs_mouse," ",CISBP_TFs_remained," > ",TFs_PFMs_output,"TFs_withheader.meme")
system(merge_script)
#
TFs_RLBPs_PFMs_output<-"/home/fangzj/Workdir/Basset_Unblocked/data/meme_db/TFs_RLBPs/"
merge_script_all<-paste0("cat ",header," ", JASPAR_TFs_human," ",JASPAR_TFs_mouse," ",CISBP_TFs_remained," ",RLBP_denovo_human," > ",TFs_RLBPs_PFMs_output,"TFs_RLBPs_withheader.meme")
system(merge_script_all)
