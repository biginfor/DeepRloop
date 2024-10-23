library("dplyr")
library("tidyverse")
library("ggplot2")
library("ggpubr")
library("data.table")
source("/home/fangzj/Workdir/Basset_Unblocked/scripts/project/dn_analysis/get_seq_loc_utils.R")
library("rlist")
#
ext_dire="both_ext";split_dire="up_split";ext_bps_nu="85"
Phylop_loc<-"/media/fangzj/data/ext_data/phyloP100way/hg38.phyloP100way.bedGraph"
file_dir<-paste0('/home/fangzj/Workdir/Basset_Unblocked/data/cleandata/flank_seqs/all_classes_flank/',ext_dire,'/',split_dire,'/1e-90/','200_600_',ext_bps_nu)
#get_Phylop_mean<-function(ext_dire,split_dire,ext_bps_nu)
query_file<-paste0(file_dir,'/bin_200_600_85_act.txt')
act_file<-read.table(query_file,sep="\t",header = T,row.names = 1)
peak_count<-as.data.frame(apply(act_file,1,sum))
peak_count <- cbind(rownames(peak_count),peak_count)
colnames(peak_count)<-c("region","count")
peak_count<-dplyr::arrange(peak_count,desc(count))
fivenums<-fivenum(peak_count$count)
peak_count[peak_count$count<fivenums[3],3]<-"low_frequency"# <2
peak_count[peak_count$count>=fivenums[3]&peak_count$count<fivenums[4],3]<-"mid_frequency"# >=2 <5
peak_count[peak_count$count>=fivenums[4],3]<-"high_frequency"#>=5
colnames(peak_count)[3]<-"class"
peak_count<-separate(peak_count,col = "region",into = c("chr","start","end","_"),sep = ":|-|\\(\\+\\)")

#
for(i in unique(peak_count$class)){
  print(i)
  tmp<-peak_count[peak_count$class==i,]
  region<-tmp[,1:3]
  region$start<-as.numeric(region$start)
  region$end<-as.numeric(region$end)
  #chr_levels<-c(paste0("chr",seq(1,22)),"chrX","chrY")
  #region$chr<-factor(region$chr,levels = chr_levels)
  region<-arrange(region,chr,start,end)
  bed_loc<-paste0(file_dir,"/",i,"_peak_region.bed")
  write.table(region,bed_loc,sep = "\t",quote = F,row.names = F,col.names = F)
  Phylop_mean_out<-paste0(file_dir,"/",i,"_peak_phylop.bed")
  bedtools<-"/home/fangzj/Workdir/software/anaconda/envs/basset/bin/bedtools"
  code<-paste0(bedtools," map -a ",bed_loc," -b ",Phylop_loc," -c 4 -o mean > ",Phylop_mean_out)
  system(code)
}
#
low_group_Phylop<-read.table(paste0(file_dir,"/low_frequency_peak_phylop.bed"))
low_group_Phylop$V4[low_group_Phylop$V4=="."]<-0
low_group_Phylop$V4<-as.numeric(low_group_Phylop$V4)
low_group_Phylop$"V5"<-"low_frequency_peak_phylop"
#
mid_group_Phylop<-read.table(paste0(file_dir,"/mid_frequency_peak_phylop.bed"))
mid_group_Phylop$V4[mid_group_Phylop$V4=="."]<-0
mid_group_Phylop$V4<-as.numeric(mid_group_Phylop$V4)
mid_group_Phylop$"V5"<-"mid_frequency_peak_phylop"
#
high_group_Phylop<-read.table(paste0(file_dir,"/high_frequency_peak_phylop.bed"))
high_group_Phylop$V4[high_group_Phylop$V4=="."]<-0
high_group_Phylop$V4<-as.numeric(high_group_Phylop$V4)
high_group_Phylop$"V5"<-"high_frequency_peak_phylop"

#
Phylop_df<-bind_rows(low_group_Phylop,mid_group_Phylop,high_group_Phylop)
colnames(Phylop_df)<-c("chr","start","end","mean_Phylop","class")
Phylop_df$class<-factor(Phylop_df$class,levels = c("low_frequency_peak_phylop","mid_frequency_peak_phylop","high_frequency_peak_phylop"))
write.table(Phylop_df,paste0(file_dir,"/Phylop_df.txt"),sep="\t",quote = F,col.names = T,row.names = F)


#
my_comparisons<-list(c("low_frequency_peak_phylop","mid_frequency_peak_phylop"),c("mid_frequency_peak_phylop","high_frequency_peak_phylop"),c("low_frequency_peak_phylop","high_frequency_peak_phylop"))
mean_Phylop<-aggregate(mean_Phylop~class,Phylop_df,mean)
mid_Phylop<-aggregate(mean_Phylop~class,Phylop_df,median)
anno_text<-paste0("mean:\n","low_group:",mean_Phylop[1,2],"\nmid_group:",mean_Phylop[2,2],"\nhigh_group:",mean_Phylop[3,2],
                  "\nmid:\n","low_group:",mid_Phylop[1,2],"\nmid_group:",mid_Phylop[2,2],"\nhigh_group:",mid_Phylop[3,2])

p1<-ggplot(Phylop_df,aes(x=class,y=mean_Phylop))+geom_violin()+geom_boxplot(width=0.2,position = position_dodge(0.9),outlier.size = 0.1)+stat_compare_means(comparisons = my_comparisons,method = "t.test",size=5)+theme_classic()+theme(text = element_text(size = 15))
p2<-p1+annotate("text",x=0.7,y=7,label=anno_text)
ggsave(paste0(file_dir,"/act_file_box.pdf"),p2,width = 15,height = 7.5)

#RANDOM SEQs
random_fa<-"/home/fangzj/Workdir/Basset_Unblocked/data/random_Seqs/fa/chr1_seed_1_start_66608964_no_20000.fa"
random_fa_txt<-readLines(random_fa)
random_fa_bed<-as.data.frame(random_fa_txt[grep("chr",random_fa_txt)])
colnames(random_fa_bed)<-"region"
random_fa_bed<-separate(random_fa_bed,region,into = c("_","chr","start","end"),sep = ">|:|\\-")
random_fa_bed<-random_fa_bed[,2:4]
random_fa_bed$start<-as.numeric(random_fa_bed$start)
random_fa_bed$end<-as.numeric(random_fa_bed$end)
write.table(random_fa_bed,"/home/fangzj/Workdir/Basset_Unblocked/data/random_Seqs/fa/chr1_seed_1_start_66608964_no_20000.bed")
random_fa_bed_sort<-arrange(random_fa_bed,chr,start,end)
bed_loc<-"/home/fangzj/Workdir/Basset_Unblocked/data/random_Seqs/fa/Phylop_loc/chr1_seed_1_start_66608964_no_20000_sort.bed"
write.table(random_fa_bed_sort,bed_loc,sep="\t",col.names = F,row.names = F,quote = F)
Phylop_mean_out<-paste0("/home/fangzj/Workdir/Basset_Unblocked/data/random_Seqs/fa/Phylop_loc/chr1_seed_1_start_66608964_no_20000_Phylop.bed")
bedtools<-"/home/fangzj/Workdir/software/anaconda/envs/basset/bin/bedtools"
code<-paste0(bedtools," map -a ",bed_loc," -b ",Phylop_loc," -c 4 -o mean > ",Phylop_mean_out)
system(code)
Phylop_random_seqs <- read.table("/home/fangzj/Workdir/Basset_Unblocked/data/random_Seqs/fa/Phylop_loc/chr1_seed_1_start_66608964_no_20000_Phylop.bed")
Phylop_random_seqs$V4[Phylop_random_seqs$V4=="."]<-0
Phylop_random_seqs$V4<-as.numeric(Phylop_random_seqs$V4)
Phylop_random_seqs$"class"<-"random_seqs_phylop"
colnames(Phylop_random_seqs)<-c("chr","start","end","mean_Phylop","class")
mean(Phylop_random_seqs$mean_Phylop)
#MOTIF
motif_loc<-paste0(file_dir,"/quick_check/m1d_optrmsprop_drop0.39_lr0.00002_step14_gama0.28_mom0.97_batch512_ks1_12_conv1_300_norm_1e9,0_run_3/samples_1000/TF_RLBP_comb/Rloop_motifs/")
motif_file_loc<-paste0(motif_loc,"Rloop_motifs.variant_function")
motif_file<-fread(input = motif_file_loc,nThread = 25,sep="\t")
motif_file<-as.data.frame(motif_file[,3:5])
motif_file<-arrange(motif_file,V3,V4,V5)
motif_bed_loc<-paste0(motif_loc,"Rloop_motifs.bed")
write.table(motif_file,motif_bed_loc,sep = "\t",quote = F,row.names = F,col.names = F)
bedtools<-"/home/fangzj/Workdir/software/anaconda/envs/basset/bin/bedtools"
Phylop_mean_out_motif<-paste0(motif_loc,"motif_phylop.bed")
code<-paste0(bedtools," map -a ",motif_bed_loc," -b ",Phylop_loc," -c 4 -o mean > ",Phylop_mean_out_motif)
system(code)
motif_Phylop<-read.table(paste0(motif_loc,"/motif_phylop.bed"))
motif_Phylop$V4[motif_Phylop$V4=='.']<-0
motif_Phylop$V4<-as.numeric(motif_Phylop$V4)
motif_Phylop$"calss"<-"motif_phylop"
colnames(motif_Phylop)<-c("chr","start","end","mean_Phylop","class")
mean(motif_Phylop$mean_Phylop)
#high exp genes
highexp_bed_loc<-"/home/fangzj/Workdir/Basset_Unblocked/data/GTEx/GTEx_V8/compare_high_exp_genes/"
get_Phylop<-function(ratio){
  highexp_bed <- paste0(highexp_bed_loc,"high_exp_600ps_exp_ratio_",ratio,".bed")
  high_exp_ratio<-read.table(highexp_bed)
  high_exp_ratio<-plyr::arrange(high_exp_ratio,V1,V2,V3)
  high_exp_ratio_sorted<-paste0(highexp_bed_loc,"/Phylop_score/high_exp_600ps_exp_ratio_",ratio,"_sorted.bed")
  write.table(x = high_exp_ratio,high_exp_ratio_sorted,row.names = F,col.names = F,sep = "\t",quote = F)
  bedtools<-"/home/fangzj/Workdir/software/anaconda/envs/basset/bin/bedtools"
  Phylop_loc<-"/media/fangzj/data/ext_data/phyloP100way/hg38.phyloP100way.bedGraph"
  Phylop_mean_out_highexp<-paste0(highexp_bed_loc,"Phylop_score/highexp_phylop_",ratio,".bed")
  code<-paste0(bedtools," map -a ",high_exp_ratio_sorted," -b ",Phylop_loc," -c 4 -o mean > ",Phylop_mean_out_highexp)
  if(!file.exists(Phylop_mean_out_highexp)){
    system(code)
  }
  highexp_Phylop<-read.table(Phylop_mean_out_highexp)
  highexp_Phylop$V4[highexp_Phylop$V4=='.']<-0
  highexp_Phylop$V4<-as.numeric(highexp_Phylop$V4)
  highexp_Phylop$"class"<-paste0("exp_ratio_",ratio)
  colnames(highexp_Phylop)<-c("chr","start","end","mean_Phylop","class")
  return(highexp_Phylop)
}
raio_001<-get_Phylop(0.01)
raio_005<-get_Phylop(0.05)
ratio_02<-get_Phylop(0.2)

#PLOT

Phylop_df_all<-rbind(Phylop_df,Phylop_random_seqs,raio_001,raio_005,ratio_02)
mean_Phylop<-aggregate(mean_Phylop~class,Phylop_df_all,mean)
mid_Phylop<-aggregate(mean_Phylop~class,Phylop_df_all,median)
anno_text<-paste0("mean:\n","low_group:",mean_Phylop[1,2],"\nmid_group:",mean_Phylop[2,2],"\nhigh_group:",mean_Phylop[3,2],"\nrandom_seqs:",mean_Phylop[4,2],"\nmotifs:",mean_Phylop[5,2],
                  "\nmid:\n","low_group:",mid_Phylop[1,2],"\nmid_group:",mid_Phylop[2,2],"\nhigh_group:",mid_Phylop[3,2],
                  "\nrandom_seqs:",mid_Phylop[4,2],"\nmotifs:",mid_Phylop[5,2]
                  )
Phylop_df_all$class<-factor(Phylop_df_all$class,levels =c("random_seqs_phylop","low_frequency_peak_phylop","mid_frequency_peak_phylop","high_frequency_peak_phylop","exp_ratio_0.2","exp_ratio_0.05","exp_ratio_0.01"))
my_comparisons<-list(c("random_seqs_phylop","low_frequency_peak_phylop"),c("random_seqs_phylop","mid_frequency_peak_phylop"),c("random_seqs_phylop","high_frequency_peak_phylop"),c("random_seqs_phylop","motif_phylop"))
#箱线图
p3<-ggplot(Phylop_df_all,aes(x=class,y=mean_Phylop))+geom_violin()+geom_boxplot(width=0.2,position = position_dodge(0.9),outlier.size = 0)+stat_compare_means(comparisons = my_comparisons,method = "t.test",size=5)+theme_classic()+theme(text = element_text(size = 15))#+annotate("text",x=5,y=7,label=anno_text)
ggsave(paste0(file_dir,"/Phylop_box.png"),p3,width = 15,height = 10)
#密度图
colors<-c("black","#FFBE7A","purple","red","#82B0D2","#2E8B57","#DAA520")
names(colors)<-levels(Phylop_df_all$class)
p4<-ggplot(Phylop_df_all,aes(x=mean_Phylop))+geom_density(aes(color=class))+scale_color_manual(values = colors)+geom_vline(data = mean_Phylop,aes(xintercept = mean_Phylop,color=class),linetype="dashed")+theme_classic()
ggsave(paste0(file_dir,"/Phylop_density.pdf"),p4,width = 20,height = 5)
#累积分布图
p5<-ggplot(Phylop_df_all,aes(x=mean_Phylop,col=class))+stat_ecdf(geom = "smooth",linewidth=0.8)+geom_vline(data = mean_Phylop,aes(xintercept = mean_Phylop,color=class),linetype="dashed")+scale_color_manual(values = colors)+theme_classic()
ggsave(paste0(file_dir,"/Phylop_ecdf.pdf"),p5,width = 20,height = 5)


#topN<-10
model_loc_meta<-read.table("/home/fangzj/Workdir/Basset_Unblocked/data/cleandata/flank_seqs/all_classes_flank/best_pt_loc.txt")
scores_df<-"/home/fangzj/Workdir/Basset_Unblocked/data/cleandata/flank_seqs/all_classes_flank/cluster_tab_scores.txt"
cluster_tab_scores<-read.table(scores_df,sep = "\t",header = T)
filter_ids<-sapply(cluster_tab_scores$best_filters,function(x){
  filter_id<-unlist(strsplit(gsub(replacement = "",x = x,pattern = "(rootcluster|subcluster_\\d+)="),";"))
})
names(filter_ids)<-cluster_tab_scores$clsuter_name
filter_ids_sel<-filter_ids
filter_bed_ls<-list()
lapply(names(filter_ids_sel),function(x){
  filters<-filter_ids_sel[[x]]
  for(i in filters){
    print(i)
    prefix<-str_extract(i,".*_split")
    filter<-str_extract(i,"filter\\d+")
    model_out_loc<-model_loc_meta$V4[match(prefix,model_loc_meta$V1)]
    filter_loc<-paste0(model_out_loc,"samples_1000/TF_RLBP_comb/",filter,"_logo.fa")
    h5formodel_loc<-paste0(unlist(strsplit(filter_loc,"\\/"))[1:13],collapse = "/")
    h5formodel<-paste0(h5formodel_loc,"/",dir(h5formodel_loc)[grep("_uniq_act\\.h5",dir(h5formodel_loc))])
    res<-get_seq_loc(fasta = filter_loc,h5formodel = h5formodel)[[1]]
    res$"filter_id"<-i
    filter_bed<-res[,c("chr","seq_start","seq_end","filter_id")]
    filter_bed_ls[[x]][[i]]<<-filter_bed
  }
})
filter_representive_loc<-"/home/fangzj/Workdir/Basset_Unblocked/data/filter_location/"
save(filter_bed_ls,file = paste0(filter_representive_loc,"filter_loc.RData"))
#

lapply(names(filter_bed_ls),function(x){
  print(x)
  tmp <- list.rbind(filter_bed_ls[[x]])
  tmp<-arrange(tmp,chr,seq_start,seq_end)
  bed_file<-paste0(filter_representive_loc,x,".bed")
  write.table(tmp[,1:3],file = bed_file,sep = "\t",row.names = F,col.names = F,quote = F)
  bedtools<-"/home/fangzj/Workdir/software/anaconda/envs/basset/bin/bedtools"
  Phylop_mean_out_motif<-paste0(filter_representive_loc,x,"_motif_phylop.bed")
  code<-paste0(bedtools," map -a ",bed_file," -b ",Phylop_loc," -c 4 -o mean > ",Phylop_mean_out_motif)
  system(code)
  #没跑完
})

topN<-10
filter_ids_sel<-filter_ids[1:topN]
Phylop_top10_motif<-list()
for(i in names(filter_ids_sel)){
  tmp <- read.table(paste0("/home/fangzj/Workdir/Basset_Unblocked/data/filter_location/",i,"_motif_phylop.bed"))
  tmp$V4[tmp$V4=="."]<-0
  tmp$V4<-as.numeric(tmp$V4)
  Phylop_top10_motif[[i]]<-tmp
}
Phylop_top10_motif_df<-list.rbind(Phylop_top10_motif)
