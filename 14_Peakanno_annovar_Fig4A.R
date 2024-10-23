library("plyr")
library("data.table")
library("stringr")
library("tidyr")
library("rlist")
library("ggplot2")
library("dplyr")
source("./myutils/get_seq_loc_utils.R")
source("./fromactfiletoanovarformat_utils.R")
library("GenomicFeatures")
options(scipen = 1000)
conda<-"/path/to/conda"
input_loc<-"/path/to/annovar/local_input/"
output_loc<-"/path/to/annovar/local_output/"
anovar_loc<-"/path/to/annovar/"
humandb<-"/path/to/annovar/humandb/"

##raw peaks
get_ratio_from_annovar<-function(fasta_file_loc,raw_anno=FALSE,pro_anno=FALSE,motif_anno=FALSE,h5formodel,actfile=F,pic_save_loc,peak_sta_df_supp,annovar_out_save=NULL){
if(raw_anno){
  prefix<-"Rloop_raw_peaks"
  peak_loc<-"/path/to/narrowpeaks/"
  get_peak_sta<-function(peak_loc){
    allpeak_list<-list()
    files<-dir(peak_loc,full.names = T)
    files<-files[grep("peaks_fdr_001_rm_black",files)]
    allpeak_list<-lapply(files, function(x){
      sample_name<-unlist(strsplit(x,split = "\\/"))
      sample_name<-sample_name[length(sample_name)]
      tmp<-read.table(x)
      colnames(tmp)<-c("chr","start","end","name","score","strand","sig_Val","p_Val","qVal","peak")
      tmp$sample<-sample_name
      tmp
    })
    peak_sta_df<-list.rbind(allpeak_list)
    peak_sta_df$width<-peak_sta_df$end-peak_sta_df$start+1
    return(list(allpeak_list,peak_sta_df))
  }
  peak_sta_df<-get_peak_sta(peak_loc)[[2]]
  #1-based
  peak_sta_df_1_based<-peak_sta_df
  peak_sta_df_1_based$start<-peak_sta_df_1_based$start+1
  peak_sta_df_1_based$end<-peak_sta_df_1_based$end+1
  peak_sta_df_1_based<-peak_sta_df_1_based[-grep("HEK293_Frag_SRP277917|JURKAT_SRP095885",peak_sta_df_1_based$sample),]
  peak_sta_df_1_based$"REF"<-0
  peak_sta_df_1_based$"ALT"<-0
  peak_sta_df_1_based$chr<-str_replace(string = peak_sta_df_1_based$chr,pattern = "chr",replacement = "")
  peak_sta_df_1_based<-dplyr::select(peak_sta_df_1_based,"chr","start","end","REF","ALT","name","width")
  fwrite(peak_sta_df_1_based,file = paste0(input_loc,prefix,".avinput"),quote = F,sep = "\t",row.names = F,col.names = F)
  get_annovar_res(fasta_file_loc=fasta_file_loc,output_loc = output_loc,prefix = prefix,annovar_out_save = annovar_out_save)
  Rloop_raw_peaks_count<-get_region_ratio(prefix = "Rloop_raw_peaks")[[1]]
  Rloop_raw_peaks_ratio<-get_region_ratio(prefix = "Rloop_raw_peaks")[[2]]
  Rloop_raw_samples_ratio<-get_sample_ratio(Rloop_raw_peaks_count,prefix = "Rloop_raw_peaks")
  #Rloop_raw_samples_ratio
}else{
  Rloop_raw_peaks_count<-get_region_ratio(prefix = "Rloop_raw_peaks")[[1]]
  Rloop_raw_peaks_ratio<-get_region_ratio(prefix = "Rloop_raw_peaks")[[2]]
  Rloop_raw_samples_ratio<-get_sample_ratio(Rloop_raw_peaks_count,prefix = "Rloop_raw_peaks")
  #Rloop_raw_samples_ratio
}
  
##processed peaks
if(pro_anno){
  prefix<-"Rloop_pro_peaks"
  if(length(actfile)!=0){
  peak_sta_df<-get_peak_df(actfile)}else{
    peak_sta_df<-peak_sta_df_supp 
  }
  peak_sta_df_ls_test_2<-lapply(unique(peak_sta_df$sample),FUN = function(x){
    tmp_mat<-peak_sta_df[peak_sta_df$sample==x,]
    tmp_mat$name<-paste0(tmp_mat$sample,paste0("_peak_",seq(1:nrow(tmp_mat))))
    tmp_mat
  })
  peak_sta_df<-as.data.frame(list.rbind(peak_sta_df_ls_test_2))
  peak_sta_df$chr<-str_replace(string = peak_sta_df$chr,pattern = "chr",replacement = "")
  #1-based
  peak_sta_df$start<-peak_sta_df$start+1
  peak_sta_df$end<-peak_sta_df$end+1
  peak_sta_df$"REF"<-0
  peak_sta_df$"ALT"<-0
  peak_sta_df<-dplyr::select(peak_sta_df,"chr","start","end","REF","ALT","name")
  fwrite(peak_sta_df,file = paste0(input_loc,prefix,".avinput"),quote = F,sep = "\t",row.names = F,col.names = F)
  get_annovar_res(fasta_file_loc=fasta_file_loc,output_loc = output_loc,prefix = prefix,annovar_out_save = annovar_out_save)
  Rloop_pro_peaks_count<-get_region_ratio(prefix = "Rloop_pro_peaks")[[1]]
  Rloop_pro_peaks_ratio<-get_region_ratio(prefix = "Rloop_pro_peaks")[[2]]
  Rloop_pro_samples_ratio<-get_sample_ratio(Rloop_pro_peaks_count,prefix = "Rloop_pro_peaks")
  Rloop_pro_samples_ratio
}else{
  Rloop_pro_peaks_count<-get_region_ratio(prefix = "Rloop_pro_peaks")[[1]]
  Rloop_pro_peaks_ratio<-get_region_ratio(prefix = "Rloop_pro_peaks")[[2]]
  Rloop_pro_samples_ratio<-get_sample_ratio(Rloop_pro_peaks_count,prefix = "Rloop_pro_peaks")
  Rloop_pro_samples_ratio
}
  
##motifs
if(motif_anno){
prefix<-"Rloop_motifs"
fasta_file<-dir(fasta_file_loc,full.names = T)[grep("filter\\d+_logo.fa",dir(fasta_file_loc,full.names = T))]
filters_source<-str_extract(fasta_file,pattern = "filter\\d+")
notnone_file <- dir(fasta_file_loc,full.names = T)[grep("filter\\d+_logo.eps",dir(fasta_file_loc,full.names = T))]
filters_sel <- str_extract(notnone_file,pattern = "filter\\d+")
sel_idx<-file.info(fasta_file[filters_source%in%filters_sel])$size!=0
filters_sel<-filters_sel[sel_idx]
#motif_list
motif_feat_list<-list()
motif_feat_list<-lapply(filters_sel,function(x){#
  print(x)
  fa_sel <- fasta_file[filters_source%in%x]
  test<-get_seq_loc(fa_sel,h5formodel = h5formodel)
  test_fa_df<-test[[1]]
  #1-based
  test_fa_df$seq_end<-test_fa_df$seq_end+1
  test_fa_df<-test_fa_df[,c("chr","seq_start","seq_end")]
  test_fa_df$REF<-0
  test_fa_df$ALT<-0
  test_fa_df$name<-x
  colnames(test_fa_df)<-c("chr","start","end","REF","ALT","name")
  test_fa_df
})
motif_feat_df<-list.rbind(motif_feat_list)
fwrite(motif_feat_df,file = paste0(input_loc,prefix,".avinput"),quote = F,sep = "\t",row.names = F,col.names = F)
get_annovar_res(fasta_file_loc=fasta_file_loc,output_loc = output_loc,prefix=prefix,annovar_out_save = annovar_out_save)
Rloop_motifs_count<-get_region_ratio(prefix = "Rloop_motifs")[[1]]
Rloop_motifs_ratio<-get_region_ratio(prefix = "Rloop_motifs")[[2]]
Rloop_motifs_samples_ratio<-get_sample_ratio(Rloop_motifs_count,prefix = "Rloop_motifs")
Rloop_motifs_samples_ratio
}else{
Rloop_motifs_count<-get_region_ratio(prefix = "Rloop_motifs")[[1]]
Rloop_motifs_ratio<-get_region_ratio(prefix = "Rloop_motifs")[[2]]
Rloop_motifs_samples_ratio<-get_sample_ratio(Rloop_motifs_count,prefix = "Rloop_motifs")
Rloop_motifs_samples_ratio
}
  
###plot
ratio_ls<-list(Rloop_raw_peaks_ratio,Rloop_pro_peaks_ratio)
names(ratio_ls)<-c("Rloop_raw_peaks_ratio","Rloop_pro_peaks_ratio")
lapply(names(ratio_ls),function(x){
  p<-ggplot(ratio_ls[[x]],aes(x = "",y =ratio,fill=Var1))+geom_bar(stat = 'identity',color='white',width = 0.5)+
    coord_polar("y",start=0)+theme_void()+theme(legend.position = 'left' )+scale_fill_manual(values = c(
      "genebody"='#53A85F',
      "intergenic"='#808080',
      "UTR3"="#FF00FF",
      "ncRNA_intronic"="#FF0000",
      "intronic"="#8FBC8F",
      "promoter"='#4169E1',
      "upstream"='#4169E1',
      "UTR5"="#7B68EE",
      "terminator"='#F1BB72',
      "downstream"='#F1BB72',
      "ncRNA_exonic"='#00FFFF',
      "exonic"='#DDA0DD',
      "splicing"='#FFD700',
      "promoter_terminator"='#F3B1A0',
      "UTR5;UTR3"='#8B008B',
      "ncRNA_splicing"="#FF8C00",
      "upstream;downstream"="#0000CD"
    ))
  p1<-p+geom_text(aes(label = paste0(round(ratio,2),"%")),position=position_stack(vjust=0.5))
  if(length(pic_save_loc)!=0){
  ggsave(filename = paste0(pic_save_loc,x,".pdf"),p1,width = 10,height = 10)}
})

final_res<-list(Rloop_raw_peaks_count,Rloop_raw_peaks_ratio,Rloop_raw_samples_ratio,Rloop_pro_peaks_count,Rloop_pro_peaks_ratio,Rloop_pro_samples_ratio,Rloop_motifs_count,Rloop_motifs_ratio,Rloop_motifs_samples_ratio)
names(final_res)<-c("Rloop_raw_peaks_count","Rloop_raw_peaks_ratio","Rloop_raw_samples_ratio","Rloop_pro_peaks_count","Rloop_pro_peaks_ratio","Rloop_pro_samples_ratio","Rloop_motifs_count","Rloop_motifs_ratio","Rloop_motifs_samples_ratio")
final_res
}


#annovar
get_annovar_res<-function(fasta_file_loc,output_loc,prefix,annovar_out_save){
script<-paste0(conda," run -n base perl ",anovar_loc,"annotate_variation.pl --outfile ", output_loc,prefix, " --geneanno --dbtype refGene --buildver hg38 --neargene 2000 ",input_loc,prefix,".avinput ",humandb," --chromosome 1-22,X --memtotal 41943040 --thread 30 ")
system(script,intern = TRUE)
##copy the result
raw_out_put<-paste0(output_loc,prefix,".variant_function")
if(!is.null(fasta_file_loc)){
  annovar_out_save<-fasta_file_loc
  }else{
    annovar_out_save<-annovar_out_save
  }
target_loc<-paste0(annovar_out_save,prefix,"/")
if(!dir.exists(target_loc)){
  dir.create(target_loc)
}
target_file<-paste0(target_loc,prefix,".variant_function")
file.copy(raw_out_put,target_file)
}



#result
get_region_ratio<-function(prefix){
  region_hits<-as.data.frame(fread(paste0(output_loc,prefix,".variant_function")))
  region_hits_sta<-as.data.frame(table(region_hits$V1),stringsAsFactors = FALSE)
  match_tab<-as.data.frame(matrix(data = c("exonic","genebody","intronic","genebody","upstream","promoter","downstream", "terminator"),ncol = 2,byrow = T))
  region_hits_sta$"class2"<-mapvalues(x = region_hits_sta$Var1,from = match_tab$V1,to = match_tab$V2)
  region_hits_sta$"class2"[region_hits_sta$"class2"%in%"upstream;downstream"]<-"promoter_terminator"
  all_peaks<-sum(region_hits_sta$Freq)
  region_hits_sta$ratio<-region_hits_sta$Freq/all_peaks*100
  region_hits_sta$Var1<-factor(region_hits_sta$Var1,levels = region_hits_sta$Var1[order(region_hits_sta$ratio,decreasing = T)])
  #
  for_pie_data<-region_hits_sta[,c("class2","ratio")]
  for_pie_data<-aggregate(ratio~class2,for_pie_data,FUN='sum')
  for_pie_data$class2<-factor(for_pie_data$class2,levels = for_pie_data$class2[order(for_pie_data$ratio,decreasing = T)])
  list(region_hits,region_hits_sta)
}


get_sample_ratio<-function(count_file,prefix){
  colnames(count_file)<-c("region","gene_anno","chr","start","end","ALT","REF","peak_ID")
  count_file$class<-str_replace(count_file$peak_ID,pattern = "_peak_\\d+",replacement = "")
  count_file<-dplyr::select(count_file,"region","class")
  region_feat_df_ls<-list()
  for(i in unique(count_file$class)){
    print(i)
    tmp<-count_file[count_file$class==i,]
    tmp<-as.data.frame(table(tmp)/nrow(tmp)*100)
    class<-paste0(prefix,"_",i)
    tmp$"class"<-class
    region_feat_df_ls[[class]]<-tmp
  }
  region_feat_df<-list.rbind(region_feat_df_ls)
}