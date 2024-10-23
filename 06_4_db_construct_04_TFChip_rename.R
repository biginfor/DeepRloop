library("stringr")
library("plyr")
library("dplyr")
library("data.table")
#
peaks_loc<-"/path/to/encode/chipsep/narrowpeaks"
TF_chip_meta<-read.table("./data/picked_BPs_downlist_23_1212_1438_fil_meta.tsv",sep = "\t",header = T,quote = "")
peaks<-dir(peaks_loc)
peaks_df<-as.data.frame(cbind(str_replace(peaks,pattern = ".narrowPeak",replacement = ""),peaks))
colnames(peaks_df)<-c("File.accession","file")
peaks_df$"TF"<-mapvalues(x = peaks_df$File.accession,from = TF_chip_meta$File.accession,to = TF_chip_meta$Experiment.target)
peaks_df$TF<-str_replace( peaks_df$TF,pattern= "-human",replacement = "")
#in case of duplicate samples
peaks_df <- peaks_df %>% group_by(TF) %>% dplyr::mutate(rep=row_number())
#
for_jaspar_loc<-"/path/to/encode/chipsep/narrowpeaks"
for( i in peaks_df$File.accession){
  TF_name<-peaks_df[peaks_df$File.accession==i,"TF"]
  new_name<-paste0("Encode-hg38_",TF_name,"_",peaks_df[peaks_df$File.accession==i ,"rep"],".narrowPeak")
  cmd<-paste0("mv ", for_jaspar_loc,peaks_df[peaks_df$File.accession==i,"file"], " ", for_jaspar_loc,new_name)
  system(cmd)
  target_dir<-paste0(for_jaspar_loc,peaks_df[peaks_df$File.accession==i,"TF"])
  if(!dir.exists(target_dir)){
    dir.create(target_dir)
    }
  system(paste0("mv ",for_jaspar_loc,new_name," ",target_dir))
}

peaks_filed_loc<-"./data/Encode_filtered/"
allpeakfiles<-dir(peaks_loc,recursive = T,full.names = T)
chrom_sizes<-read.table("./data/hg38.chrom.sizes")
all_chrs<-chrom_sizes$V1
odd_peak_files<-list()
counter<-0
for(peak_file in allpeakfiles){
  filename<-str_replace(peak_file,pattern = ".*\\/","")
  tmp<-read.table(peak_file)
  norm_chr_idx<-tmp[,1]%in%all_chrs
  odd_chr_idx<-!tmp[,1]%in%all_chrs
  peak_filed<-tmp[norm_chr_idx,]
  if(nrow(tmp[odd_chr_idx,])>0){
  odd_peak_files[[peak_file]]<-tmp[odd_chr_idx,]
  write.table(peak_filed,file = paste0(peaks_filed_loc,filename),sep = "\t",col.names = F,row.names = F,quote = F)
  }
  counter<-counter+1
  print((counter/length(allpeakfiles))*100)
  gc()
}