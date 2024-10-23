library("data.table")
library("tidyr")
library('tibble')
library("stringr")
library("rlist")
library("parallel")
#
get_peak_df<-function(actfile){
  save_file<-str_replace(actfile,'\\.txt',"_lnformat.RDS")
  if(file.exists(save_file)){
    peak_sta_df<-readRDS(save_file)
  }else{
  print("fromactfiletoanovarformat_utils start,may cost many time!!!")
  pro_peaks<-fread(actfile)
  pro_peaks<-column_to_rownames(pro_peaks,"V1")
  peak_sta_df_ls<-apply(X = pro_peaks,MARGIN = 1,FUN = function(x){
    sample<-names(x)[x==1]
    sample
  })
  #
  num_cores1<-floor(0.9*detectCores())
  cl1<- makeCluster(num_cores1)
  #parallel::clusterExport(cl1 ,c("peak_sta_df_ls"))
  peak_sta_df_ls_test<-parLapply(cl = cl1,X=names(peak_sta_df_ls),fun = function(x){
    sample<-x
    x<-peak_sta_df_ls[[x]]
    tmp<-matrix(data = cbind(rep(sample,length(x)),x),nrow = length(x),ncol = 2,dimnames = list(NULL,c("loc","sample")))
    tmp
  })
  stopCluster(cl1) 
  gc()
  peak_sta_df<-as.data.frame(list.rbind(peak_sta_df_ls_test))
  peak_sta_df<-separate(data = peak_sta_df,col = loc,into = c("chr","start","end"),sep = ":|\\-")
  peak_sta_df$strand<-str_extract(peak_sta_df$end,"(\\+)|(\\-)|(\\*)",)
  peak_sta_df$end<-str_extract(peak_sta_df$end,"\\d+")
  peak_sta_df$start<-as.numeric(peak_sta_df$start)
  peak_sta_df$end<-as.numeric(peak_sta_df$end)
  saveRDS(peak_sta_df,file=save_file)
  }
  return(peak_sta_df)
}

#actfile<-"/home/fangzj/Workdir/Basset_Unblocked/data/cleandata/flank_seqs/all_classes_flank/both_ext/up_split/1e-90/200_600_85/bin_200_600_85_act.txt"