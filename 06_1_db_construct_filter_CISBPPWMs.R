library("stringr")
#
PWM_files<-dir("/path/to/CISBP_PWMs",full.names = T)
#
CISBP_meta<-read.table("/path/to/CISBP_TF_Information.txt",sep = "\t",check.names = F,comment.char = "",header = T)
id_vec<-c()
for(PWMs in PWM_files){
  id<- regmatches(PWMs, regexpr("M0.*?(?=\\.txt)", PWMs, perl = TRUE))
  id_vec<<-c(id,id_vec)
}

PWM_files_df<-data.frame(PWM_files)
PWM_files_df$Motif_ID<-regmatches(PWM_files_df$PWM_files, regexpr("M\\d+.*?(?=\\.txt)", PWM_files_df$PWM_files, perl = TRUE))
res<-c()
for(i in PWM_files_df$PWM_files){
  nums<-system(command = paste0("tail -1  ", i, "| awk {'print $1'}"),intern = T)
  res<-c(res,nums)
}
PWM_files_df$motif_len<-res
PWM_files_df<-PWM_files_df[grep("\\d+",PWM_files_df$motif_len),]
PWM_files_df$motif_len<-as.numeric(PWM_files_df$motif_len)
#
CISBP_meta_fil<-CISBP_meta[CISBP_meta$Motif_ID%in%PWM_files_df$Motif_ID,]
#
library("plyr")
CISBP_meta_fil$File_loc<-mapvalues(x = CISBP_meta_fil$Motif_ID,from = PWM_files_df$Motif_ID,to = PWM_files_df$PWM_files)
CISBP_meta_fil$motif_len<-mapvalues(x = CISBP_meta_fil$Motif_ID,from = PWM_files_df$Motif_ID,to = PWM_files_df$motif_len)
CISBP_meta_fil$motif_len<-as.numeric(CISBP_meta_fil$motif_len)
write.table(CISBP_meta_fil,"./data/CISBP_meta_fil.txt",sep = "\t",quote = F,col.names = T,row.names = F)