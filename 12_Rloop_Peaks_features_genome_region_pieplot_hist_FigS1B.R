#----raw peak distribution----#
library("rlist")
library(ggplot2)
library(ggforce)
library("GenomicFeatures")
library("data.table")
library("tidyr")
library('tibble')
library("stringr")
peak_loc<-"/home/fangzj/Workdir/Basset_Unblocked/data/Rloop/filter_narrow_peaks"
#
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
  peak_sta_df$width<-peak_sta_df$end-peak_sta_df$start
  return(list(allpeak_list,peak_sta_df))
}
peak_sta_df<-get_peak_sta(peak_loc)[[2]]#
#exclude c("HEK293_Frag_SRP277917","JURKAT_SRP095885")
peak_sta_df<-peak_sta_df[-grep("HEK293_Frag_SRP277917|JURKAT_SRP095885",peak_sta_df$sample),]
#width of raw peaks
draw_bar<-function(peak_sta,range,breaks){
  width<-peak_sta$width
  cuts<-cut(width,breaks = c(seq(0,range,breaks),Inf),labels = c(paste(seq(0,(range-breaks),breaks),seq(breaks,range,breaks),sep="_"),paste0(">",range,"bp")),right = FALSE)
  forpic<-as.data.frame(table(cuts))
  forpic$`Peak count (x10^3)`<-forpic$Freq/1000
  p1<-ggplot(data = forpic, mapping = aes(x = cuts, y = `Peak count (x10^3)`)) + geom_bar(stat = 'identity')+theme_classic()+scale_y_continuous(limits = c(0,max(forpic$`Peak count (x10^3)`)),breaks = seq(0,max(forpic$`Peak count (x10^3)`),50))+theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))
  p1
}
#
p1<-draw_bar(peak_sta = peak_sta_df,range=3000,breaks=100)
