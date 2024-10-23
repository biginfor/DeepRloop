library(tidyverse)
library(reshape2)
setwd("/home/fangzj/Workdir/Basset_Unblocked/all_the_needed_files_v2/submit_for_github/")
#read table
meta_info <-read.csv("/home/fangzj/Workdir/Basset_Unblocked/scripts/optuna/all_weights/weights_optuna_bin_200_peak_600_uniqcell_nojurkat_optim_test4/run_metainfo.csv")
meta_info<-meta_info[!is.na(meta_info$AUC),]
#statics
sta<-aggregate(x = meta_info$AUC,by=list(meta_info$KS1_len),FUN="mean")
colnames(sta)<-c("KS1_len","mean_AUC")
p1<-ggplot(meta_info,aes(x=KS1_len,y=AUC))+
  geom_jitter(size = 1.5, width = 0.05, show.legend = FALSE) +
  stat_summary(geom = 'line',fun='mean',linewidth=1)+
  stat_summary(aes(color="mean_point"),fun = 'mean',geom = 'point')+
  stat_summary(geom = 'errorbar',linewidth=0.8,fun.data = 'mean_se',fun.args = list(mult = 1))+
  scale_y_continuous(limits = c(0.68,0.74),expand = c(0,0))+
  scale_x_continuous(breaks = seq(2,50,2))+
  scale_color_manual(values = c(mean_point='#DA4E33'))+
  geom_hline(yintercept = 0.72,colour="blue",linetype="dashed")+
  theme_classic()
if(!dir.exists("./plot")){dir.create("./plot")}
ggsave("./plot/KS1_len_statics.pdf",p1,height = 6,width = 10)
