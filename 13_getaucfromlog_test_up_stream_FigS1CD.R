library('stringr')
library('rlist')
library('ggplot2')
library('dplyr')
library('plyr')
library('tidyr')
library('GenomicRanges')
#
top_groups<-c('both_ext_up_split','both_ext_dn_split','dn_ext_up_split','up_ext_dn_split')
mid_groups<-c('up_ext_up_split','dn_ext_dn_split')
bot_groups<-c('dn_ext_both_split','up_ext_both_split','both_ext_both_split')
#
auc_result_loc<-'/path/to/run/log/'
all_files<-dir(auc_result_loc,full.names = T,recursive = T)
auc_result_file<-all_files[grep('quick_checkflank_\\d+.log',all_files)]
auc_result_file_supply<-all_files[grep('aucs_all.log',all_files)]
auc_result_file[match(str_replace(auc_result_file_supply,"[^/]*$",""),str_replace(auc_result_file,"[^/]*$",""))]<-auc_result_file_supply
#
get_auc_res<-function(filename){
  auc_result_con<-file(filename,'r')
  wait_2_extract<-readLines(auc_result_con)
  close(auc_result_con)
  auc_times<-wait_2_extract[grep('current trial best auc',wait_2_extract)]
  aucs<-unlist(lapply(strsplit(auc_times," "),function(x){x[length(x)]}))
  #注意顺序
  auc_mat<-as.data.frame(matrix(data = c(paste0('AUC',seq(1,length(aucs))),aucs),ncol = 2,dimnames = list(NULL,c('rep','value'))))
  auc_mat$value<-as.numeric(auc_mat$value)
  return(auc_mat)
  }
#
counter<-0
aus_res_ls<-lapply(auc_result_file,function(x){
all_fields<-unlist(strsplit(x,'\\/'))
ext_direct<-all_fields[grep('ext',all_fields)]
split_direct<-all_fields[grep('split',all_fields)]
ext_bps<-all_fields[grep('200_600_\\d+',all_fields)]
aucs<-get_auc_res(x)
aucs$'ext_direct'<-ext_direct
aucs$'split_direct'<-split_direct
aucs$'ext_bps'<-ext_bps
counter<<-counter+1
print(paste0('Progress in ',counter/length(auc_result_file)*100,' %'))
return(aucs)
})
auc_res_df<-list.rbind(aus_res_ls)
auc_res_df<-na.omit(auc_res_df)
auc_res_df$ext_bps_nu<-as.numeric(str_replace(auc_res_df$ext_bps,"200_600_",''))
auc_res_df$group<-paste(auc_res_df$ext_direct,auc_res_df$split_direct,sep = "_")
auc_res_df$lm_group<-mapvalues(auc_res_df$group,from = c(top_groups,mid_groups,bot_groups),to = c(rep('top_groups',4),rep('mid_groups',2),rep('bot_groups',3)))
#write.table(auc_res_df,paste0(auc_result_loc,'/important_backup.txt'),sep = "\t",row.names = F,col.names = T,quote = F)
auc_res_df<-read.table(paste0(auc_result_loc,'/important_backup.txt'),header = T)
auc_res_df_group1<-auc_res_df[auc_res_df$lm_group=='top_groups',]
auc_res_group1_sta<-as.data.frame.matrix(table(auc_res_df_group1$ext_bps_nu,auc_res_df_group1$group))
auc_res_group1_sta<-cbind(rownames(auc_res_group1_sta),auc_res_group1_sta)
colnames(auc_res_group1_sta)[1]<-"ext_bps_nu"
write.table(auc_res_group1_sta,paste0(auc_result_loc,'/auc_res_group1_sta.txt'),sep = "\t",quote = F,col.names = T,row.names = F)
#comes from KS1_len_statics_linechart.R
forlm1<-auc_res_df[auc_res_df$lm_group=='top_groups',]
forlm2<-auc_res_df[auc_res_df$lm_group=="mid_groups",]
forlm3<-auc_res_df[auc_res_df$lm_group=="bot_groups",]
get_lm_sum<-function(x,y){
  return(summary(lm(y~x)))
}
get_lm_sum(x = forlm1[forlm1$ext_bps_nu>=85,]$ext_bps_nu,y = forlm1[forlm1$ext_bps_nu>=85,]$value)
get_lm_sum(x = forlm1[forlm1$ext_bps_nu<=85,]$ext_bps_nu,y = forlm1[forlm1$ext_bps_nu<=85,]$value)
get_lm_sum(x=forlm2$ext_bps_nu,y=forlm2$value)
get_lm_sum(x=forlm3$ext_bps_nu,y=forlm3$value)
p1<-ggplot(auc_res_df,aes(x=ext_bps_nu,y=value,group=group))+
  geom_jitter(size = 1.5, width = 0.05, show.legend = FALSE,aes(color=group)) +
  stat_summary(aes(color=group),geom = 'line',fun='mean',cex=1)+
  stat_summary(fun = 'mean',geom = 'point',shape=3,size=1)+
  scale_x_continuous(breaks = seq(0,100,5))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 0,vjust = 0.85,hjust = 0.1))+
  geom_smooth(data = forlm1,mapping = aes(x = ext_bps_nu, y = value,group=lm_group), 
              method = "lm", formula = y ~ x + I((x - 85) * (x > 85)),
              colour = c("blue"),linewidth=0.8,linetype=2)+
  geom_smooth(data = forlm2,mapping = aes(x = ext_bps_nu, y = value,group=lm_group), 
              method = "lm", formula = y ~ x ,
              colour = c("#006400"),linewidth=0.8,linetype=2)+
  geom_smooth(data = forlm3,mapping = aes(x = ext_bps_nu, y = value,group=lm_group), 
              method = "lm", formula = y ~ x ,
              colour = c("red"),linewidth=0.8,linetype=2)
  

#
bartlett.test(value ~ group, data = auc_res_df)
#
aggregate(value~group,data = auc_res_df,FUN = function(x){shapiro.test(x)$p.value})

pdf("./plot/ext_split_qqplot.pdf",height = 10,width = 10)
par(mfrow = c(3,3))
for(i in unique(auc_res_df$group)){
  qqnorm(auc_res_df[auc_res_df$group==i,'value'], ylab="Sample Quantiles", main=i)
  qqline(auc_res_df[auc_res_df$group==i,'value']) 
}
dev.off()
#
kruskal.test(value ~ group, data=auc_res_df)
#
p_res<-with(auc_res_df,pairwise.wilcox.test(value, group, p.adjust.method="holm"))
#
levs<-arrange(aggregate(value ~ group,auc_res_df,median),dplyr::desc(value))$'group'
auc_res_df$group<-factor(auc_res_df$group,levels = levs)
p2<-ggplot(data = auc_res_df, mapping = aes(
  x = group, y = value)) +
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.9))
ggsave("./plot/ext_split_box.pdf",p2,height = 6,width = 8)


#找到最佳的AUC的bp取值处
auc_res_df_forcluster<-auc_res_df[auc_res_df$group%in%top_groups,]
auc_res_df_forcluster$ext_bps<-factor(auc_res_df_forcluster$ext_bps,levels = paste0('200_600_',seq(0,100)))
p3<-ggplot(data = auc_res_df_forcluster, mapping = aes(
  x = ext_bps, y = value,fill=group)) +
  scale_y_continuous(limits = c(0.70,0.75),breaks = seq(0.70,0.75,0.01))+
  geom_boxplot(position=position_dodge(width=1),width=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.9))
ggsave("./plot/best_group_statistics.pdf",p3,height = 9,width = 12)
#
res_cluster_ls<-list()
for(i in as.character(unique(auc_res_df_forcluster$group))){
  tmp<-auc_res_df_forcluster[auc_res_df_forcluster$group==i&auc_res_df_forcluster$ext_bps_nu=='85',]
  res_cluster_ls[[i]]<-tmp[which.max(tmp$value),c('rep','value')]
}
res_cluster_df<-list.rbind(res_cluster_ls)
#generate a tab-split file for motif clustering
gen_meta_tab<-function(groupname,ext_bps_nu,locs,poolmeme=T){
  #
  inputloc='/path/to/forclusterings_input/'
  mirror_loc<-'/path/to/forclusterings_input/'
  docker_loc='/path/to/forclusterings_input/'
  #
  prefix<-'m1d_optrmsprop_drop0.39_lr0.00002_step14_gama0.28_mom0.97_batch512_ks1_12_conv1_300_norm_1e9,0_run_'
  ext<-str_extract(groupname,".*_ext")
  split<-str_replace(groupname,paste0(ext,'_'),'')
  query_meme<-paste0('/path/to/run/log/',ext,'/',split,'/1e-90/','200_600_',ext_bps_nu,'/quick_check/',prefix,locs-1,'/samples_1000/TF_RLBP_comb/filters_meme.txt')
  meme_prefix<-paste0(groupname,'_',ext_bps_nu,"_AUC_",locs)
  query_meme_mirror<-paste0(mirror_loc,meme_prefix,'.meme')
  if(!file.exists(query_meme)){
    #print(query_meme)
    print(paste0('Please run downstream_ana ',meme_prefix," at AUC ",locs))
    }
  if(!dir.exists(query_meme_mirror)){
    system(paste0('cp ',query_meme,' ',query_meme_mirror))
  }
  query_meme_docker<-paste0(docker_loc,meme_prefix,'.meme')
  ifelse(poolmeme,poolmeme<-paste0(docker_loc,'TFs_RLBPs_withheader_IDs.meme'),NULL)
  meta_tab<-matrix(data = c(query_meme_docker,meme_prefix,'meme',poolmeme,'TFs_RLBPs_withheader_IDs','meme'),ncol = 3,byrow = T)
  return(meta_tab)
}
#
meta_tab<-list()
for(i in 1:nrow(res_cluster_df)){
  groupname<-rownames(res_cluster_df[i,])
  ext_bps_nu<-85
  locs<-as.numeric(str_extract(res_cluster_df[i,'rep'],'\\d+'))
  meta_tab[[groupname]]<-gen_meta_tab(groupname = groupname,ext_bps_nu = ext_bps_nu ,locs =locs,poolmeme = T)
}
meta_tab_df<-as.data.frame(unique(list.rbind(meta_tab)))
#
mirror_loc<-'/path/to/forclusterings_input/'
write.table(meta_tab_df,paste0(mirror_loc,'best_auc_with_local_db.tab'),sep = "\t",quote = F,col.names = F,row.names = F)
write.table(meta_tab_df[meta_tab_df$V2!="TFs_RLBPs_withheader_IDs",],paste0(mirror_loc,'best_auc.tab'),sep = "\t",quote = F,col.names = F,row.names = F)
#now goto docker image and run matrix-clusters
