library('stringr')
library('dplyr')
library(ggimage)
library('data.table')
#
cluster_info<-'/path/to/clusters.tab'
compa_info<-'/path/to/pairwise_compa_matrix_descriptions.tab'
tree_info_loc<-'/path/to/motifs_only_trees/'
cor_info<-'/path/to/pairwise_compa.tab'
#
rad_2_ang<-function(radian){
  return( radian / (pi / 180))
}
ang_2_rad<-function(angle){
  return(angle*(pi/180))
}
#
cor_info_tab<-fread(cor_info,fill = T,nThread = 20)
cor_info_tab<-as.data.frame(cor_info_tab[-c(1,length(cor_info_tab)),])
cor_info_tab[1,1]<-'id1'
colnames(cor_info_tab)<-cor_info_tab[1,]
cor_info_tab<-cor_info_tab[-1,]
#
compa_tab<-read.table(compa_info,sep = "\t",header = F)
colnames(compa_tab)<-c('nu','id','name','width','consensus','rc_consensus','IC','nb_sites')
#
cluster_tab<-read.table(cluster_info)
colnames(cluster_tab)<-c('clsuter_name','itmes')
cluster_tab$itmes_nu<-str_count(cluster_tab$itmes,',')+1
#
get_scores<-function(x){
  #Frequency of a cluster in four groups
  ids<-unlist(strsplit(x,","))
  groups<-str_replace(unlist(ids),"_AUC_.*","")
  #The average value of IC,nbsites
  ids_sta<-cbind(groups,ids,compa_tab[compa_tab[,'id']%in%ids,c('IC','nb_sites')])
  mean_IC<-mean(ids_sta$IC)
  mean_nb_sites<-mean(ids_sta$nb_sites)
  #information entropy
  groups_count<-as.numeric(table(ids_sta$groups))
  if(length(groups_count)<4){
    groups_count<-c(groups_count,rep(1e-4))
  }
  p_i <-groups_count/sum(groups_count)
  entropy <- (-sum(p_i * log2(p_i)))
  return(c(mean_IC,mean_nb_sites,entropy))
}
#
socres<-t(sapply(cluster_tab$clsuter_name,function(x){
  #
  tmp_filters<-cluster_tab[cluster_tab$clsuter_name==x,'itmes']
  scores_cluster<-get_scores(tmp_filters)
  filters<-unlist(strsplit(tmp_filters,','))
  cor_tmp<-cor_info_tab[cor_info_tab$id1%in%filters&cor_info_tab$id2%in%filters,]
  #the best number of of cuts ，the most representative motifs
  if(length(filters)<2){
    min_cor_df=data.frame(cluster='rootcluster',best_id=tmp_filters,min_cor=1)
  }else{
  load(paste0(tree_info_loc,'tree_',x,'.Rdata'))
  min_cor<-0
  counter<-1
  while(min_cor<0.6){
    tmp_tree<-cutree(tree,k = counter)
    tmp_tree_df<-data.frame(cluster=paste0('subcluster_',tmp_tree),filter=names(tmp_tree))
    #
    min_cor_df<-aggregate(filter~cluster,tmp_tree_df,function(x){
      if(length(x)>=2){
      compa_tmp<-compa_tab[compa_tab$id%in%x,c('id','IC','nb_sites')]
      res_norm<-list()
      
      for(i in c('IC','nb_sites')){
        res_norm[[i]]<-scale(compa_tmp[,i])
      }
      res_norm_df<-as.data.frame(list.cbind(res_norm))
      res_norm_df$'scores'<-apply(res_norm_df,1,sum)
      res_norm_df$'id'<-compa_tmp$id
      max_id<-res_norm_df$id[which.max(res_norm_df$scores)]
      min_cor<-min(cor_tmp[cor_tmp$id1%in%x&cor_tmp$id2%in%x,'Ncor'])}else{
        max_id<-x
        min_cor<-1
      }
      return(c(max_id,min_cor))
    })
    min_cor_df<-as.matrix(min_cor_df)%>%as.data.frame
    colnames(min_cor_df)[2:3]<-c('best_id','min_cor')
    #
    min_cor<-min(as.numeric(min_cor_df$min_cor))
    counter<-counter+1
  }}
  #
  return(c(scores_cluster,nrow(min_cor_df),paste0(paste(min_cor_df$cluster,min_cor_df$best_id,sep= '='),collapse = ';')))
  }))

cluster_tab[,c('mean_IC','mean_nb_sites','entropy','best_cluster_nu','best_filters')]<-socres
cluster_tab[,3:7]<-apply(cluster_tab[,3:7], 2, as.numeric)

res_norm<-list()
for(i in c('itmes_nu','mean_IC','mean_nb_sites','entropy')){
  res_norm[[i]]<-scale(as.numeric(cluster_tab[,i]))
}
res_norm_df<-as.data.frame(list.cbind(res_norm))
res_norm_df$'scores'<-apply(res_norm_df,1,sum)
cluster_tab$'scores'<-res_norm_df$'scores'
cluster_tab<-cluster_tab%>%arrange(dplyr::desc(scores))
cluster_tab$'ranks'<-seq(1,nrow(cluster_tab))
write.table(cluster_tab,'/path/to/cluster_tab_scores.txt',row.names = F,col.names = T,sep='\t',quote = F)


#######-------画个树型图-------#######
cluster_tab<-read.table('/path/to/cluster_tab_scores.txt',header = T)
library(ggtree)
get_tree_plots<-function(topN=25,clusters=F,gaps_angle=3){
  if(topN){
  clusters=cluster_tab[1:topN,"clsuter_name"]}else if(any(grepl('cluster_\\d+',clusters))){
  clusters=cluster_tab$clsuter_name[cluster_tab$clsuter_name%in%clusters]
  }
#预留缝隙
 gaps<-length(clusters)
 gaps_angle<-gaps_angle
 remained_angle<-360-gaps*gaps_angle
 #
  cluster_tab_df_tmp<-cluster_tab[cluster_tab$clsuter_name%in%clusters,c('clsuter_name','itmes_nu')]
  cluster_tab_df_tmp$'open_angle'<-360-(remained_angle*cluster_tab_df_tmp%>%summarise(itmes_nu/sum(itmes_nu)))
    pic_ls<-list()
    for(i in cluster_tab_df_tmp$clsuter_name){
      #tree_info_loc
      tree_info<-paste0('tree_',i)
      load(paste0(tree_info_loc,tree_info,'.Rdata'))
      assign(tree_info,tree)
      pic_ls[[i]]<-ggtree(get(tree_info),layout = "fan",open.angle =as.numeric(cluster_tab_df_tmp[cluster_tab_df_tmp$clsuter_name==i,'open_angle']),ladderize=FALSE,branch.length="none")
      rm(list=ls()[grep('tree_cluster',ls())])
    }
 merge_pic_data<-lapply(names(pic_ls),function(x){
   tmp<-pic_ls[[x]]
   cbind(x,tmp$data)
   })
 merge_pic_data_df<-list.rbind(merge_pic_data)
 colnames(merge_pic_data_df)[1]<-'cluster_id'
 merge_pic_data_df$cluster_id<-factor(merge_pic_data_df$cluster_id,levels = unique(merge_pic_data_df$cluster_id))
 merge_pic_data_df<-merge_pic_data_df%>%arrange(cluster_id,y)
 levs_cluster<-levels(merge_pic_data_df$cluster_id)
 init_pic_df<-merge_pic_data_df[merge_pic_data_df$cluster_id==levs_cluster[1],]
 init_root_node<-init_pic_df$parent[which(init_pic_df$branch.length==0)]
 #
   for(i in 2:length(levs_cluster)){
   #极坐标与角度
     for(j in c('y','angle')){
       orig<-merge_pic_data_df[merge_pic_data_df$cluster_id==levs_cluster[i],j]
       if(j == 'y'){
         per_gap<-init_pic_df%>%summarise((max(angle)-min(angle))/(max(y)-1))
         orig<-orig+as.numeric(gaps_angle/per_gap)
       }
       if(j == 'angle'){
         orig<-orig+gaps_angle
       }
       merge_pic_data_df[merge_pic_data_df$cluster_id==levs_cluster[i],j]<-max(merge_pic_data_df[merge_pic_data_df$cluster_id==levs_cluster[i-1],j])+orig
       }
   #node起始点
   merge_pic_data_df[merge_pic_data_df$cluster_id==levs_cluster[i],c('parent','node')]<-max(merge_pic_data_df[merge_pic_data_df$cluster_id==levs_cluster[i-1],'node'])+merge_pic_data_df[merge_pic_data_df$cluster_id==levs_cluster[i],c('parent','node')]
   #纵深
   fold<-max(abs(init_pic_df$x))/max(abs(merge_pic_data_df[merge_pic_data_df$cluster_id==levs_cluster[i],'x']))
   merge_pic_data_df[merge_pic_data_df$cluster_id==levs_cluster[i],'x']<- fold*merge_pic_data_df[merge_pic_data_df$cluster_id==levs_cluster[i],'x']
 }
 
p_all<-pic_ls$cluster_4
p_all$data<- merge_pic_data_df[merge_pic_data_df$cluster_id%in%c('cluster_44'),]
p_all$data<- merge_pic_data_df
#注释树状图
cluster_tab_new<-cluster_tab
colnames(cluster_tab_new)[1]<-'cluster_id'
p_all$data$"best_cluster_nu"<-mapvalues(p_all$data$cluster_id,from = cluster_tab_new$cluster_id,to = cluster_tab_new$best_cluster_nu)
match_res<-unlist(strsplit(cluster_tab_new$best_filters,'\\=|;'))
all_best_filters<-match_res[!grepl('rootcluster|subcluster',match_res)]
#
p_all$data$'filter_class'<-sapply(p_all$data$label,function(x){
   if (is.na(x)){return(NA)}
   if(x %in% all_best_filters){
     return('best_filter')
   }else{
     return('redundant')
   }
})
#
load('./data/colours25.Rdata')
p_all$data$label2<-p_all$data$label
p_all$data$label[!grepl('best_filter',p_all$data$filter_class)]<-''
p_all$data$'color_for_dot'<-mapvalues(p_all$data$filter_class,from = c('best_filter','redundant',NA),to = c('red','grey',NA))
color_for_dot<-na.omit(p_all$data$'color_for_dot')
names(colors25)<-cluster_tab_new$cluster_id[1:25]
p_all$data$colour<-mapvalues(p_all$data$cluster_id,from = names(colors25),to =colors25 )#
#
logos_loc<-sapply(p_all$data$label[nchar(p_all$data$label)!=0],function(x){
  prefix<-'m1d_optrmsprop_drop0.39_lr0.00002_step14_gama0.28_mom0.97_batch512_ks1_12_conv1_300_norm_1e9,0_run_'
  tmp<-unlist(strsplit(x,'_'))
  ext_dire<-paste(tmp[1],tmp[2],sep = "_")
  split_dire<-paste(tmp[3],tmp[4],sep = "_")
  ext_bps<-tmp[5]
  ranks<-as.numeric(tmp[7])-1
  filters<-tmp[9]
  eps_loc<-paste0('/path/to/',ext_dire,'/',split_dire,'/1e-90/','200_600_',ext_bps,'/quick_check/',prefix,ranks,'/samples_1000/TF_RLBP_comb/',filters,'_logo.eps')
  #print(file.exists(eps_loc))
  return(eps_loc)
})
p_all$data$logos_loc<-NA
p_all$data$logos_loc[match(names(logos_loc),p_all$data$label)]<-logos_loc

best_logos_loc<-'/path/to/motifs_clustering/logos/'
for(i in unique(p_all$data$cluster_id)){
  cluster_loc<-paste0(best_logos_loc,i,'/')
  ifelse(dir.exists(cluster_loc),'ok',dir.create(cluster_loc))
  for(j in na.omit(p_all$data$logos_loc[p_all$data$cluster_id==i])){
    label<-p_all$data$label[which(p_all$data$logos_loc==j)]
    script_cp<-paste0('cp ',j,' ',cluster_loc)
    system(script_cp)
    str_vec<-unlist(strsplit(j,'\\/'))
    file_name<-str_vec[grep('\\.eps',str_vec)]
    script_rename<-paste0('mv ',cluster_loc,file_name,' ',cluster_loc,label,'.eps')
    system(script_rename)
  }
  
}
#
p1<-p_all+geom_tiplab2(hjust=0)+
  aes(col=cluster_id)+
  scale_colour_manual(values =colors25)+
  geom_tippoint(size=0.5,color=na.omit(p_all$data$color_for_dot))+
  theme(legend.title=element_text(face="bold"), legend.position="left", legend.box="horizontal", legend.text=element_text(size=rel(0.8)))
  #+geom_tiplab(aes(image= logos_loc), geom="image", offset=2, align=F, size=.16, hjust=0) 
ggsave(filename = paste0('/path/to/tree_top',topN,'.pdf'),p1,width = 35,height = 35)
p1
}
get_tree_plots(topN=10)
