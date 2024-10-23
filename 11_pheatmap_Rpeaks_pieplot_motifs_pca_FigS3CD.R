library("pheatmap")
library("tibble")
library("ggplot2")
library("ggrepel")
library("rlist")
library("dplyr")
library("readxl")
library("plyr")
library("tidyverse")
library(ConsensusClusterPlus)
set.seed(1)
#
colors_<-readxl::read_xlsx('./data/colors.xlsx')
cell_meta_info<-read.table("./data/cell_bam_match.txt",header = T)
auc_result<-read.table('./data//important_backup.txt',header = T)
top_groups<-c('both_ext_up_split','both_ext_dn_split','dn_ext_up_split','up_ext_dn_split')
res_cluster_ls<-list()
for(i in top_groups){
  tmp<-auc_result[auc_result$group==i&auc_result$ext_bps_nu==85,]
  res_cluster_ls[[i]]<-tmp[which.max(tmp$value),c('rep','value')]
}
res_cluster_df<-list.rbind(res_cluster_ls)
#
scores_df<-'/path/to/cluster_tab_scores.txt'
cluster_tab_scores<-read.table(scores_df,sep = "\t",header = T)
load('./data/colours25.Rdata')
cols_for_rank<-colors25
names(cols_for_rank)<-cluster_tab_scores$clsuter_name[1:25]
cols_sec_cluster<-colors25[1:24]
names(cols_sec_cluster)<-c(paste0("subcluster_",seq(1,22)),"rootcluster","not_important")
#
get_score_meta<-function(filter_name){
  prefix<-"(subcluster_\\d+|rootcluster)"
  map_res<-str_extract(cluster_tab_scores$best_filters,paste0("\\b",prefix,"=",filter_name,"\\b"))
  secondary_cluster<-str_extract(na.omit(map_res),prefix)
  if(length(secondary_cluster)==0){
    secondary_cluster<-NA
    map_res<-str_extract(cluster_tab_scores$itmes,paste0("\\b",filter_name,"\\b"))
  }
  row_idx<-which(!is.na(map_res))
  meta_info<-cluster_tab_scores[row_idx,c(1,3:7,9:10)]
  meta_info$"sec_cluster"<-secondary_cluster
  return(meta_info)
}

for(i in 1:nrow(res_cluster_df)){
  groupname<-rownames(res_cluster_df[i,])
  ext_dire<-str_extract(groupname,".*ext")
  split_dire<-str_extract(groupname,"(both|up|dn)_split")
  ext_bps_nu<-85
  locs<-as.numeric(str_extract(res_cluster_df[i,'rep'],'\\d+'))
  prefix<-'m1d_optrmsprop_drop0.39_lr0.00002_step14_gama0.28_mom0.97_batch512_ks1_12_conv1_300_norm_1e9,0_run_'
  pic_save_loc<-paste0("./pic/",groupname,"_AUC_",locs,"_",ext_bps_nu)
  forheat_mat<-read.csv(paste0("/home/fangzj/Workdir/Basset_Unblocked/data/cleandata/flank_seqs/all_classes_flank/",ext_dire,"/",split_dire,"/1e-90/","200_600_",ext_bps_nu,"/quick_check/",prefix,locs-1,"/samples_1000/TF_RLBP_comb_inflout/unit_target_deltas_norm.csv"))
  data<-column_to_rownames(forheat_mat,"X")
  data_fil<-data[apply(data,1,function(x){sum(abs(x))>1e-3*ncol(data)}),]
  colnames(data_fil)[grep("X786.O_ERP120322",colnames(data_fil))]<-'786-O_ERP120322'
  colnames(data_fil)[grep("HeLa_SRP285209",colnames(data_fil))]<-'Hela_SRP285209'
  p1<-pheatmap(t(data_fil),scale = "none")
  #ggsave(filename = paste0(pic_save_loc,"_raw_pic.pdf"),p1,height = 10,width = 20)
  ###
  data_fil_new<-data_fil
  rownames(data_fil_new)<-str_replace(str_replace(rownames(data_fil_new),pattern = "_.*",""),pattern = "^f","filter")
  meta_ls<-list()
  for(i in rownames(data_fil_new)){
  query_item<-paste0(ext_dire,"_",split_dire,"_",ext_bps_nu,"_AUC_",locs,"_m\\d+_",i)
  meta_ls[[i]]<-get_score_meta(query_item)
  }
  meta_for_col<-list.rbind(meta_ls)
  meta_for_col$"color"<-mapvalues(meta_for_col$clsuter_name,from = names(cols_for_rank),to = cols_for_rank)
  meta_for_col$"color"[grep("cluster_",meta_for_col$"color")]<-"white"
  meta_for_col$clsuter_name<-factor(meta_for_col$clsuter_name,levels = cluster_tab_scores$clsuter_name[cluster_tab_scores$clsuter_name%in%unique(meta_for_col$clsuter_name)])
  #
  anno_col<-dplyr::select(meta_for_col,"clsuter_name","sec_cluster")
  anno_col$sec_cluster[is.na(anno_col$sec_cluster)]<-"not_important"
  meta_for_col<-arrange(meta_for_col,clsuter_name)
  anno_col_col<-meta_for_col$color
  names(anno_col_col)<-meta_for_col$clsuter_name
  anno_col_col<-anno_col_col[!duplicated(names(anno_col_col))]
  #
  sub_clusters<-unique(anno_col$sec_cluster)
  all_nu<-sort(as.numeric(na.omit(str_extract(sub_clusters,"\\d+"))))
  anno_col_col_sub<-c(paste0("subcluster_",all_nu),"rootcluster","not_important")
  sec_cluster<-mapvalues(anno_col_col_sub,from = names(cols_sec_cluster),to=cols_sec_cluster)
  names(sec_cluster)<-anno_col_col_sub
  #
  anno_row<-mapvalues(colnames(data_fil_new),from = cell_meta_info$tissue,to = cell_meta_info$category)
  anno_row<-data.frame("category"=anno_row)
  rownames(anno_row)<-colnames(data_fil_new)
  anno_row_col<-mapvalues(unique(anno_row$category),from = colors_$`#单独样本`,to =  colors_$`16进制`)
  names(anno_row_col)<-unique(anno_row$category)
  cols<-list(anno_col_col,anno_row_col,sec_cluster)
  names(cols)<-c("clsuter_name","category","sec_cluster")
  #
  setwd("/home/fangzj/Workdir/Basset_Unblocked/picture/pheatmap/consensus_clust")
  res <- ConsensusClusterPlus(as.matrix(data_fil_new),maxK = 30,reps = 5000,pItem = 0.8,pFeature = 1,clusterAlg = "hc",distance = "euclidean",title = paste0(groupname,"_AUC_",locs,"_",ext_bps_nu),seed = 1,plot = "pdf",writeTable = TRUE)
  icl<-calcICL(res,title = "both_up_85_run3",plot = "pdf")
  best_clust<-res[[11]]
  meta<-as.data.frame(best_clust$consensusClass)
  colnames(meta)<-"groups"
  meta$groups<-mapvalues(meta$groups,from = seq(1,11),to = LETTERS[1:11])
  anno_row$"groups"<-meta$groups
  col_group<-colors25[1:11]
  names(col_group)<-LETTERS[1:11]
  cols$"groups"<-col_group
  p2<-pheatmap(t(data_fil_new),scale = "none",annotation_col =anno_col,annotation_row =  anno_row,annotation_colors = cols)
  ggsave(filename = paste0(pic_save_loc,"_anno_pic.pdf"),p2,height = 10,width = 20)
  #motifs pca
  pca<-prcomp(t(data_fil_new),scale. = F,center = F)
  pca_scores<-as.data.frame(pca$x)
  pca_scores$"ID"<-rownames(pca_scores)
  #
  pca_scores$"batch"<-plyr::mapvalues(x = pca_scores$ID,from = cell_meta_info$tissue,to = cell_meta_info$batch)
  color_sel<-read_xlsx("/home/fangzj/Workdir/Basset_Unblocked/scripts/project/dn_analysis/WGCNA/配色.xlsx")
  col_cell_line<-mapvalues(unique(pca_scores$"batch"),from = color_sel$`#单独样本`,to=color_sel$`16进制`)
  names(col_cell_line)<-unique(pca_scores$"batch")
  pca_scores$mark_labels<-mapvalues(rownames(pca_scores),from = cell_meta_info$tissue,to=cell_meta_info$cell_line)
  pca_scores$"Category"<-mapvalues(rownames(pca_scores),from = cell_meta_info$tissue,to=cell_meta_info$category)
  #
  p3<-ggplot(pca_scores,aes(x=PC1,y=PC2))+geom_point(aes(color=pca_scores$batch))+xlab("PC1")+ylab("PC2")+ggtitle("PCAplot")+scale_color_manual(values = col_cell_line)+theme_classic()+geom_text_repel(aes(label=pca_scores$mark_labels),max.overlaps = 100,size=5)
  #+stat_ellipse(aes(fill=pca_scores$Category),type="norm",geom="polygon",alpha=0.2,color=NA,level = 0.95)
  pca_save_loc<-paste0("/home/fangzj/Workdir/Basset_Unblocked/picture/PCA/infl_pca/",groupname,"_AUC_",locs,"_",ext_bps_nu)
  #ggsave(filename = paste0(pca_save_loc,"_raw_pic.pdf"),p3,width = 8.5,height = 6.5)
  }