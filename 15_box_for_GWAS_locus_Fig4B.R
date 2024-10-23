library("rlist")
library("dplyr")
library("ggplot2")
library("stringr")
library("RColorBrewer")
library("ggpubr")
library('data.table')
#
heat_data_loc <- "/path/to/prediction_result/"
fig_loc<-str_replace(heat_data_loc,"GWAS_locus","figs_for_locus_test")
binary_data<-dir(heat_data_loc)[!grepl("raw",dir(heat_data_loc))]
#
genomic_region_loc<-"/path/to/functional/genomic/region/"
# binary_downstream等的来源
for( i in dir(genomic_region_loc)){
  names_for_region <- str_replace(i,"_.*","")
  env_name<-paste0("binary_",names_for_region)
  tmp <- binary_data[str_replace(binary_data,"\\.csv","")%in%i]
  assign(env_name,tmp)
}
#
rand_region<-dir("/path/to/random/genomic/region/")
binary_rand<-binary_data[str_replace(binary_data,"\\.csv","")%in%rand_region]
binary_data_new<-c(binary_downstream,binary_exonic,binary_intergenic,binary_intronic,binary_upstream,binary_UTR3,binary_UTR5,binary_rand)


files_ls<-list()
for(files in binary_data_new){
  file_loc<-paste0(heat_data_loc,files)
  tmp<-fread(file_loc,sep = "\t",nThread = 30)
  files_ls[[files]]<- as.data.frame(tmp[,-1])
}

res<-lapply(files_ls,function(x){
  table_ls<-apply(x,2,table)
  if(!is.list(table_ls)){
    names_ls<-colnames(table_ls)
    table_ls<-lapply(1:ncol(table_ls), function(x){table_ls[,x]})
    names(table_ls)<-names_ls
  }
  
  table_ls<-lapply(table_ls,function(y){
    differ<-setdiff(c("0","1"),names(y))
    if(length(differ)!=0){
    tmp1<-0
    names(tmp1)<-differ
    y<-as.table(c(y,tmp1))}
    y<-as.matrix(t(y))
    colnames(y)<-c("close","open")
    y<-cbind(y,apply(y,1,function(x){x[2]/sum(x)}))
    colnames(y)[ncol(y)]<-"access_ratio"
    y
    })
   table_df<-list.rbind(table_ls)
   rownames(table_df)<-names(table_ls)
   table_df<-as.data.frame(table_df)
   table_df
}) %>% list.rbind()

res$"ID"<-rownames(res)
res$"class"<-str_replace(res$ID,"\\.csv.*","")
res$"cell_line"<-str_replace(res$ID,".*\\.csv\\.","")
median_sta<-aggregate(access_ratio ~ class, res, median)
res$class <- factor(res$class,levels = arrange(median_sta,access_ratio)$"class")

all_conditions<-combn(levels(res$class),m = 2)
my_comparisons <- lapply(1:ncol(all_conditions),function(x){ all_conditions[,x]} )
set.seed(1)
p1<-ggplot(data=res,aes(x=class,y=access_ratio))+
  geom_boxplot(aes(fill=class))+
  geom_point(position = position_jitter(0.1))+theme_classic()
  #+scale_fill_manual(values = brewer.pal(length(levels(res$class)),"Set1"))+
  #+stat_compare_means(comparisons = my_comparisons,label="p.signif")
stat_of_data <- as.data.frame(compare_means(access_ratio ~ class, data = res,method = "t.test"))

save_prefix<-"genomic_regions"
#
ggsave(paste0(fig_loc,save_prefix,".pdf"),width = 15,height = 10,p1)
#write.table(stat_of_data,file = paste0(fig_loc,save_prefix,"_test_sta.txt"),sep = "\t",quote = F,row.names = F,col.names = T)
#write.table(res,file = paste0(fig_loc,save_prefix,"plot_data.txt"),sep = "\t",quote = F,row.names = F,col.names = T)