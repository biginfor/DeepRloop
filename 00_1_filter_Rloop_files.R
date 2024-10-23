library("dplyr")
library(corrplot)
library("data.table")
library("rlist")
library("dplyr")
library("stringr")
library("tibble")
##
setwd("/home/fangzj/Workdir/Basset_Unblocked/all_the_needed_files_v2/submit_for_github/")

# meta info of samples
load("./data/rlsamples.rda")
# manully annotated sample meta info
anno<-read.csv("./data/rlsamples_ann.csv",header = T)
# QC
rlsamples$decide<-plyr::mapvalues(x = rlsamples$other,from = anno$treatment,to=anno$decide)
sample_sel<-subset(rlsamples,label=="POS"&condition=="S96"&mode=="DRIP"&genotype=="WT"&family=="DRIP"&ip_type=="S9.6"&strand_specific=="FALSE"&moeity=="DNA"&bisulfite_seq=="FALSE"&genome=="hg38"&prediction=="POS"&decide=="pick")
sample_sel_bam<-dplyr::select(sample_sel,rlsample,control)
mark_rep<-dplyr::select(sample_sel,"study","control","rlsample","tissue","eff_genome_size")

#
bam_loc<-"/path/to/bamfiles/"
allfiles<-dir(bam_loc,full.names = T,recursive = T)
bai_file<-allfiles[grep("bai",allfiles)]
bam_file<-substr(x = bai_file,start = 1,stop = nchar(bai_file)-4)
bam_file<-cbind(bam_file,unlist(lapply(X = strsplit(x = bam_file,split = "\\/",perl = T),FUN = function(x){x[length(x)-1]})))
write.table(bam_file,"./data/allbam.txt",quote = F,sep = "\t",row.names = F,col.names = F)


#
mark_rep$bamloc <- plyr::mapvalues(x=mark_rep$rlsample,from = bam_file[,2],to = bam_file[,1])
write.table(mark_rep,"./data/sample_meta_info.txt",col.names = T,row.names = F,quote = F)


#calculate bin reads coverage
system(" bash ./00_2_genome_bin_read_coverage.sh ")

#calculate the similarity between BAM samples 
getcorplot<-function(SRP_ID){
RC_files <- dir("./data/genome_coverage",full.names = T)
file_ID<-str_extract(pattern = "(ERX|SRX)\\d+",string  = RC_files)
RC_files<-unlist(RC_files[file_ID%in%unlist(sample_sel[sample_sel$study==SRP_ID,c("rlsample","control")])])

file_lists<-lapply(RC_files, function(x){
  tmp<-as.data.frame(fread(x))
  tmp$"Seqnames"<-paste(tmp$"V1",tmp$"V2",tmp$"V3",sep  = "_")
  tmp<-tmp[,-c(1,2,3)]
  file_ID<-str_extract(pattern = "(ERX|SRX)\\d+",string  = x)
  colnames(tmp)[1]<-file_ID
  tmp
})
join_fun<-function(df1,df2){
  joined_df <- left_join(df1,df2,by="Seqnames")
  return(joined_df)
}
result <- Reduce(join_fun,file_lists)
result<-column_to_rownames(.data = result,var = "Seqnames")
cor_res<-cor(result)
cor_res_list[[SRP_ID]][["cor"]]<<-cor_res
cor_res_sig<-cor_res
cor_res_sig[cor_res_sig<0.9]<-0
if(!dir.exists("./data/cor_bam/")){dir.create("./data/cor_bam/")}
pdf(file = paste("./data/cor_bam/",SRP_ID,".pdf"),width = 5,height = 5)
if(ncol(cor_res)>3){
corr_test<-cor.mtest(cor_res)
cor_res_list[[SRP_ID]][["p_value"]]<<-corr_test
corrplot(cor_res_sig,method = "color",order = "hclust",p.mat = corr_test$p,sig.level = 0.05,insig = "pch",pch.cex = 0.1,tl.cex = 1)
dev.off()
}else{
pdf(file = paste("./data/cor_bam/",SRP_ID,".pdf"),width = 5,height = 5)
corrplot(cor_res_sig,method = "color",order = "hclust",tl.cex = 1)
dev.off()
}
}
#
multi_sample<-names(table(sample_sel$study)[table(sample_sel$study)>1])
cor_res_list<-list()
for(i in multi_sample){
  print(i)
  getcorplot(SRP_ID = i)
}



#Mark the samples that can be merged together based on the inter-sample correlation coefficient,
#then reanme the samples manually! 
