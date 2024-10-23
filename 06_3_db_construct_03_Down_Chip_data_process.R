library("readxl")
library("dplyr")
library("rlist")
library("stringr")
library("rtracklayer")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg38")
#downloaded RLBP-TFChip-Seq  source data  ref. PMID:35666183  Table S1 S3

#
RLBPs_all_meta_loc<-"./data/Table_S1.xlsx"
#
RLBPs_IP_MS<-read_xlsx(RLBPs_all_meta_loc,sheet = 1)
RLBPs_IP_MS<-RLBPs_IP_MS[4:nrow(RLBPs_IP_MS),1:3]
colnames(RLBPs_IP_MS)<-c("Cristini_2018_IPMS","Wang_2018_IPMS","Wu_2021_IPMS")
#
RLBPs_Prox_MS<-read_xlsx(RLBPs_all_meta_loc,sheet = 2)
RLBPs_Prox_MS<-RLBPs_Prox_MS[4:nrow(RLBPs_Prox_MS),1:2]
colnames(RLBPs_Prox_MS)<-c("Mosler_ProxMS","Yan_ProxMS")
#
All_RLBPs<-bind_rows(RLBPs_IP_MS,RLBPs_Prox_MS)
All_RLBPs_vec<-unique(unlist(apply(All_RLBPs,2,na.omit)))

# download files
encode_Func_meta_loc<-"./data/experiment_report_2023_12_12_1h_43m.tsv"
encode_proj_meta<-read.table(encode_Func_meta_loc,sep = "\t",skip = 1,header = T,check.names = F,quote = "",comment.char = "")
encode_proj_meta_human_TFChIP<-encode_proj_meta[encode_proj_meta$Organism=="Homo sapiens"&encode_proj_meta$`Assay title`=="TF ChIP-seq",]
RLBP_genes<-All_RLBPs_vec
encode_genes<-c()
for(gene in encode_proj_meta_human_TFChIP$`Target gene symbol`){
  if(length(grep(",",gene))>0){
    gene<-unlist(strsplit(gene,","))
  }else
    gene<-gene
  encode_genes<<-c(encode_genes,gene)  
}
encode_genes<-as.character(na.omit(unique(encode_genes)))
#
inter_BPs<-intersect(encode_genes,RLBP_genes)# number: 276
#setdiff with local db
MEME_info<-read.table('./data/TFs_withheader.meme',fill = T)
MEME_info<-MEME_info[grep("MOTIF",MEME_info$V1),]
genes_in_MEME<-str_replace(MEME_info$V3,pattern = ".*\\.","")
#
drop_BPs<-inter_BPs[inter_BPs%in%genes_in_MEME]
pick_BPs<-inter_BPs
pick_BPs_final<-c(inter_BPs[!inter_BPs%in%genes_in_MEME],"DLX6","ZSCAN29","ATF7","HEY1","NR2C2","RBPJ")

#filter the downlist
meta_info<-read.table("./data/picked_BPs_downlist_23_1212_1438_meta.tsv",sep = "\t",header = T,quote = "")
orign_dn_list<-read.table("./data/picked_BPs_downlist_23_1212_1438.txt")
#check 
table(pick_BPs%in%unique(str_replace(meta_info$Experiment.target,"-human","")))
only_one_format<-as.data.frame(table(meta_info$Experiment.accession)==1)
only_one_format<-meta_info[meta_info$Experiment.accession%in%rownames(only_one_format)[only_one_format[,1]],]
#
doubled_format<-meta_info[(!meta_info$Experiment.accession%in%only_one_format$Experiment.accession)&(meta_info$File.type=="bed"),]
meta_info_filtered<-rbind(doubled_format,only_one_format)
meta_info_filtered<-meta_info_filtered[meta_info_filtered$Experiment.target%in%paste0(pick_BPs_final,"-human"),]
#
library("stringr")
dn_files_id<-str_replace(string = str_replace(string = orign_dn_list$V1,pattern = ".*\\/",replacement = ""),pattern = ".bed.gz|.bigBed",replacement = "")
pick_idx<-dn_files_id%in%meta_info_filtered$File.accession
final_dn_list<-orign_dn_list[pick_idx,]
#
write.table(final_dn_list,file = "./data/picked_BPs_downlist_23_1212_1438_fil.txt",sep = "\t",quote = F,col.names = F,row.names = F)
write.table(meta_info_filtered,file = "./data/picked_BPs_downlist_23_1212_1438_fil_meta.tsv",sep = "\t",quote = F,col.names = T,row.names = F)