library("stringr")

#
get_uniq_id<-function(filename){
  file_con <- file(filename,'r')
  tmp<-readLines(file_con)
  gene_IDs<- unlist(lapply(strsplit(tmp[grep("MOTIF",tmp)]," "),function(x){
    ID<-x[length(x)]
    ID<-str_extract(ID,"(?<=\\.)[^.]+$")
    ID
  }))
  close(file_con)
  gene_IDs_uniq<-unique(gene_IDs)
}

#CIS-BP
CISBP_meta_fil<-read.table("./data/CISBP_meta_fil.txt",sep = "\t",comment.char = "",header = T)
gene_IDs_CISBP_uniq<-unique(CISBP_meta_fil$TF_Name)
#JASPAR-human
jaspar_meme_human<-"/path/to/JASPAR_human_meme.txt"
gene_IDs_JASPAR_uniq_human<-get_uniq_id(jaspar_meme_human)
#JASPAR-mouse
jaspar_meme_mouse<-"/path/to/JASPAR_mouse_meme.txt"
gene_IDs_JASPAR_uniq_mouse<-get_uniq_id(jaspar_meme_mouse)
#CIS-BP specific TFs
CISBP_ID_spec <- setdiff(gene_IDs_CISBP_uniq,gene_IDs_JASPAR_uniq_human)
CISBP_uniq_meta <- CISBP_meta_fil[CISBP_meta_fil$TF_Name%in%CISBP_ID_spec,]
TF_freq<-as.data.frame(table(CISBP_uniq_meta$TF_Name))
TF_need2trans<-as.character(TF_freq$Var1[TF_freq$Freq>1])
TF_noneed2trans<-as.character(TF_freq$Var1[TF_freq$Freq==1])
#homogene,human:9606,mouse:10090
library("homologene")
TF_need2trans_mouse<-homologene(TF_need2trans, inTax = 9606, outTax = 10090)$`10090`
matched_genes<- gene_IDs_JASPAR_uniq_mouse[gene_IDs_JASPAR_uniq_mouse%in%TF_need2trans_mouse]
#
unmatched_genes<-TF_need2trans[!TF_need2trans%in%homologene(matched_genes, inTax = 10090, outTax = 9606)$`9606`]
CISBP_meta_fil_for_meme<-CISBP_meta_fil[CISBP_meta_fil$TF_Name%in%c(unmatched_genes,TF_noneed2trans),]
write.table(CISBP_meta_fil_for_meme,"./data/CISBP_meta_fil_for_meme.txt",sep = "\t",quote = F,col.names = T,row.names = F)


get_meme_from_CISBP<-function(meta_info,file_prefix){
  apply(X = meta_info,MARGIN = 1,FUN = function(x){
    tmp_file<-read.table(x[29],header = T,sep = "\t")
    tmp_len<-as.numeric(x[30])
    tmp_id<-x[4]
    tmp_name<-x[7]
    L1<-paste("MOTIF",  tmp_id, tmp_name)
    L2<-paste0("letter-probability matrix: alength= 4 w= ",tmp_len," nsites= 1 E= 0")
    L3<-c()
    for(i in 1:nrow(tmp_file)){
      options(digits = 6)
      L3_1<-paste0("  ",paste0(sprintf("%0.6f",tmp_file[i,2:5]),collapse = "    "),collapse = "")
      L3<-c(L3, L3_1)
      L3
    }
    L4<-paste("URL",x[29],sep = " ")
    ###wite to file 
    out_dir<-"./data/Homo_sapiens_CIS_BP_2023_10_25/"
    if(!dir.exists(out_dir)){dir.create(out_dir)}
    file_path <- paste0(out_dir,file_prefix,"_withoutheader.meme")
    file_conn <- file(file_path,"a")
    lines<-c(L1,L2,L3,L4)
    for( j in 1:2){
      writeLines(lines[j],file_conn)
      if(j == 1 )
        writeLines("",file_conn)
    }
    for( k in 3:(3+tmp_len-1)){
      writeLines(lines[k],file_conn)
    }
    writeLines("",file_conn)
    writeLines(lines[length(lines)],file_conn)
    writeLines("",file_conn)
    close(file_conn)
  })
  file_path <- paste0("./data/Homo_sapiens_CIS_BP_2023_10_25/",file_prefix,"_withoutheader.meme")
  header <- "./data/header_for_tf_chip.txt"
  setwd("./data/Homo_sapiens_CIS_BP_2023_10_25/")
  system(command = paste("cat ",header," ", file_path," > ","./",file_prefix,"_withheader.meme",sep=""))
}
get_meme_from_CISBP(meta_info = CISBP_meta_fil_for_meme,file_prefix = "remained_CISBP_TFs")

##JASPAR_PFMs merge
jaspar_mouse_fil_meme<-"./data/JASPAR_23_11_28/JASPAR_mouse_fil.meme"
JASPAR_mouse_con<-file(jaspar_meme_mouse,"r")
JASPAR_mouse_fil_con<-file(jaspar_mouse_fil_meme,"a")

JASPAR_mouse_meme_edit<-readLines(JASPAR_mouse_con)
wait2write<-NA
del_mark<-NA
for(i in 10:length(JASPAR_mouse_meme_edit)){
  if(grepl("MOTIF",JASPAR_mouse_meme_edit[i])){
    tmp<-JASPAR_mouse_meme_edit[i]
    gene_id<-str_extract(tmp,"(?<=\\.)[^.]+$")
    if(gene_id%in%matched_genes){
      print(gene_id)
      writeLines(text = tmp,con = JASPAR_mouse_fil_con)
    }else{
      del_mark<-i+1
    }
    }
  if(grepl("letter-probability",JASPAR_mouse_meme_edit[i])){
    tmp2<-JASPAR_mouse_meme_edit[i]
    tmp2<-unlist(strsplit(tmp2," "))
    motif_len<-as.numeric(tmp2[which(tmp2=="w=")+1])
    if(i %in% del_mark){
    wait2del<-(i+1):(i+motif_len+2)
    }else{
    writeLines(text = JASPAR_mouse_meme_edit[i],con = JASPAR_mouse_fil_con)
    wait2write<-(i+1):(i+motif_len+2)}
    print(paste0(wait2write[1],"_",wait2write[length(wait2write)]))
  }
  if(i %in% wait2write){
    print(i)
    writeLines(text = JASPAR_mouse_meme_edit[i],con = JASPAR_mouse_fil_con) 
  }
}
close(JASPAR_mouse_con)
close(JASPAR_mouse_fil_con)
