library("reticulate")
library("data.table")
library("stringr")
library("plyr")
library("dplyr")

#
setwd("/home/fangzj/Workdir/Basset_Unblocked/all_the_needed_files_v2/submit_for_github/")
merged_bam_loc= "./data/merged_bam/"
bamCoverage_loc<-"./data/normliazed_bw/"
if(!dir.exists(paths = merged_bam_loc)){dir.create(path = merged_bam_loc)}
if(!dir.exists(paths = bamCoverage_loc)){dir.create(path = bamCoverage_loc)}


#
final_data_meta<-read.table("./data/sample_meta_info.txt",sep = ",",header = T)
samtools<-"/path/to/samtools"
conda<-"/path/to/conda"
for (i in unique(final_data_meta$tissue)) {
  print(i)
  tmp <-final_data_meta[final_data_meta$tissue==i,]
  GenomeSize<-unique(tmp$eff_genome_size)
  if(nrow(tmp)>1){
  if(!file.exists(paste0(merged_bam_loc,i,"merged.bam"))){
    print("merge duplicated files!")
  system(paste(samtools," merge -f -@ 25 ./data/merged_bam/",i,"_merged.bam ",paste(tmp$bamloc,collapse =   " "),sep = ""))
  system(paste(samtools," index -@ 25 ./data/merged_bam/",i,"_merged.bam",sep = ""))
  system(paste(conda, " run -n deeptools bamCoverage -b ./data/merged_bam/",i,"_merged.bam"," -o ",bamCoverage_loc,i,".bw"," -of bigwig --binSize 50 --normalizeUsing CPM  --minMappingQuality 20 --ignoreDuplicates --effectiveGenomeSize ", GenomeSize, " -p 25 2> ",bamCoverage_loc,i,".log",sep = ""))
  }}else{
    if(!file.exists(paste0(bamCoverage_loc,i,".bw"))){
    bamloc<-tmp$bamloc
      print("just run bamcoverge")
    system(paste(conda, " run -n deeptools bamCoverage -b ",bamloc," -o ",bamCoverage_loc,i,".bw"," -of bigwig --binSize 50 --normalizeUsing CPM --minMappingQuality 20 --ignoreDuplicates --effectiveGenomeSize ", GenomeSize, " -p 25 2> ",bamCoverage_loc,i,".log",sep = ""))
    }else{
    print("bw file exists!")
  }
  }
}

#merge ctrl SRX2187015 SRX2187016
i<-"SHS5Y5_SRP090328_ctrl"
ctrl_ID<-c("SRX2187015","SRX2187016")
ctrl_files<-paste0("/path/to/bamfiles/",ctrl_ID,"/",ctrl_ID,"_hg38",".bam")
system(paste(samtools," merge -f -@ 27 ./data/merged_bam/",i,"_merged.bam ",paste(ctrl_files,collapse =   " "),sep = ""))
system(paste(samtools," index -@ 27 ./data/merged_bam/",i,"_merged.bam",sep = ""))
system(paste(conda, " run -n deeptools bamCoverage -b ./data/merged_bam/",i,"_merged.bam"," -o ",bamCoverage_loc,i,".bw"," -of bigwig --binSize 50 --normalizeUsing CPM --minMappingQuality 20 --ignoreDuplicates --effectiveGenomeSize ", GenomeSize, " -p 27 2> ",bamCoverage_loc,i,".log",sep = ""))


match_info<-read.table("./data/allbam.txt")
plot_bigwig<-function(Acc_ID){
  if(!Acc_ID%in%c("SRX2187015","SRX2187016")){
  bam_file<-match_info[match_info$V2==Acc_ID,"V1"]
  genomsize<-unique(final_data_meta$eff_genome_size[final_data_meta$control==Acc_ID])
  tissue<-unique(final_data_meta$tissue[final_data_meta$control==Acc_ID])
system(paste(conda, " run -n deeptools bamCoverage -b ", bam_file," -o ",bamCoverage_loc,tissue,"-ctrl.bw"," -of bigwig --binSize 50 --normalizeUsing CPM  --minMappingQuality 20 --ignoreDuplicates --effectiveGenomeSize ", genomsize, " -p 27 2> ",bamCoverage_loc,tissue,"-ctrl.log",sep = ""))}
}

for(i in unique(final_data_meta$control)[unique(final_data_meta$control)!=""]){
  print(i)
  plot_bigwig(Acc_ID =i)  
}


#peak calling
peak_loc<-"./data/filter_narrow_peaks/"
load("./data/rlsamples.rda")
final_data_meta$pe_bam<-plyr::mapvalues(final_data_meta$rlsample,from= rlsamples$rlsample,to=rlsamples$paired_end)
final_data_meta$pe_bam[final_data_meta$pe_bam==TRUE]<-" -f BAMPE "
final_data_meta$pe_bam[final_data_meta$pe_bam==FALSE]<-""
merged_bams<-dir(merged_bam_loc,full.names = T)
##This may cost a large amount of time 
for(i in unique(final_data_meta$tissue)){
  print(i)
  treat_file<-merged_bams[grep(pattern =paste0(i,"_merged.bam$"),x = merged_bams)]
  if(length(treat_file)==0){
    print("one treat")
    treat_file<-final_data_meta$bamloc[final_data_meta$tissue==i]
  }else{
    print("multi treat")
  }
  print(treat_file)
  input_file_ID<-unique(final_data_meta[final_data_meta$tissue==i,"control"])
  if(length(input_file_ID)>1){
    print("find merged input")
    input_file<-paste0(merged_bam_loc,i,"_ctrl_merged.bam")
    print(input_file)
  }else{
  input_file<-match_info[match_info$V2==input_file_ID,"V1"]}
  genome_size<-unique(final_data_meta$eff_genome_size[final_data_meta$tissue==i])
  pe_bam <- unique(final_data_meta$pe_bam[final_data_meta$tissue==i])
  print(genome_size)
  print(paste0("paired_end : ",pe_bam))
  if(length(input_file)==0){
    print("no input")
    system(paste0(conda, " run -n macs3 macs3 callpeak -g ", genome_size, pe_bam," -t ", treat_file ," -n ", i," --outdir ", peak_loc," 2> ",peak_loc,i,"_callnarrowpeaks.log",sep = ""))
  }else{
    system(paste0(conda, " run -n macs3 macs3 callpeak -g ", genome_size, pe_bam," -c " ,input_file," -t ", treat_file ," -n ", i," --outdir ", peak_loc," 2> ",peak_loc,i,"_callnarrowpeaks.log",sep = ""))
  print(paste0("input file is: ",input_file))
    }
}

#QC for narrowPeak files
peak_loc<-"./data/filter_narrow_peaks"
blacklist<-"./data/hg38-blacklist.v2.bed"
bedtools<-"/path/to/bedtools"
#
all_files<-dir(peak_loc,full.names = T,recursive = F) 
narrow_peak_file<-all_files[grep("peaks\\.narrowPeak",all_files)]
lapply(X = narrow_peak_file,FUN = function(x){
  tmp<- unlist(strsplit(x,"\\/"))
  tissue<-gsub("_peaks\\.narrowPeak","",tmp[length(tmp)])
  if( file.exists(paste0(peak_loc,"/",tissue,"_peaks_fdr_001_rm_black.narrowPeak"))){
    print(paste0("QC for ", tissue," has been down already"))
  }else{
    print(paste0("QC for ", tissue,"  is going"))
  tmp_peak<-fread(x)
  colnames(tmp_peak)<-c("chrom","chromStart","chromEnd","name","score","strand","signalValue","pValue","qValue","peak")
  #1.fdr<0.01/[-log10(fdr)/qValue>2]
  tmp_peak<-tmp_peak[tmp_peak$qValue>2,]
  #chr1-22,X,Y
  tmp_peak<-tmp_peak[tmp_peak$chrom%in%c(paste("chr",seq(1,22),sep = ""),"chrX","chrY"),]
  fdr_001_file<-gsub(pattern = "\\.",replacement = "_fdr_001\\.",x = x)
  fwrite(x = tmp_peak,file = fdr_001_file,row.names = F,col.names = F,sep = "\t")
  #3.remove black list
  rm_black_file<-gsub(pattern = "\\.",replacement = "_rm_black\\.",x = fdr_001_file)
  command<-paste0(bedtools," intersect -v -a ",fdr_001_file, " -b ", blacklist ," > ",rm_black_file,sep = "")
  system(command = command)}
  #4.discard peaks with a lengh less than 100bps.
  tmp <- fread(x)
  colnames(x)<-c("chr","start","end","name","score","strand","signalValue",'pValue','qValue','peak')
  tmp$width<-tmp$end-tmp$start
  tmp<-as.data.frame(tmp[tmp$width>=100,])
  tmp<-tmp[,-ncol(tmp)]
  QCed_file<-gsub(pattern = "\\.",replacement = "_drop100\\.",x = x)
  fwrite(tmp,file=QCed_file,sep = "\t",quote = F,col.names = F,row.names = F)
})

# E_R_meta for preprocess_features
all_files<-dir(peak_loc,full.names = T,recursive = F) 
E_R_meta_loc<-"./data/"
QC_ed_peaks<-all_files[grep("fdr_001_rm_black",all_files)]
cellline<-unlist(lapply(X = strsplit(QC_ed_peaks,split = "\\/"),FUN = function(x){
   gsub(pattern = "_peaks_fdr_001_rm_black_drop100\\.narrowPeak","",x[length(x)])
}))
E_R_meta<-as.data.frame(cbind(cellline,QC_ed_peaks))
#Exclude outlier HEK293_Frag_SRP277917 and JURAKT
E_R_meta<-E_R_meta[!E_R_meta$cellline%in%c("HEK293_Frag_SRP277917","JURKAT_SRP095885"),]
write.table(E_R_meta,paste0(E_R_meta_loc,"Rloop_E_R_meta.txt",sep=""),col.names = F,row.names = F,quote = F,sep = "\t")

#final sample with bam loc
merged_bam_loc <- "./data/merged_bam"
final_data_meta<-read.table("./data/sample_meta_info.txt",sep = ",",header = T)
bam_treat_idx<-grep(".*(SRP|ERP)\\d+.*(?<!ctrl)_merged.bam$",dir(merged_bam_loc),perl = T)
mergedbam_treat<-dir(merged_bam_loc,full.names = T)[bam_treat_idx]
tissue<-str_replace(string = dir(merged_bam_loc)[bam_treat_idx],pattern = "_merged.bam",replacement = "")
final_data_meta$"treat_bam"<-mapvalues(x = final_data_meta$tissue,from = tissue,to = mergedbam_treat)
final_data_meta$"treat_bam"[-grep("merged.bam",final_data_meta$"treat_bam")]<-final_data_meta$bamloc[-grep("merged.bam",final_data_meta$"treat_bam")]
match_info<-read.table("./data/allbam.txt")
final_data_meta$"ctrl_bam"<-mapvalues(final_data_meta$control,from = match_info$V2,to = match_info$V1)
final_data_meta$"ctrl_bam"[final_data_meta$tissue%in%"SHS5Y5_SRP090328"]<-"./data/merged_bam/SHS5Y5_SRP090328_ctrl_merged.bam"
final_data_meta<-select(final_data_meta,tissue,treat_bam,ctrl_bam)
final_data_meta<-unique(final_data_meta)
write.table(final_data_meta,"./data/cell_bam_match.txt",sep = "\t",col.names = T,row.names = F,quote = F)
