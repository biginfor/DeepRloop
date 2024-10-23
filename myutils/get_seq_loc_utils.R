library("rhdf5")
library(GenomicRanges)
get_seq_loc<-function(fasta,h5formodel){
  input_file=h5formodel
  test_headers <- h5read(input_file,'test_headers')
  res<-list()
  fafile<-read.table(fasta)
  odd_rows = fafile[seq(1,nrow(fafile),by=2),]
  even_rows = fafile[seq(2,nrow(fafile),by=2),]
  fafile<-cbind(odd_rows,even_rows)
  colnames(fafile)<-c("ID","Seq")
  fafile<-as.data.frame(fafile)
  fafile$"ID"<-gsub(pattern = ">",replacement = "",x = fafile$"ID")
  fafile<-separate(data = fafile, col = ID, into = c("seq_no", "loc_no"), sep = "_")
  fafile$seq_no<-as.numeric(fafile$seq_no)+1#python 和R的idx方式不同
  fafile$loc_no<-as.numeric(fafile$loc_no)+1
  fafile$"chr_region"<-test_headers[fafile$seq_no]
  fafile$"chr_region_bk"<-fafile$"chr_region"
  fafile<-separate(data = fafile, col = chr_region_bk, into = c("chr", "chr_start","chr_end","strand"), sep = ":|\\-|\\(\\+\\)")
  fafile$chr_start<-as.numeric(fafile$chr_start)
  fafile$chr_end<-as.numeric(fafile$chr_end)
  fafile<-fafile[,-ncol(fafile)]
  fafile$chr_start<-fafile$chr_start+1  #custom, 注意和bedtools getfasta 基因组定位方式的不同 #后者忽略了定位的第一个bp
  #注意这里的向下取整来自于motifs_unblocked.py
  filter_length<-unique(nchar(fafile$Seq))
  fafile$"seq_start_tmp"<-fafile$loc_no-floor((filter_length-1)/2)
  fafile$"seq_end_tmp"<-fafile$"seq_start"+filter_length-1
  
  fafile$"seq_start"<-fafile$chr_start+fafile$seq_start_tmp-1
  fafile$"seq_end"<-fafile$seq_start+filter_length-1
  #
  grangeobj<-GRanges(seqnames = fafile$chr,ranges = IRanges(start=fafile$seq_start,end = fafile$seq_end))
  GenomeInfoDb::seqlevels(grangeobj) <- c(paste0("chr",seq(1,22)),"chrX")
  res[[1]]<-fafile
  res[[2]]<-sort(grangeobj)
  return(res)
  #return(sort(grangeobj))
}
