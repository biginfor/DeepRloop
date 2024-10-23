library("data.table")
library('tidyverse')
library("tibble")
library("pheatmap")
library("stringr")
library("tibble")
library("plyr")
library('logisticPCA')
library("reticulate")
library("ggrepel")
set.seed(1)

#binary-matrix
colors_<-readxl::read_xlsx('./data/colors.xlsx')
clean_data_loc<-"/path/to/clean/data/"
ext_direcs<-c('both_ext')
split_direcs<-c('up_split')
combs<-expand.grid(ext_direcs,split_direcs)
for( i in 1:nrow(combs)){
  print(i)
  ext_direc<-combs[i,1]
  split_direc<-combs[i,2]
  cell_peak<-fread(paste(clean_data_loc,ext_direc,"/",split_direc,"/1e-90/200_600_85/bin_200_600_85_act.txt",sep = ""),check.names = F,header = T)
  cell_peak<-separate(cell_peak,"V1",sep = "[:-]",into = c("chr","start","end"))
  cell_peak$strand<-str_extract(cell_peak$end,"(\\+)|(\\-)|(\\*)",)
  cell_peak$end<-str_extract(cell_peak$end,"\\d+",)
  cell_peak$chr<-factor(cell_peak$chr,levels = c(paste("chr",seq(1,22),sep = ""),c("chrX","chrY")))
  cell_peak$start<-as.numeric(cell_peak$start)
  cell_peak$end<-as.numeric(cell_peak$end)
  cell_peak<-arrange(cell_peak,chr,start,end)
  cell_peak$ID<-paste(cell_peak$chr,cell_peak$start,cell_peak$end,cell_peak$strand,sep="_")
  cell_peak<-cell_peak%>%column_to_rownames("ID")
  cell_peak<-cell_peak[,-c(1,2,3,ncol(cell_peak))]
  #
  ann_row<-as.data.frame(cbind(rownames(cell_peak),NA))
  ann_row$V2<-str_extract(string = ann_row$V1,pattern = "chr(\\d+|[[:alpha:]])")
  ann_row$V2<-factor(ann_row$V2,levels = c(paste("chr",seq(1,22),sep = ""),c("chrX","chrY")))
  ann_row<-column_to_rownames(ann_row,"V1")
  colnames(ann_row)<-"chr"
  #
  idx<-sort(sample(nrow(cell_peak),10000))
  #col
  ann_col<-as.data.frame(cbind(colnames(cell_peak),str_extract(string = colnames(cell_peak),pattern = "(SRP|ERP)\\d+")))
  ann_col<-column_to_rownames(ann_col,"V1")
  colnames(ann_col)<-"batch"
  ann_col$"cell_line"<-str_replace(rownames(ann_col),pattern = "_.*","")
  ann_col[ann_col$cell_line=="HeLa","cell_line"]<-"Hela"
  batch_col<-mapvalues(unique(ann_col$batch),from = colors_$sample,to =  colors_$code)
  names(batch_col)<-unique(ann_col$batch)
  celline_col<-mapvalues(unique(ann_col$cell_line),from = colors_$sample,to =  colors_$code)
  names(celline_col)<-unique(ann_col$cell_line)
  chr_col<-mapvalues(as.character(unique(ann_row$chr)),from = colors_$sample,to =  colors_$code)
  names(chr_col)<-unique(ann_row$chr)
  ann_colors<-list(batch = batch_col,cell_line=celline_col,chr=chr_col)
  #
  getasvdpca(cell_peak = cell_peak,ann_col=ann_col,color_sel = colors_,ext_direc = ext_direc,split_direc = split_direc)
  # #all peaks
  p3<-pheatmap::pheatmap(cell_peak,show_rownames = F,cluster_rows = F,cluster_cols = T,color = c("#FDF5E6","#8B4513"),annotation_row = ann_row,annotation_col = ann_col,annotation_colors = ann_colors)
  save_file<-paste0("./plot/",ext_direc,"_",split_direc,"/")
  ggsave(paste0(save_file,"_default_order_all.png"),p3,width = 16,height = 19)
}

#
getasvdpca<-function(cell_peak,ann_col,color_sel,ext_direc,split_direc,note){
    save_file<-paste0("./plot/logsvd_",ext_direc,"_",split_direc)
  if(file.exists(paste0(save_file,".rds"))){
    logsvd_model_all<-readRDS(paste0(save_file,".rds"))
  }else{
    logsvd_model_all<-logisticSVD(t(cell_peak),k = 2,quiet = F)
    saveRDS(object = logsvd_model_all, file = paste0(save_file,".rds"))
  }
  batch <- ann_col
  batch<- cbind(rownames(batch),batch)
  colnames(batch)[1]<-"sample"
  batch$col<-mapvalues(x = batch$batch,from = color_sel$sample,to=color_sel$code)
  col_pick<-batch$col
  names(col_pick)<-batch$batch
  batch$cell_line<-sub(pattern =  "_.*","",batch$sample) 
  p1<-plot(logsvd_model_all,type="scores")+geom_point(aes(colour = batch$batch))+ggtitle("logisticSVD")+scale_color_manual(values = col_pick)+geom_text_repel(aes(label=batch$cell_line),max.overlaps = 100,size=5)+theme_classic()
  ggsave(paste0(save_file,".pdf"),plot = p1,width = 8.5,height = 6.5)
}
