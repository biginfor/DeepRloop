library('data.table')
library("tidyr")
library("dplyr")
library("plyr")
library("stringr")
library("rlist")
library("tibble")
library("ggplot2")
library("ggpubr")
##ref. doi. 10.1038/s41422-019-0259-z
sel_tissues_1<-c("Skin...Not.Sun.Exposed..Suprapubic.","Skin...Sun.Exposed..Lower.leg.","Lung","Breast...Mammary.Tissue","Kidney...Cortex","Kidney...Medulla")
GTEx_v8_loc<-"/home/fangzj/Workdir/Basset_Unblocked/data/GTEx/GTEx_V8/"
tmp <- read.table(paste0(GTEx_v8_loc,"GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"),skip = 2,header = T,sep = "\t",stringsAsFactors = F)
Name_Description <- tmp[,c("Name","Description")]
rownames(tmp) <- tmp$Name
#scaled_tmp<-cbind(tmp[,c(1,2)],apply(tmp[,c(-1,-2)],2,scale))
scaled_tmp<-cbind(tmp[,c(1,2)],log2(tmp[,c(-1,-2)]+1))
gene_median_tpm <- tmp[,c(-1,-2)]
tissue_gene <- t(gene_median_tpm)

#judge matrix
judge_matrix <- matrix(nrow = length(colnames(gene_median_tpm)),ncol=length(rownames(gene_median_tpm)))
#quantile
vec_quan <- c()
for (tissue_num in seq(1,length(rownames(tissue_gene)),1)) {
  quan=quantile(tissue_gene[tissue_num,],0.9)
  vec_quan <- append(vec_quan,as.numeric(quan))
}

#judge gene and tissue
for (col_i in seq(1,length(colnames(tissue_gene)),1)) {
  is_top5=rank(tissue_gene[,col_i])>=length(rownames(tissue_gene))-(5-1)
  for (row_j in seq(1,length(colnames(gene_median_tpm)),1)) {
    if (is_top5[row_j]) {
      quan=vec_quan[row_j]
      if (tissue_gene[row_j,col_i] > quan) {
        judge_matrix[row_j,col_i]=1
      } else {
        judge_matrix[row_j,col_i]=0
      }
    } else {
      judge_matrix[row_j,col_i]=0
    }
  }
}
judge_matrix_df<-as.data.frame(judge_matrix)
colnames(judge_matrix_df)<-colnames(tissue_gene)
rownames(judge_matrix_df)<-rownames(tissue_gene)
gene_tissue_specific <- as.data.frame(t(judge_matrix_df))
gene_tissue_specific <- gene_tissue_specific[!(rowSums(gene_tissue_specific)==0),]
#
get_a_DEGs<-function(tissue_name){
  tissue_DEGs <- as.data.frame(gene_tissue_specific[,tissue_name])
  colnames(tissue_DEGs) <- c("specific")
  rownames(tissue_DEGs) <- rownames(gene_tissue_specific)
  tissue_DEGs$Name <- rownames(tissue_DEGs)
  tissue_DEGs <- filter(tissue_DEGs,tissue_DEGs$specific == 1)
  tissue_DEGs <- inner_join(tissue_DEGs,Name_Description,by="Name")
  tissue_DEGs <- tissue_DEGs[,-1]
  return(tissue_DEGs)
}
#head(get_a_DEGs(tissue_name = "Skin...Sun.Exposed..Lower.leg."))
#colnames(gene_tissue_specific)[grep("kidney",colnames(gene_tissue_specific),ignore.case = T)]

sel_tissue_DEGs<-list()
for(i in sel_tissues_1){
  sel_tissue_DEGs[[i]]<-get_a_DEGs(tissue_name = i)
}
save(sel_tissue_DEGs,file = paste0(GTEx_v8_loc,"sel_tissue_DEGs.RData"))




#########细胞特异性, no use#########
library("GenomicFeatures")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <-TxDb.Hsapiens.UCSC.hg38.knownGene
seqnames_sel <- c(paste0("chr",seq(1,22)),"chrX")
seqlevels(txdb) <- seqnames_sel
#
#columns(txdb);keytypes(txdb)
#from gene symbol to gene id 
library(org.Hs.eg.db)
library(clusterProfiler)
#
keytypes(org.Hs.eg.db)
load("/home/fangzj/Workdir/Basset_Unblocked/data/random_region_seqs/FGRs.RData")
genic <- genes(txdb)

get_FGRs<-function(ids){
  tmp_gene_id<-str_replace(ids,pattern = "\\.\\d+","")
  mapped_gene_id<-bitr(tmp_gene_id,fromType ='ENSEMBL' ,toType ="ENTREZID",OrgDb = 'org.Hs.eg.db')
  #
  intronic_region<-sel_range$intronic[unlist(sel_range$intronic$gene_id)%in%mapped_gene_id$ENTREZID,]
  exonic_region<-sel_range$exonic[unlist(sel_range$exonic$gene_id)%in%mapped_gene_id$ENTREZID,]
  #
  mapped_tx_name<-AnnotationDbi::select(txdb, keys = mapped_gene_id$ENTREZID, keytype="GENEID", columns="TXNAME")$"TXNAME"
  upstream_region<-sel_range$upstream[sel_range$upstream$tx_name%in%mapped_tx_name,]
  #
  mapped_exon_id<-AnnotationDbi::select(txdb, keys = mapped_gene_id$ENTREZID, keytype="GENEID", columns="EXONID")$"EXONID"
  UTR5_region<-sel_range$UTR5[sel_range$UTR5$exon_id%in%mapped_exon_id,]
  UTR3_region<-sel_range$UTR3[sel_range$UTR3$exon_id%in%mapped_exon_id,]
  downstream_region<-sel_range$downstream[sel_range$downstream$tx_name%in%mapped_tx_name,]
  genic_region_raw <- genic[genic$gene_id%in%mapped_gene_id$ENTREZID,]
  genic_region <- GenomicRanges::reduce(genic_region_raw,ignore.strand=T)
  intergenic_region <- GenomicRanges::gaps(genic_region)
  intergenic_region <- intergenic_region[strand(intergenic_region) == "*"]
  if(intergenic_region@ranges[1]@start==1){
    intergenic_region<-intergenic_region[-1,]
  }else{
    intergenic_region<-intergenic_region
  }
  tmp_res<-list(intronic_region,exonic_region,upstream_region,UTR5_region,UTR3_region,downstream_region,intergenic_region,genic_region_raw)
  names(tmp_res)<-c("intronic","exonic","upstream","UTR5","UTR3","downstream","intergenic","genic_raw")
  return(tmp_res)
}

#get_FGRs
FGRs_res_all<-lapply(X =names(sel_tissue_DEGs),FUN = function(x){
  ids<-sel_tissue_DEGs[[x]][[1]]
  get_FGRs(ids = ids)
} )
names(FGRs_res_all)<-names(sel_tissue_DEGs)
save(FGRs_res_all,file=paste0(GTEx_v8_loc,"FGRs_res_all.RData"))
FGRs_res<-lapply(FGRs_res_all,function(x){x<-x[-length(x)]})
#
save_loc<-"/home/fangzj/Workdir/Basset_Unblocked/data/GTEx/GTEx_V8/bed_files/"
get_bed_fa<-function(tissue_name,FGR_name,save_loc){
  print(FGR_name)
  print(tissue_name)
  grange_wait<-regions[[FGR_name]]
  n_samples<-length(grange_wait)
  wait_mid<-mid(grange_wait)
  wait_start<-wait_mid-300
  wait_end<-wait_mid+300
  region<-as.data.frame(cbind(as.character(seqnames(grange_wait)),wait_start,wait_end))
  colnames(region)<-c("chr","start","end")
  save_name<-paste0(tissue_name,"_",FGR_name,"_",n_samples,".bed")
  fwrite(region,file = paste0(save_loc,save_name),col.names = F,sep = "\t",nThread = 30)
  setwd(save_loc)
  hg38_fa<-"/home/fangzj/Workdir/Basset_Unblocked/data/genomes/hg38/chroms/hg38.fa"
  bedtools <- '/home/fangzj/Workdir/software/anaconda/envs/basset/bin/bedtools'
  get_seqs<-paste0(bedtools,' getfasta -fi ',hg38_fa, " -bed ",save_name, " -fo ", paste0(save_name,'.fa'))
  system(get_seqs)
  gc()
}
#
lapply(names(FGRs_res),FUN = function(x){
  print(x)
  regions<-FGRs_res[[x]]
  lapply(names(regions),FUN = function(y){get_bed_fa(tissue_name = x,FGR_name = y)})
})
###house keeping gene
house_keeping_gene_loc<-"/home/fangzj/Workdir/Basset_Unblocked/data/GeneSets/Housekeeping/HOUNKPE_HOUSEKEEPING_GENES.v2023.2.Hs.gmt"
house_keeping_sets_ls<-clusterProfiler::read.gmt(house_keeping_gene_loc)
house_keeping_sets_1<-clusterProfiler::bitr(house_keeping_sets_ls$gene,fromType = "SYMBOL",toType = "ENSEMBL",OrgDb = 'org.Hs.eg.db')
house_keeping_sets_2<-clusterProfiler::bitr(house_keeping_sets_ls$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = 'org.Hs.eg.db')
house_keeping_sets<-merge(house_keeping_sets_1,house_keeping_sets_2)
save(house_keeping_sets,file=paste0(house_keeping_gene_loc,".RData"))
house_keeping_FGRs<-get_FGRs(ids = house_keeping_sets$ENSEMBL)
house_keeping_FGRs<-house_keeping_FGRs[-length(house_keeping_FGRs)]
save_loc="/home/fangzj/Workdir/Basset_Unblocked/data/GeneSets/Housekeeping/"
lapply(names(house_keeping_FGRs),function(x){
  get_bed_fa(tissue_name = "house_keeping",FGR_name = x,save_loc = save_loc) 
})
#now go to /home/fangzj/Workdir/Basset_Unblocked/scripts/optuna/get_GWAS_out/test_GWAS_locus.py
#c("intronic","exonic","upstream","UTR5","UTR3","downstream","intergenic")
best_model_loc = '/home/fangzj/Workdir/Basset_Unblocked/data/cleandata/flank_seqs/all_classes_flank/both_ext/up_split/1e-90/200_600_85/quick_check/m1d_optrmsprop_drop0.39_lr0.00002_step14_gama0.28_mom0.97_batch512_ks1_12_conv1_300_norm_1e9,0_run_3/GWAS_locus/'
allfiles<-dir(best_model_loc,full.names = T)###严重警告！！在for循环中，dir读取不全？？？
#FGR<-"UTR5"
pic_loc<-"/home/fangzj/Workdir/Basset_Unblocked/picture/GTEx_out/accessibility_ratio/"
get_predict_heat<-function(FGR){
  load("/home/fangzj/Workdir/Basset_Unblocked/data/GTEx/GTEx_V8/meta_info_match.RData")
  FGRs_sel<-paste0("(",paste(sel_tissues_1,collapse = "|"),")","_",FGR,"_\\d+.bed.fa.csv")
  #print(FGRs_sel)
  all_FGRs_predict<-allfiles[grep(FGRs_sel,allfiles)]
  #print(grep(FGRs_sel,dir(best_model_loc),perl = T))
  #print(all_FGRs_predict)
  FGR_predict_ls<-list()
  for(predict_res in all_FGRs_predict){
    tissue_name<-tail(unlist(str_split(predict_res,"\\/")),n=1)
    tissue_name<-str_replace(tissue_name,pattern = "_.*","")
    test_out<-read.table(predict_res)
    test_out<-test_out[,meta_info$cell_line]
    sta_test_out<-as.data.frame(apply(test_out,2,sum)/nrow(test_out))
    sta_test_out<-cbind(rownames(sta_test_out),sta_test_out)
    colnames(sta_test_out)<-c("cell_line","accessibility_ratio")
    FGR_predict_ls[[tissue_name]]<-sta_test_out
  }
  save(FGR_predict_ls,file = paste0("/home/fangzj/Workdir/Basset_Unblocked/data/GTEx/GTEx_V8/FGR_predict_ls_",FGR,".RData"))
  FGR_predict_df<-list.rbind(FGR_predict_ls)
  FGR_predict_df$"tissue_source"<-str_replace(rownames(FGR_predict_df),pattern = FGR_predict_df$cell_line,replacement = "")
  #print(FGR)
  FGR_predict_df_wide<-pivot_wider(data = FGR_predict_df,names_from = cell_line,values_from = "accessibility_ratio")
  FGR_predict_df_wide<-column_to_rownames(FGR_predict_df_wide,var = "tissue_source")
  FGR_predict_df_wide<-t(FGR_predict_df_wide)
  row_anno<-column_to_rownames(meta_info,"cell_line")
  p1<-pheatmap::pheatmap(FGR_predict_df_wide,scale = "column",main = FGR,annotation_row = row_anno)
  ggsave(paste0(pic_loc,FGR,"accessibility_ratio_heat.pdf"),p1,width = 7,height = 6)
  #
  lapply(names(FGR_predict_ls),function(x){
    #print(paste0(FGR,"_",x))
    tmp<-FGR_predict_ls[[x]]
    tmp$"matched_tissue"<-mapvalues(tmp$cell_line,from = meta_info$cell_line,to=meta_info$matched_tissue)
    my_comparisons<-list(c("Breast","Kidney"),c("Breast","Skin"))
    
    p<-ggplot(tmp,aes(x=matched_tissue,y=accessibility_ratio))+geom_boxplot(width=0.2,position = position_dodge(0.9),outlier.size = 0)+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",size=5)+geom_jitter(position = position_dodge(width = 0.1))+theme_classic()+theme(text = element_text(size = 15))+ggtitle(paste0(FGR,"_",x))
    ggsave(paste0(pic_loc,FGR,"_",x,"_accessibility_ratio_box.pdf"),p,width = 6,height = 6)
  })
}
all_FGRs<-c("intronic","exonic","upstream","UTR5","UTR3","downstream","intergenic")
for(FGR in all_FGRs){
  #print(FGR)
  get_predict_heat(FGR)
}
#house keeping
house_keeping_prediction<-allfiles[grep("house_keeping.*\\.bed\\.fa\\.csv",allfiles)]
FGR_predict_ls_house_keeping<-list()
for(i in house_keeping_prediction){
  FGR<-str_extract(string = i,pattern = "house_keeping_.*_\\d+")
  tmp<-read.table(i)
  tmp<-tmp[,meta_info$cell_line]
  sta_tmp<-as.data.frame(apply(tmp,2,sum)/nrow(tmp))
  sta_tmp<-cbind(rownames(sta_tmp),sta_tmp)
  colnames(sta_tmp)<-c("cell_line","accessibility_ratio")
  FGR_predict_ls_house_keeping[[FGR]]<-sta_tmp
}
FGR_predict_df_house_keeping<-list.rbind(FGR_predict_ls_house_keeping)
FGR_predict_df_house_keeping$"tissue_source"<-str_replace(rownames(FGR_predict_df_house_keeping),pattern = FGR_predict_df_house_keeping$cell_line,replacement = "")

FGR_predict_df_wide_house_keeping<-pivot_wider(data = FGR_predict_df_house_keeping,names_from = cell_line,values_from = "accessibility_ratio")
FGR_predict_df_wide_house_keeping<-column_to_rownames(FGR_predict_df_wide_house_keeping,"tissue_source")
FGR_predict_df_wide_house_keeping<-t(FGR_predict_df_wide_house_keeping)
row_anno<-column_to_rownames(meta_info,"cell_line")
p1<-pheatmap::pheatmap(FGR_predict_df_wide_house_keeping,scale = "column",annotation_row = row_anno)
ggsave(filename = "/home/fangzj/Workdir/Basset_Unblocked/picture/Housekeeping_out/accessibility_ratio/house_keeping_FGR_scale_FGR.pdf",p1,width = 10,height = 8)
p2<-pheatmap::pheatmap(FGR_predict_df_wide_house_keeping,scale = "row",annotation_row = row_anno)
ggsave(filename = "/home/fangzj/Workdir/Basset_Unblocked/picture/Housekeeping_out/accessibility_ratio/house_keeping_FGR_scale_sample.pdf",p2,width = 10,height = 8)

#和act_file取交集
library(GenomicFeatures)
load(paste0(GTEx_v8_loc,"FGRs_res_all.RData"))
tissue_genic<-lapply(FGRs_res_all,function(x){
  return(x[["genic_raw"]])
})
ext_dire="both_ext";split_dire="up_split";ext_bps_nu="85"
file_dir<-paste0('/home/fangzj/Workdir/Basset_Unblocked/data/cleandata/flank_seqs/all_classes_flank/',ext_dire,'/',split_dire,'/1e-90/','200_600_',ext_bps_nu)
#get_Phylop_mean<-function(ext_dire,split_dire,ext_bps_nu)
query_file<-paste0(file_dir,'/bin_200_600_85_act.txt')
act_file<-read.table(query_file,sep="\t",header = T,row.names = 1)
act_file_sel<-act_file[,grep("Fib|IMR90|SRP101967|HEK293|MCF10A",colnames(act_file))]
genic_peak_ratio<-list()
for(i in colnames(act_file_sel)){
  print(i)
  peak_region<-rownames(act_file_sel)[act_file_sel[,i]==1]
  peak_region_df<-list.rbind(strsplit(peak_region,split = ":|-|\\(\\+\\)"))%>%as.data.frame()
  colnames(peak_region_df)<-c("chr","start","end")
  peak_range<-makeGRangesFromDataFrame(peak_region_df)
  tmp_ratio<-lapply(tissue_genic,function(x){
    overlaps<-findOverlaps(query = x,subject = peak_range,minoverlap = 300)
    #拿each组织特异性DEG去衡量细胞系样本（细胞系间可比）
    query_nu<-length(x)
    subject_nu<-length(peak_range)
    Hits_nu<-length(overlaps)
    #
    subjectHits_ratio<-Hits_nu/subject_nu #R-loop peaks上，each tissue-DEG的分布比例,cell-line间可比
    #
    queryHits_ratio<-Hits_nu/query_nu #平均1个tissue-DEG 区域对应着几条peaks，tissue-DEGs间可比
    #
    return(c(query_nu,subject_nu,Hits_nu,subjectHits_ratio,queryHits_ratio))
  })
  tmp_ratio<-list.rbind(tmp_ratio)
  tmp_ratio<-cbind(names(tissue_genic),unlist(tmp_ratio),i)%>%as.data.frame()
  colnames(tmp_ratio)<-c("tissue","query_nu","subject_nu","Hits_nu","subjectHits_ratio","queryHits_ratio","cell_line")
  tmp_ratio[,2:6]<-apply(tmp_ratio[,2:6],2,as.numeric)
  genic_peak_ratio[[i]]<-tmp_ratio
}

genic_peak_ratio_df<-list.rbind(genic_peak_ratio)
meta_info<-as.data.frame(unique(genic_peak_ratio_df$cell_line))
colnames(meta_info)[1]<-"cell_line"
meta_info$"matched_tissue"<-NA
meta_info$"matched_tissue"[grep("SRP101967|MCF10A",meta_info$cell_line)]<-"Breast"
meta_info$"matched_tissue"[grep("SRP041718",meta_info$cell_line)]<-"Skin"
meta_info$"matched_tissue"[grep("HEK293",meta_info$cell_line)]<-"Kidney"
meta_info$"matched_tissue"[grep("IMR90",meta_info$cell_line)]<-"Lung"
save(meta_info,file = "/home/fangzj/Workdir/Basset_Unblocked/data/GTEx/GTEx_V8/meta_info_match.RData")
#
genic_peak_ratio_df$"match_tissue"<-mapvalues(genic_peak_ratio_df$cell_line,from = meta_info$cell_line,to = meta_info$matched_tissue)
#
pic_loc2<-"/home/fangzj/Workdir/Basset_Unblocked/picture/GTEx_out/Hits_ratio/"
pic_ls_cellines<-list()
for(i in unique(genic_peak_ratio_df$tissue)){
  box_for_cell_line<-genic_peak_ratio_df[genic_peak_ratio_df$tissue==i,]
  my_comparisons<-list(c("Breast","Kidney"),c("Breast","Skin"))
  p<-ggplot(box_for_cell_line,aes(x=match_tissue,y=subjectHits_ratio))+geom_boxplot(width=0.2,position = position_dodge(0.9),outlier.size = 0)+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",size=5)+geom_jitter(position = position_dodge(width = 0.1))+theme_classic()+theme(text = element_text(size = 15))+ggtitle(i)
  ggsave(paste0(pic_loc2,"/",i,"_subjectHits_ratio.pdf"),p,height = 6,width = 7)
  pic_ls_cellines[[i]]<-p
}

#
queryHits_df<-genic_peak_ratio_df[,c("tissue","queryHits_ratio","cell_line")]
DEGs_heat<-pivot_wider(
  names_from = cell_line,
  values_from = queryHits_ratio,
  data = queryHits_df
)
DEGs_heat<-column_to_rownames(DEGs_heat,"tissue")

col_anno<-column_to_rownames(meta_info,"cell_line")
p2<-pheatmap::pheatmap(DEGs_heat,scale = "column",annotation_col = col_anno)
p3<-pheatmap::pheatmap(DEGs_heat,scale = "row",annotation_col = col_anno)
ggsave(paste0(pic_loc2,"queryHits_ratio_scaled_on_sample.pdf"),p2,height = 6,width = 8)
ggsave(paste0(pic_loc2,"queryHits_ratio_scaled_on_DEG.pdf"),p3,height = 6,width = 8)
#########细胞特异性, no use#########



#gene Exp 
#####no use#####
sel_tissues_1<-c("Skin...Not.Sun.Exposed..Suprapubic.","Skin...Sun.Exposed..Lower.leg.","Lung","Breast...Mammary.Tissue","Kidney...Cortex","Kidney...Medulla")
GTEx_v8_loc<-"/home/fangzj/Workdir/Basset_Unblocked/data/GTEx/GTEx_V8/"
GTEx <-fread(paste0(GTEx_v8_loc,"GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct"),skip = 2,header = T,sep = "\t",nThread = 30)
GTEx<-as.data.frame(GTEx)
save(GTEx,file=paste0(GTEx_v8_loc,"GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.Rdata"))
load(paste0(GTEx_v8_loc,"GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.Rdata"))
meta_info1<-read.table(paste0(GTEx_v8_loc,'GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'),header = T,sep = "\t",quote = '')
#####no use#####


#########gene exp############
genesymbol_by_tissue_exp<-tmp[,-1]
genesymbol_by_tissue_exp<-genesymbol_by_tissue_exp[,c("Description",sel_tissues_1)]
genesymbol_by_tissue_exp_mean<-as.data.frame(cbind(genesymbol_by_tissue_exp$Description,apply(genesymbol_by_tissue_exp[,-1],1,mean)))
colnames(genesymbol_by_tissue_exp_mean)<-c("gene_symbol","mean_Exp")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <-TxDb.Hsapiens.UCSC.hg38.knownGene
genic <- genes(txdb)
load("/home/fangzj/Workdir/Basset_Unblocked/data/GTEx/GTEx_V8/meta_info_match.RData")
ext_dire="both_ext";split_dire="up_split";ext_bps_nu="85"
file_dir<-paste0('/home/fangzj/Workdir/Basset_Unblocked/data/cleandata/flank_seqs/all_classes_flank/',ext_dire,'/',split_dire,'/1e-90/','200_600_',ext_bps_nu)
#get_Phylop_mean<-function(ext_dire,split_dire,ext_bps_nu)
query_file<-paste0(file_dir,'/bin_200_600_85_act.txt')
act_file<-fread(query_file,sep="\t",header = T,nThread = 25)%>%as.data.frame()
act_file<-column_to_rownames(act_file,"V1")
act_file<-act_file[,meta_info$cell_line]
peak_count<-as.data.frame(apply(act_file,1,sum))
peak_count <- cbind(rownames(peak_count),peak_count)
colnames(peak_count)<-c("region","count")
peak_count<-peak_count[!peak_count$count==0,]
peak_count<-dplyr::arrange(peak_count,desc(count))
peak_count<-separate(peak_count,col = "region",into = c("chr","start","end","_"),sep = ":|-|\\(\\+\\)")
peak_count<-peak_count[,-4]
#
peak_range_for_count<-makeGRangesFromDataFrame(peak_count[,1:3])
overlaps<-findOverlaps(query = peak_range_for_count,subject = genic,minoverlap = 300)%>%as.data.frame()
overlaps$"gene_id"<-genic$gene_id[overlaps$subjectHits]
overlap_result_df<-overlaps%>%group_by(queryHits)%>%dplyr::summarise(gene_ids=paste(gene_id,collapse = ";"))
peak_count$"gene_ids"<-NA
peak_count$"gene_ids"[overlap_result_df$queryHits]<-overlap_result_df$gene_ids
#
overlaps$"peak_id"<-rownames(peak_count)[overlaps$queryHits]
mapped_gene<-clusterProfiler::bitr(overlaps$gene_id,fromType ="ENTREZID",toType ="SYMBOL",OrgDb = 'org.Hs.eg.db')
overlaps$"gene_symbol"<-mapvalues(overlaps$gene_id,from = mapped_gene$ENTREZID,to = mapped_gene$SYMBOL)
overlaps<-overlaps[overlaps$gene_symbol%in%genesymbol_by_tissue_exp_mean$gene_symbol,]
overlaps<-left_join(x = overlaps,y = genesymbol_by_tissue_exp,by = c("gene_symbol" = "Description"),relationship = "many-to-many")
overlaps$"mean_exp"<-rowMeans(overlaps[,colnames(overlaps)%in%sel_tissues_1])
#dropped
#scaled_exp<-apply(overlaps[,colnames(overlaps)%in%sel_tissues_1],2,scale)
#colnames(scaled_exp)<-paste0("scaled_",colnames(scaled_exp))
#dropped
scaled_exp<-scaled_tmp[match(overlaps$gene_symbol,scaled_tmp$Description),sel_tissues_1]
colnames(scaled_exp)<-paste0("scaled_",colnames(scaled_exp))
overlaps<-cbind(overlaps,scaled_exp)
overlaps$"scaled_mean_exp"<-rowMeans(overlaps[,colnames(overlaps)%in%colnames(scaled_exp)])
fold_change_exp<-overlaps[,colnames(overlaps)%in%sel_tissues_1]/apply(genesymbol_by_tissue_exp[,sel_tissues_1],2,mean)
colnames(fold_change_exp)<-paste0("fold_chg_",colnames(fold_change_exp))
overlaps<-cbind(overlaps,fold_change_exp)
overlaps$"mean_exp"<-genesymbol_by_tissue_exp_mean$mean_Exp[match(overlaps$gene_symbol,genesymbol_by_tissue_exp_mean$gene_symbol)]
#
overlaps$"mean_foldchg_exp"<-rowMeans(overlaps[,colnames(overlaps)%in%colnames(fold_change_exp)])
overlaps$"peak_count"<-peak_count$count[match(overlaps$peak_id,rownames(peak_count))]
save(overlaps,peak_count,file = paste0(GTEx_v8_loc,"overlaps_and_peakcount.RData"))

GTEx_v8_loc<-"/home/fangzj/Workdir/Basset_Unblocked/data/GTEx/GTEx_V8/"
load(paste0(GTEx_v8_loc,"overlaps_and_peakcount.RData"))

test_overlaps<-overlaps[overlaps$peak_count>=0,c(paste0("scaled_",sel_tissues_1),"peak_count")]
test_overlaps<-test_overlaps%>%dplyr::group_by(peak_count)%>%dplyr::summarise(across(everything(),mean))
test_overlaps<-pivot_longer(data = test_overlaps,cols = starts_with("scaled_"),names_to = "tissue_mean_exp",values_to = "scaled_exp")
test_peak_mean_exp<-aggregate(scaled_exp~peak_count,mean,data=test_overlaps)
test_overlaps$tissue_mean_exp<-factor(test_overlaps$tissue_mean_exp,levels = c('scaled_Kidney...Medulla','scaled_Lung','scaled_Breast...Mammary.Tissue','scaled_Kidney...Cortex','scaled_Skin...Sun.Exposed..Lower.leg.','scaled_Skin...Not.Sun.Exposed..Suprapubic.'))
colors<-c("black","#FFBE7A","purple","red","#82B0D2","#2E8B57")
p1<-ggplot(data = test_overlaps,aes(x=peak_count,y=scaled_exp))+geom_boxplot(aes(factor(peak_count)),outlier.size = 0,outlier.color = NA)+geom_jitter(width = 0.1,size=1,aes(color=tissue_mean_exp))+scale_color_manual(values = colors)+geom_smooth(method = "lm",formula = y~x,color="red",fill="pink")+stat_cor(digits = 4,method = "pearson",size=5)+theme_classic()+theme(text = element_text(size = 20))
ggsave(filename = "/home/fangzj/Workdir/Basset_Unblocked/picture/GTEx_out/Exp_box_cor/scaled_mean_exp_with_tissue.pdf",p1,width = 12,height = 6)
p1<-ggplot(data = test_peak_mean_exp,aes(x=peak_count,y=scaled_exp))+geom_point()+geom_smooth(method = "lm",formula = y~x,color="red",fill="pink")+stat_cor(digits = 4,method = "pearson",size=5)+theme_classic()+theme(text = element_text(size = 20))
ggsave(filename = "/home/fangzj/Workdir/Basset_Unblocked/picture/GTEx_out/Exp_box_cor/scaled_mean_exp.pdf",p1,width = 12,height = 6)


#one class of tissue with house keeping and random genes
GTEx_v8_loc<-"/home/fangzj/Workdir/Basset_Unblocked/data/GTEx/GTEx_V8/"
load(paste0(GTEx_v8_loc,"overlaps_and_peakcount.RData"))
house_keeping_gene_loc<-"/home/fangzj/Workdir/Basset_Unblocked/data/GeneSets/Housekeeping/HOUNKPE_HOUSEKEEPING_GENES.v2023.2.Hs.gmt"
load(paste0(house_keeping_gene_loc,".RData"))
house_keeping_sets<-house_keeping_sets[house_keeping_sets$SYMBOL%in%(scaled_tmp$Description),"SYMBOL"]
#tissue<-"Kidney...Cortex"
#judge_type<-"scaled_exp"
#iterations<-1000
sel_cols<-c("queryHits","subjectHits","gene_id","peak_id","gene_symbol","mean_exp","scaled_mean_exp","mean_foldchg_exp","peak_count")
get_one_tissue_box<-function(tissue,judge_type,iterations){
  if(judge_type=="scaled_exp"){
    tissue_exp<-overlaps[,c(sel_cols,paste0("scaled_",tissue))]
    p1<-ggplot(data = tissue_exp,aes(x=peak_count,y=tissue_exp[,paste0("scaled_",tissue)]))+geom_boxplot(aes(factor(peak_count)),outlier.size = 0,outlier.color = NA)+geom_jitter(width = 0.1,size=1)+geom_smooth(method = "lm",formula = y~x,color="red",fill="pink")+stat_cor(digits = 4,method = "pearson",size=5)+theme_classic()+theme(text = element_text(size = 20))
    ggsave(filename = paste0("/home/fangzj/Workdir/Basset_Unblocked/picture/GTEx_out/Exp_box_cor/with_random_housekeeping/",tissue,"_",judge_type,"_",iterations,"box.png"),p1,width = 8,height = 5)
    # #
    # idx_match<-cbind(rownames(tissue_exp),tissue_exp$peak_count)
    # get_random<-function(x){
    #   if(length(x)>=3000){
    #     set.seed(4)
    #     return(sample(x,3000))
    #   }
    # }
    # random_res<-as.data.frame(aggregate(V1~V2,get_random,data=idx_match))
    # sel_idx<-as.numeric(unlist(random_res$V1))
    # keep<-random_res$"V2"[which(unlist(lapply(random_res$V1,is.null)))]
    # sel_idx<-c(which(tissue_exp$peak_count%in%keep),sel_idx)
    #tissue_exp_sel<-tissue_exp[sel_idx,]
    tissue_exp_sel<-tissue_exp
    random_genes_exp_pool<-list()
    housekeeping_genes_exp_pool<-list()
    file_save<-paste0(GTEx_v8_loc,"random_and_housekeeping_",iterations,"_",tissue,".RData")
    if(!file.exists(file_save)){
      for(i in 1:iterations){
        print(paste0("Progress in ",i/iterations*100,"%"))
        set.seed(i)
        random_genes<-scaled_tmp$Description[sample(nrow(scaled_tmp),nrow(tissue_exp_sel),replace = T)]
        random_genes_exp_pool[[i]]<-scaled_tmp[match(random_genes,scaled_tmp$Description),tissue]
        set.seed(i)
        housekeeping_genes<-sample(house_keeping_sets,nrow(tissue_exp_sel),replace = T)
        housekeeping_genes_exp_pool[[i]]<-scaled_tmp[match(housekeeping_genes,scaled_tmp$Description),tissue]
      }
      random_genes_exp_pool_df<-list.cbind(random_genes_exp_pool)
      housekeeping_genes_exp_pool_df<-list.cbind(housekeeping_genes_exp_pool)
      save(random_genes_exp_pool_df,housekeeping_genes_exp_pool_df,file = file_save)
      rm(random_genes_exp_pool,housekeeping_genes_exp_pool)
      gc()
    }else{
      load(file_save)
    }
    tissue_exp_sel$"random_genes_exp_mean"<-apply(random_genes_exp_pool_df,1,mean)
    tissue_exp_sel$"housekeeping_genes_exp_mean"<-apply(housekeeping_genes_exp_pool_df,1,mean)
    rm(random_genes_exp_pool_df,housekeeping_genes_exp_pool_df)
    gc()
    tissue_exp_sel_long<-pivot_longer(data = tissue_exp_sel,names_to = "groups",values_to = "exp",cols = c(all_of(paste0("scaled_",tissue)),random_genes_exp_mean,housekeeping_genes_exp_mean))
    tissue_exp_sel_long<-tissue_exp_sel_long[,c("peak_count","groups","exp")]
    tissue_exp_sel_long$groups<-factor(tissue_exp_sel_long$groups,levels = c("random_genes_exp_mean" ,"housekeeping_genes_exp_mean",paste0("scaled_",tissue)))
    p2<-ggplot(data = tissue_exp_sel_long,aes(x=peak_count,y=exp,fill=groups))+geom_boxplot(aes(factor(peak_count)),outlier.size = 0,outlier.color = NA)+geom_jitter(size=1,aes(x=peak_count,shape=groups),position = position_dodge(0.75))+theme_classic()
    ggsave(filename = paste0("/home/fangzj/Workdir/Basset_Unblocked/picture/GTEx_out/Exp_box_cor/with_random_housekeeping/",tissue,"_",judge_type,"_",iterations,"box_with_random_housekeeping.png"),p2,width = 15,height = 10)
    library("Rmisc")
    tgc<-summarySE(tissue_exp_sel_long,measurevar = "exp",groupvars = c("peak_count","groups"))
    p3<-ggplot(tgc,aes(x=peak_count,y=exp,colour=groups))+geom_errorbar(aes(ymin=exp-se,ymax=exp+se),width=0.1)+geom_line()+geom_point()+theme_classic()
    ggsave(filename = paste0("/home/fangzj/Workdir/Basset_Unblocked/picture/GTEx_out/Exp_box_cor/with_random_housekeeping/",tissue,"_",judge_type,"_",iterations,".pdf"),p3,width = 8,height = 5)
    }
}
for(j in sel_tissues_1){
  print(j)
  get_one_tissue_box(tissue = j,judge_type = "scaled_exp",iterations = 1000)
  gc()
}

#fold change
peak_foldchg_exp<-aggregate(mean_foldchg_exp~peak_count,mean,data=overlaps)
cor.test(as.numeric(peak_foldchg_exp$peak_count),peak_foldchg_exp$mean_foldchg_exp)
p2<-ggplot(data = peak_foldchg_exp,aes(x=peak_count,y=mean_foldchg_exp))+geom_point()+geom_smooth(method = "lm",formula = y~x,color="red",fill="pink")+stat_cor(digits = 4,method = "pearson")+theme_classic()
ggsave(filename = "/home/fangzj/Workdir/Basset_Unblocked/picture/GTEx_out/Exp_box_cor/peak_foldchg_exp.pdf",p2,width = 9,height = 6)


overlaps_for_box<-overlaps[overlaps$peak_count>=0,]
overlaps_for_box<-pivot_longer(data = overlaps_for_box, cols = matches("^scaled_(?!mean_exp)",perl = T), names_to = "name", values_to = "value")
#ggplot(data = overlaps_for_box,aes(x=factor(peak_count),y=value))+geom_boxplot(outlier.colour = NA)+geom_jitter(width = 0.1,size=0.001)+theme_classic()
p3<-ggplot(data = overlaps_for_box,aes(x=peak_count,y=value))+geom_jitter(width = 1,size=0.001)+geom_smooth(method = "lm",formula = y~x,color="red",fill="pink")+stat_cor(digits = 4,method = "pearson")+theme_classic()
ggsave(filename = "/home/fangzj/Workdir/Basset_Unblocked/picture/GTEx_out/Exp_box_cor/mean_exp_with_tissue_box_allpoint.png",p3,width = 9,height = 6)
#p3<-ggplot(data = overlaps_for_box,aes(x=factor(peak_count),y=value))+geom_boxplot(outlier.colour = NA)+geom_count(position =position_jitter(width=0.05))+scale_size(range = c(0.1,3))
cor.test(as.numeric(overlaps$peak_count),as.numeric(overlaps$mean_exp))

###查看每个组织特有的peaks no use
tissue_specifc_region_ls<-list()
for(i in colnames(act_file)){
  tissue_specifc_region_ls[[i]]<-rownames(act_file)[act_file[,i]==1]
}




##
library("Seurat")
library("ggpubr")
all_tissue_rds<-readRDS("/media/fangzj/Expansion/ATAC_GWAS_DB_bk/ATAC_GWAS_DB/HUSCH/Output/alltissues/rds/alltissues.rds")
sel_rds<-subset(all_tissue_rds,orig.ident%in%c("Breast","Lung","Kidney","Skin"))
rm(all_tissue_rds)
gc()
#
#set.seed(1)
#idx<-sample(ncol(sel_rds),2000)
#sel_rds<-sel_rds[,idx]

#VlnPlot(sel_rds,c("NEAT1","EIF1","TPM4"),group.by = "orig.ident")
#averages<-AverageExpression(sel_rds,features = c("NEAT1","EIF1","TPM4","NFKBIA"),group.by = "orig.ident")
uniqued_gene_symbol<-overlaps%>%dplyr::group_by(peak_count)%>%dplyr::reframe(unique(gene_symbol))
colnames(uniqued_gene_symbol)[2]<-"gene_symbol"
uniqued_gene_symbol_ls<-split(uniqued_gene_symbol,uniqued_gene_symbol$peak_count)
uniqued_gene_symbol_ls<-lapply(uniqued_gene_symbol_ls,function(x){
  x<-x[,2]
  x<-as.vector(unlist(x))
  })
save(uniqued_gene_symbol_ls,file = paste0(GTEx_v8_loc,"uniqued_gene_symbol_ls.RData"))
#house_keeping_gene_loc<-"/home/fangzj/Workdir/Basset_Unblocked/data/GeneSets/Housekeeping/HOUNKPE_HOUSEKEEPING_GENES.v2023.2.Hs.gmt"
#house_keeping_sets<-clusterProfiler::read.gmt(house_keeping_gene_loc)
#uniqued_gene_symbol_ls$"house_keeping"<-house_keeping_sets$gene


#GO enrich
#fivenum(as.numeric(overlaps$peak_count))[4]
load(paste0(GTEx_v8_loc,"uniqued_gene_symbol_ls.RData"))
GO_res_raw<-list()
GO_res_simp<-list()
lapply(seq(19,10,-1), function(x){
  print(x)
  high_freq_peak_gene<-c()
  for(i in as.character(which(as.numeric(names(uniqued_gene_symbol_ls))>=x))){
    high_freq_peak_gene<-c(high_freq_peak_gene,uniqued_gene_symbol_ls[[i]]) 
  }
  high_freq_peak_gene<-unique(high_freq_peak_gene)
  go_enrich<-clusterProfiler::enrichGO(gene = high_freq_peak_gene,
                                       ont = "all",
                                       keyType = "SYMBOL",
                                       OrgDb = org.Hs.eg.db,
                                       pAdjustMethod = "BH",
                                       pvalueCutoff = 0.05,
                                       qvalueCutoff = 0.05)
  if(nrow(as.data.frame(go_enrich))!=0){
    go_simp<-clusterProfiler::simplify(go_enrich,cutoff=0.7,by="p.adjust",select_fun=min)
    GO_res_simp[[paste0("gene_set_at_least",x)]]<<-go_simp
  }
  GO_res_raw[[paste0("gene_set_at_least",x)]]<<-go_enrich
})
save(GO_res_raw,GO_res_simp,file = "/home/fangzj/Workdir/Basset_Unblocked/data/GTEx/GTEx_V8/GO_res/GO_res_with10genesets.RData")

lapply(names(GO_res_simp),function(x){
  p4<-barplot(GO_res_simp[[x]],
              x="GeneRatio",
              color="p.adjust",
              showCategory=10,
              split="ONTOLOGY",
              label_format=Inf)+facet_grid(ONTOLOGY~.,space='free_y',scale='free_y')
  ggsave(filename = paste0("/home/fangzj/Workdir/Basset_Unblocked/picture/GTEx_out/GO_res/",x,".pdf"),p4,width = 15,height = 10)
})

sce_score<-AddModuleScore(sel_rds,features = uniqued_gene_symbol_ls,name = "features",ctrl = 1000)
save(sce_score,file="/home/fangzj/Workdir/Basset_Unblocked/data/GTEx/GTEx_V8/sce_score_with_ctrl1000.RData")
#VlnPlot(sce_score,features = c("features1","features19"),group.by = "orig.ident",pt.size = 0)
res_scores_long<-pivot_longer(data = sce_score@meta.data,names_to = "features",values_to = "gene_set_score",cols = starts_with("features"))
res_scores_long$features2<-as.numeric(mapvalues(res_scores_long$features,from = paste0("features",seq(1,19)),to = seq(1,19)))

ggviolin(data = res_scores_long,
         x="features2",y="gene_set_score",
         width = 0.8,color = "black",
         xlab = F,
         add = "boxplot")
#
#score_for_plot<-res_scores_long%>%dplyr::group_by(features2,orig.ident)%>%dplyr::summarise(mean_score=mean(gene_set_score))
#ggplot(data = score_for_plot,aes(x=features2,y=mean_score))+geom_boxplot(aes(factor(features2)),outlier.size = 0,outlier.color = NA)+geom_smooth(method = "lm",formula = y~x,color="red",fill="pink")+geom_jitter(width = 0.1,size=1)+stat_cor(digits = 4,method = "pearson",size=5)+theme_classic()+theme(text = element_text(size = 20))
p5<-ggplot(data = res_scores_long,aes(x=features2,y=gene_set_score))+geom_boxplot(aes(factor(features2)),outlier.size = 0,outlier.color = NA)+geom_smooth(method = "lm",formula = y~x,color="red",fill="pink")+geom_jitter(width = 0.1,size=0.01)+stat_cor(digits = 4,method = "pearson",size=5)+theme_classic()+theme(text = element_text(size = 20))
ggsave(filename = "/home/fangzj/Workdir/Basset_Unblocked/picture/GTEx_out/Exp_box_cor/add_module_score.png",p5,width = 9,height = 6)
#
res_scores_long_tissue<-aggregate(gene_set_score~features2+orig.ident,mean,data=res_scores_long)
res_scores_long_tissue$orig.ident<-factor(res_scores_long_tissue$orig.ident,levels = c("Kidney","Lung","Breast","Skin"))
colors<-c("black","#FFBE7A","purple","#82B0D2")
p6<-ggplot(data = res_scores_long_tissue,aes(x=features2,y=gene_set_score))+geom_boxplot(aes(factor(features2)),outlier.size = 0,outlier.color = NA)+geom_smooth(method = "lm",formula = y~x,color="red",fill="pink")+geom_jitter(width = 0.1,size=1,aes(color=orig.ident))+scale_color_manual(values=colors)+stat_cor(digits = 4,method = "pearson",size=5)+theme_classic()+theme(text = element_text(size = 20))
ggsave(filename = "/home/fangzj/Workdir/Basset_Unblocked/picture/GTEx_out/Exp_box_cor/add_module_score_tissue.pdf",p6,width = 9,height = 6)
#
res_scores_long_tissue_overall<-aggregate(gene_set_score~features2,mean,data=res_scores_long)
p7<-ggplot(data = res_scores_long_tissue_overall,aes(x=features2,y=gene_set_score))+geom_point()+geom_smooth(method = "lm",formula = y~x,color="red",fill="pink")+stat_cor(digits = 4,method = "pearson",size=5)+theme_classic()+theme(text = element_text(size = 20))
ggsave(filename = "/home/fangzj/Workdir/Basset_Unblocked/picture/GTEx_out/Exp_box_cor/add_module_score_overall.pdf",p7,width = 9,height = 6)


#GSVA
library(GSVA)
library(GSEABase)
library("limma")
GTEx_v8_loc<-"/home/fangzj/Workdir/Basset_Unblocked/data/GTEx/GTEx_V8/"
load(paste0(GTEx_v8_loc,"uniqued_gene_symbol_ls.RData"))
test_genesets<-uniqued_gene_symbol_ls[1:19]
mydata<-as.matrix(sel_rds@assays$RNA@counts)
set.seed(123)
mydata<-mydata[,sample(x = ncol(mydata),size = floor(0.25*ncol(mydata)))]
gc()
res <- gsva(expr=mydata, gset.idx.list=test_genesets,method="gsva",kcdf="Poisson",parallel.sz=30L)
save(res,file="/home/fangzj/Workdir/Basset_Unblocked/data/GTEx/GTEx_V8/gsva_res_ratio_0.25.RData")
load("/home/fangzj/Workdir/Basset_Unblocked/data/GTEx/GTEx_V8/gsva_res_ratio_0.25.RData")

res_t<-t(res)%>%as.data.frame()
res_mean<-as.data.frame(colMeans(res_t))
res_mean<-cbind(rownames(res_mean),res_mean)
colnames(res_mean)<-c("id","value")
res_mean$id<-as.numeric(res_mean$id)
p<-ggplot(res_mean,aes(x=id,y=value))+geom_point()+geom_smooth(method = "lm",formula = y~x,color="red",fill="pink")+stat_cor(digits = 4,method = "pearson")+theme_classic()
ggsave("/home/fangzj/Workdir/Basset_Unblocked/picture/GTEx_out/Exp_box_cor/gsva_mean_score.pdf",p,height = 6,width = 9)

res_t$"tissue"<-sel_rds@meta.data$orig.ident[match(rownames(res_t),rownames(sel_rds@meta.data))]
res_t_tissue<-res_t%>%group_by(tissue)%>%dplyr::summarise(across(everything(),mean))
res_t_tissue<-pivot_longer(data = res_t_tissue,names_to = "id",values_to = "mean_score",cols = !starts_with("tissue"))
res_t_tissue$id<-as.numeric(res_t_tissue$id)
res_t_tissue$tissue<-factor(res_t_tissue$tissue,levels = c("Kidney","Lung","Breast","Skin"))
colors<-c("black","#FFBE7A","purple","#82B0D2")
p<-ggplot(res_t_tissue,aes(x=id,y=mean_score))+geom_boxplot(aes(x=factor(id)),outlier.size = 0,outlier.color = NA)+geom_jitter(width = 0.1,size=1,aes(color=tissue))+scale_color_manual(values = colors)+geom_smooth(method = "lm",formula = y~x,color="red",fill="pink")+stat_cor(digits = 4,method = "pearson")+theme_classic()
ggsave("/home/fangzj/Workdir/Basset_Unblocked/picture/GTEx_out/Exp_box_cor/gsva_mean_score_with_tissue.pdf",p,height = 6,width = 9)




##差异分析#####no use####
Diff_se_obj<-list()
for(j in as.character(levels(factor(sel_rds@meta.data$orig.ident)))){ 
  type<-dplyr::select(sel_rds@meta.data,"orig.ident")
  type$orig.ident<-as.character(type$orig.ident)
  type$orig.ident[which(type$orig.ident!=j)]<-"control"
  type$orig.ident[which(type$orig.ident==j)]<-"case"
  grouP<-as.factor(type$orig.ident)
  desigN <- model.matrix(~ grouP + 0)
  rownames(desigN)<-colnames(res)
  comparE <- makeContrasts(grouPcase-grouPcontrol,levels=desigN)
  fiT <- lmFit(res, desigN)
  fiT2 <- contrasts.fit(fiT, comparE)
  fiT3 <- eBayes(fiT2)
  print(j)
  Diff_se_obj[[j]]<-topTable(fiT3,p.value=0.05,num=Inf,adjust.method ="BH")
  type<-dplyr::select(sel_rds@meta.data,"orig.ident")
  type$orig.ident<-as.character(type$orig.ident)
}

test<-lapply(Diff_se_obj,FUN = function(x){as.data.frame(cbind(rownames(x),x$t))})
n_row<-10
down<-as.data.frame(matrix(NA,nrow =n_row,ncol = length(unique(sel_rds@meta.data$orig.ident)),dimnames = list(names(test_genesets),as.character(levels(factor(sel_rds@meta.data$orig.ident))))))
#
for (k in names(test)){
  down[match(test[[k]]$V1,rownames(down)),k]<-as.numeric(unlist(test[[k]]$V2))
}
down[which(is.na(down),arr.ind = T)]<-0
pheatmap::pheatmap(down)
res_t<-t(res)
#####no use####






##enhancer region
get_enhancer_region<-function(cell_line){
  enhancer_loc<-paste0("/home/fangzj/Workdir/Basset_Unblocked/data/Enhancer_region/",cell_line,"/")
  cell_line_enhancer<-read.table(dir(enhancer_loc,full.names = T)[grep("hg38",dir(enhancer_loc))])
  cell_line_enhancer<-cell_line_enhancer[,1:3]
  colnames(cell_line_enhancer)<-c("chr","start","end")
  cell_line_enhancer<-cell_line_enhancer[cell_line_enhancer$chr%in%c(paste0("chr",seq(1,22)),"chrX"),]
  cell_line_enhancer<-makeGRangesFromDataFrame(cell_line_enhancer)
  wait_mid<-mid(cell_line_enhancer)
  wait_start<-wait_mid-300
  wait_end<-wait_mid+300
  region<-as.data.frame(cbind(as.character(seqnames(cell_line_enhancer)),wait_start,wait_end))
  colnames(region)<-c("chr","start","end")
  save_name<-paste0(enhancer_loc,cell_line,"_enhancer_600ps.bed")
  write.table(region,save_name,quote = F,sep="\t",col.names = F,row.names = F)
  setwd(enhancer_loc)
  hg38_fa<-"/home/fangzj/Workdir/Basset_Unblocked/data/genomes/hg38/chroms/hg38.fa"
  bedtools <- '/home/fangzj/Workdir/software/anaconda/envs/basset/bin/bedtools'
  get_seqs<-paste0(bedtools,' getfasta -fi ',hg38_fa, " -bed ",save_name, " -fo ", paste0(save_name,'.fa'))
  system(get_seqs)
  gc()
}
#"HT1080","K562","U2OS","Hela",
for(i in c("NHDF")){
  get_enhancer_region(i)
}


best_model_loc = '/home/fangzj/Workdir/Basset_Unblocked/data/cleandata/flank_seqs/all_classes_flank/both_ext/up_split/1e-90/200_600_85/quick_check/m1d_optrmsprop_drop0.39_lr0.00002_step14_gama0.28_mom0.97_batch512_ks1_12_conv1_300_norm_1e9,0_run_3/GWAS_locus/'
allfiles<-dir(best_model_loc,full.names = T)###严重警告！！在for循环中，dir读取不全？？？
cell_line_res<-read.table(allfiles[grep("HT1080.*\\.bed\\.fa\\.csv",allfiles,ignore.case = T)])
res<-as.data.frame(colMeans(cell_line_res))
res<-cbind(rownames(res),res)
colnames(res)<-c("cell_line","ratio")
median(res$ratio)
idx<-grep("U2OS",res$cell_line,ignore.case = T)
t.test(res[idx,"ratio"],res[-idx,"ratio"])



###
#找到在所有细胞系中高表达的基因（非差异表达)
library(clusterProfiler)
tissue_match_tpm<-scaled_tmp[,sel_tissues_1]
gene_rank<-apply(tissue_match_tpm,2,function(x){rank(dplyr::desc(x))})
save_loc<-"/home/fangzj/Workdir/Basset_Unblocked/data/GTEx/GTEx_V8/compare_high_exp_genes/tss/"
###trash
# get_mid_region<-function(grange_region,save_loc,save_name){
#   wait_mid<-mid(grange_region)
#   wait_start<-wait_mid-300
#   wait_end<-wait_mid+300
#   region<-as.data.frame(cbind(as.character(seqnames(grange_region)),wait_start,wait_end))
#   colnames(region)<-c("chr","start","end")
#   save_name<-paste0(save_loc,save_name)
#   write.table(region,save_name,quote = F,sep="\t",col.names = F,row.names = F)
#   setwd(save_loc)
#   hg38_fa<-"/home/fangzj/Workdir/Basset_Unblocked/data/genomes/hg38/chroms/hg38.fa"
#   bedtools <- '/home/fangzj/Workdir/software/anaconda/envs/basset/bin/bedtools'
#   get_seqs<-paste0(bedtools,' getfasta -fi ',hg38_fa, " -bed ",save_name, " -fo ", paste0(save_name,'.fa'))
#   system(get_seqs)
# }
###trash

get_grange_region<-function(grange_region,save_loc,save_name){
  wait_start<-start(grange_region)
  wait_end<-end(grange_region)
  region<-as.data.frame(cbind(as.character(seqnames(grange_region)),wait_start,wait_end))
  colnames(region)<-c("chr","start","end")
  save_name<-paste0(save_loc,save_name)
  write.table(region,save_name,quote = F,sep="\t",col.names = F,row.names = F)
  setwd(save_loc)
  hg38_fa<-"/home/fangzj/Workdir/Basset_Unblocked/data/genomes/hg38/chroms/hg38.fa"
  bedtools <- '/home/fangzj/Workdir/software/anaconda/envs/basset/bin/bedtools'
  get_seqs<-paste0(bedtools,' getfasta -fi ',hg38_fa, " -bed ",save_name, " -fo ", paste0(save_name,'.fa'))
  system(get_seqs)
}

library("GenomicFeatures")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <-TxDb.Hsapiens.UCSC.hg38.knownGene
seqnames_sel <- c(paste0("chr",seq(1,22)),"chrX")
seqlevels(txdb) <- seqnames_sel
all_transcripts<-transcripts(txdb)
tss_region<-promoters(all_transcripts,upstream = 300,downstream = 301)
get_grange_region(grange_region = tss_region,save_loc = paste0(save_loc,'gene_region/'),save_name="all_genic_regions_tss_600bps.bed")
#
extract_grange<-function(ids,txdb,pool_region){
  tmp_gene_id<-str_replace(ids,pattern = "\\.\\d+","")
  mapped_gene_id<-bitr(tmp_gene_id,fromType ='ENSEMBL' ,toType ="ENTREZID",OrgDb = 'org.Hs.eg.db')
  mapped_meta<-AnnotationDbi::select(txdb, keys = mapped_gene_id$ENTREZID, keytype="GENEID", columns="TXNAME")
  mapped_meta<-mapped_meta[!is.na(mapped_meta$TXNAME),]
  mapped_tx_name<-mapped_meta$"TXNAME"
  TSS_sel_region<-pool_region[pool_region$tx_name%in%mapped_tx_name,]
  rank_class_Exp_genic_grange_df<-as.data.frame(TSS_sel_region)
  rank_class_Exp_genic_grange_df$"gene_id"<-mapvalues(x = rank_class_Exp_genic_grange_df$tx_name,from = mapped_meta$TXNAME,to = mapped_meta$GENEID)
  gene_symbol_meta<-bitr(rank_class_Exp_genic_grange_df$gene_id,fromType ='ENTREZID' ,toType ="SYMBOL",OrgDb = 'org.Hs.eg.db')
  rank_class_Exp_genic_grange_df$"genesymbol"<-mapvalues(rank_class_Exp_genic_grange_df$gene_id,from = gene_symbol_meta$ENTREZID,to=gene_symbol_meta$SYMBOL)
  return(list(rank_class_Exp_genic_grange_df,TSS_sel_region))
}
get_special_region<-function(rank_class){
  #genic_grange_ls<-list()
  #seq(0.01,0.2,0.01)
  if(rank_class=="high"){range_<-seq(0.01,0.2,0.01)}else{range_<-0.01}
  for(i in range_){
    print(i)
    if(rank_class=="high"){
      all_ranks_nu<-apply(gene_rank,2,max)
    gene_rank_res<-apply(gene_rank,1,function(x){
      all(x<=i*all_ranks_nu)
    }) %>%as.matrix()}else{
      if(rank_class=="low"){
        gene_rank_res<-apply(gene_rank,1,function(x){
          all(x>=((1-i)*all_ranks_nu))
        }) %>%as.matrix()
      }
      }
    gene_rank_res<-cbind(rownames(gene_rank_res),gene_rank_res)
    colnames(gene_rank_res)<-c("gene_id","judge")
    gene_rank_res<-as.data.frame(gene_rank_res)
    gene_rank_res$"gene_symbol"<-mapvalues(gene_rank_res$gene_id,tmp$Name,tmp$Description)
    rank_class_Exp_gene<-gene_rank_res[gene_rank_res$judge=="TRUE",]
    ids = rank_class_Exp_gene$gene_id
    res<-extract_grange(ids = ids,txdb = txdb,pool_region = tss_region)
    rank_class_Exp_genic_grange_df<-res[[1]]
    TSS_sel_region<-res[[2]]
    #genic_grange_ls[[paste0('exp_ratio_',i)]]<-rank_class_Exp_genic_grange_df
    write.table(rank_class_Exp_genic_grange_df,paste0(save_loc,'gene_list/',rank_class,'_exp_ratio_',i,".txt"),sep = "\t",quote = F,row.names = F,col.names = T)
    get_grange_region(grange_region = TSS_sel_region,save_loc = paste0(save_loc,'gene_region/'),save_name=paste0(rank_class,"_exp_600ps_exp_ratio_",i,".bed"))
  }
}
get_special_region(rank_class = 'high')
get_special_region(rank_class = 'low')
###
best_model_loc = '/home/fangzj/Workdir/Basset_Unblocked/data/cleandata/flank_seqs/all_classes_flank/both_ext/up_split/1e-90/200_600_85/quick_check/m1d_optrmsprop_drop0.39_lr0.00002_step14_gama0.28_mom0.97_batch512_ks1_12_conv1_300_norm_1e9,0_run_3/GWAS_locus/'
allfiles<-dir(best_model_loc,full.names = T)###严重警告！！在for循环中，dir读取不全？？？

load("/home/fangzj/Workdir/Basset_Unblocked/data/GTEx/GTEx_V8/meta_info_match.RData")
genic_res1<-read.table(allfiles[grep("all_genic_regions_tss.*\\.bed\\.fa\\.csv",allfiles,ignore.case = T)])
res1<-as.data.frame(colMeans(genic_res1))
res1<-cbind(rownames(res1),res1)
colnames(res1)<-c("cell_line","ratio")
res1<-res1[res1$cell_line%in%meta_info$cell_line,]
res1$"tissue"<-mapvalues(res1$cell_line,from = meta_info$cell_line,to = meta_info$matched_tissue)
res1<-aggregate(ratio~tissue,mean,data=res1)
res1$"exp_ratio_nu"<-"all_genic"
res1<-dplyr::select(res1,exp_ratio_nu,tissue,ratio)
res1$'group'<-'all_genic'

pattern1<-"high_exp_600ps_exp_ratio_"
pattern2<-"high_exp_Kidney_ratio_"
pattern3<-"high_exp_Lung_ratio_"
pattern4<-"high_exp_Breast_ratio_"
pattern5<-"high_exp_Skin_ratio_"
genic_res2_ls<-list()

for(i in seq(0.01,0.2,0.01)){
  for(pattern_str in c(pattern1,pattern2,pattern3,pattern4,pattern5)){
  #print(i)
  pattern<-paste0(best_model_loc,pattern_str,i,".bed.fa.csv")
  genic_res2<-read.table(pattern)
  res2<-as.data.frame(colMeans(genic_res2))
  res2<-cbind(rownames(res2),res2)
  colnames(res2)<-c("cell_line","ratio")
  res2<-res2[res2$cell_line%in%meta_info$cell_line,]
  res2$"exp_ratio"<-paste0(pattern_str,i)
  genic_res2_ls[[paste0(pattern_str,"ratio_",i)]]<-res2
}
}

genic_res2_df<-list.rbind(genic_res2_ls)
genic_res2_df$"tissue"<-mapvalues(genic_res2_df$cell_line,from = meta_info$cell_line,to = meta_info$matched_tissue)
genic_res2_df$"exp_ratio_nu"<-as.numeric(str_replace(genic_res2_df$exp_ratio,pattern = ".*_",replacement = ""))
genic_res2_df$'group'<-str_replace(genic_res2_df$exp_ratio,pattern = "_ratio.*",replacement = "")
genic_res2_df<-aggregate(ratio~exp_ratio_nu+tissue+group,mean,data=genic_res2_df)
pic_save_loc<-"/home/fangzj/Workdir/Basset_Unblocked/picture/GTEx_out/Exp_box_cor/model_predict/"

genic_res3<-read.table(allfiles[grep("low_exp_600ps_exp_ratio_0.01.*\\.bed\\.fa\\.csv",allfiles,ignore.case = T)])
res3<-as.data.frame(colMeans(genic_res3))
res3<-cbind(rownames(res3),res3)
colnames(res3)<-c("cell_line","ratio")
res3<-res3[res3$cell_line%in%meta_info$cell_line,]
res3$"tissue"<-mapvalues(res3$cell_line,from = meta_info$cell_line,to = meta_info$matched_tissue)
res3<-aggregate(ratio~tissue,mean,data=res3)
res3$"exp_ratio_nu"<-"lowest_exp"
res3<-dplyr::select(res3,exp_ratio_nu,tissue,ratio)
res3$'group'<-'lowest_exp'

set.seed(1)
group_colors<-c("black","#FFBE7A","purple","#82B0D2",'#F0A19A','#00A664','blue')
names(group_colors)<-c("high_exp_Kidney","high_exp_Lung","high_exp_Breast","high_exp_Skin","high_exp_600ps_exp",'all_genic',"lowest_exp")
colors<-c("black","#FFBE7A","purple","#82B0D2")
names(colors)<-c('Kidney','Lung','Breast','Skin')
df_with_all_genic<-rbind(genic_res2_df,res1,res3)
df_with_all_genic$group<-factor(df_with_all_genic$group,levels = c("high_exp_600ps_exp","high_exp_Kidney","high_exp_Skin","high_exp_Breast","high_exp_Lung",'all_genic',"lowest_exp"))
df_with_all_genic$tissue<-factor(df_with_all_genic$tissue,levels = c("Kidney","Lung","Breast","Skin"))
df_with_all_genic$colour<-mapvalues(df_with_all_genic$tissue,from = names(colors),colors)
df_with_all_genic$'judge'<-paste0(paste(df_with_all_genic$tissue,df_with_all_genic$group,sep = "_"))
df_with_all_genic_sel<-df_with_all_genic[df_with_all_genic$judge%in%c('Kidney_high_exp_Kidney','Skin_high_exp_Skin','Breast_high_exp_Breast','Lung_high_exp_Lung')|df_with_all_genic$group%in%c('high_exp_600ps_exp','all_genic','lowest_exp'),]
df_with_all_genic_sel$exp_ratio_nu[df_with_all_genic_sel$exp_ratio_nu=='all_genic']=0.21
df_with_all_genic_sel$exp_ratio_nu[df_with_all_genic_sel$exp_ratio_nu=='lowest_exp']=0.22
df_with_all_genic_sel$exp_ratio_nu<-as.numeric(df_with_all_genic_sel$exp_ratio_nu)
df_with_all_genic_sel$'plot_group'<-str_replace(string =df_with_all_genic_sel$group ,pattern = 	"high_exp_(Kidney|Lung|Breast|Skin)",replacement = "high_exp_tissue")
#df_with_all_genic_sel<-df_with_all_genic[df_with_all_genic$group%in%c('high_exp_600ps_exp'),]

df_with_all_genic_sel$group<-as.factor(df_with_all_genic_sel$group)
p1<-ggplot(df_with_all_genic,aes(x=factor(exp_ratio_nu),y=ratio,fill=group))+geom_boxplot(outlier.size = 0,outlier.color = NA,position = position_dodge(),alpha=0.6)+geom_jitter(position = position_jitterdodge(jitter.width = 0.05,seed = 1),size=1,colour=df_with_all_genic$colour)+theme_classic()+scale_fill_manual(values=group_colors)

p1_1<-ggplot(data = df_with_all_genic_sel,aes(x=as.numeric(exp_ratio_nu),y=ratio,shape=plot_group,color=tissue))+geom_line(linewidth=1)+geom_point(size=3)+theme_classic()+scale_color_manual(values = colors)


ggsave(paste0(pic_save_loc,"with_all_genic_boxplot.pdf"),p1,width = 15,height = 10)
ggsave(paste0(pic_save_loc,"with_all_genic_linechat.pdf"),p1_1,width = 15,height = 10)
set.seed(1)

genic_res2_df$tissue<-factor(genic_res2_df$tissue,levels = c("Kidney","Lung","Breast","Skin"))
genic_res2_df$colour<-mapvalues(genic_res2_df$tissue,from = names(colors),colors)
genic_res2_df$group<-factor(genic_res2_df$group,levels = c("high_exp_600ps_exp","high_exp_Kidney","high_exp_Skin","high_exp_Breast","high_exp_Lung"))
p2<-ggplot(genic_res2_df,aes(x=factor(exp_ratio_nu),y=ratio,fill=factor(group)))+geom_boxplot(outlier.size = 0,outlier.color = NA,alpha=0.6)+geom_point(position = position_jitterdodge(),size=1,color=genic_res2_df$colour)+theme_classic()+scale_fill_manual(values=group_colors)
ggsave(paste0(pic_save_loc,"no_genic_boxplot.pdf"),p2,width = 15,height = 10)
set.seed(1)
#p3<-ggplot(genic_res2_df,aes(x=exp_ratio_nu,y=ratio))+geom_jitter(width = 0.00001,size=1)+geom_smooth(method = "lm",formula = y~x,color="red",fill="pink")+stat_cor(digits = 4,method = "pearson")+theme_classic()
#ggsave(paste0(pic_save_loc,"cor_test.pdf"),p3,width = 9,height = 6)


#ratio_res<-aggregate(ratio~exp_ratio,data=genic_res2_df,mean)
#ggplot(data = ratio_res,aes(x=exp_ratio_nu,y=ratio))+geom_point()+geom_smooth(method = "lm",formula = y~x,color="red",fill="pink")+stat_cor(digits = 4,method = "pearson",size=5)


##GO enrichment
library(org.Hs.eg.db)
genelist_001<-read.table('/home/fangzj/Workdir/Basset_Unblocked/data/GTEx/GTEx_V8/compare_high_exp_genes/tss/gene_list/high_exp_ratio_0.01.txt',sep = "\t",header = T)
high_freq_peak_gene<-unique(genelist_001$genesymbol)
go_enrich<-clusterProfiler::enrichGO(gene = high_freq_peak_gene,
                                     ont = "all",
                                     keyType = "SYMBOL",
                                     OrgDb = org.Hs.eg.db,
                                     pAdjustMethod = "BH",
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.05)
go_simp<-clusterProfiler::simplify(go_enrich,cutoff=0.7,by="p.adjust",select_fun=min)

p<-barplot(go_simp,
            x="GeneRatio",
            color="p.adjust",
            showCategory=10,
            split="ONTOLOGY",
            label_format=Inf)+facet_grid(ONTOLOGY~.,space='free_y',scale='free_y')
ggsave(filename = "/home/fangzj/Workdir/Basset_Unblocked/data/GTEx/GTEx_V8/compare_high_exp_genes/GO/peak_001.pdf",plot = p,width = 8,height = 6)


#看每个组织高表达的基因（非差异表达基因）
tissue_match_tpm<-scaled_tmp[,sel_tissues_1]
gene_rank<-apply(tissue_match_tpm,2,function(x){rank(dplyr::desc(x))})%>%as.data.frame()
genelist_per_tissue<-list()
for(tissue in sel_tissues_1){
  gene_exp_ranked<-dplyr::arrange(gene_rank,gene_rank[,tissue])
  max_rank<-max(gene_exp_ranked[,tissue])
 for (ratio in seq(0.01,0.2,0.01)) {
   genesets<-rownames(gene_exp_ranked)[gene_exp_ranked[,tissue]<=ratio*max_rank]
   genesets<-Name_Description[Name_Description$Name%in%genesets,]
   genelist_per_tissue[[tissue]][[paste0("ratio_",ratio)]]<-genesets
 }
  }

skin_1<-genelist_per_tissue$Skin...Not.Sun.Exposed..Suprapubic.
skin_2<-genelist_per_tissue$Skin...Sun.Exposed..Lower.leg.
skin_union<-list()
for (ratio in seq(0.01,0.2,0.01)) {
  test1<-skin_1[[paste0("ratio_",ratio)]]
  test2<-skin_2[[paste0("ratio_",ratio)]]
  unions<-rbind(test1,test2)
  unions<-unique(unions)
  skin_union[[paste0("ratio_",ratio)]]<-unions
}
kidney_1<-genelist_per_tissue$Kidney...Cortex
kidney_2<-genelist_per_tissue$Kidney...Medulla
kidney_union<-list()
for (ratio in seq(0.01,0.2,0.01)) {
  test1<-kidney_1[[paste0("ratio_",ratio)]]
  test2<-kidney_2[[paste0("ratio_",ratio)]]
  unions<-rbind(test1,test2)
  unions<-unique(unions)
  kidney_union[[paste0("ratio_",ratio)]]<-unions
}

genelist_per_tissue_sel<-list()
genelist_per_tissue_sel[['Kidney']]<-kidney_union
genelist_per_tissue_sel[['Lung']]<-genelist_per_tissue$Lung
genelist_per_tissue_sel[['Breast']]<-genelist_per_tissue$Breast...Mammary.Tissue
genelist_per_tissue_sel[['Skin']]<-skin_union
#
library("GenomicFeatures")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <-TxDb.Hsapiens.UCSC.hg38.knownGene
seqnames_sel <- c(paste0("chr",seq(1,22)),"chrX")
seqlevels(txdb) <- seqnames_sel
all_transcripts<-transcripts(txdb)
tss_region<-promoters(all_transcripts,upstream = 300,downstream = 301)
#
save_loc<-'/home/fangzj/Workdir/Basset_Unblocked/data/GTEx/GTEx_V8/compare_high_exp_genes/tss_across_tissues/'
lapply(names(genelist_per_tissue_sel),function(x){
  tmp<-genelist_per_tissue_sel[[x]]
  lapply(names(tmp), function(y){
    tmp2<-tmp[[y]]
    gene_ids<-tmp2$"Name"
    res<-extract_grange(ids = gene_ids,txdb = txdb,pool_region = tss_region)
    rank_class_Exp_genic_grange_df<-res[[1]]
    TSS_sel_region<-res[[2]]
    write.table(rank_class_Exp_genic_grange_df,paste0(save_loc,'/gene_list/',x,'/high_exp_',x,'_',y,".txt"),sep = "\t",quote = F,row.names = F,col.names = T)
    get_grange_region(grange_region = TSS_sel_region,save_loc = paste0(save_loc,'/gene_region/',x,'/'),save_name=paste0('high_exp_',x,'_',y,".bed"))
  })
})
