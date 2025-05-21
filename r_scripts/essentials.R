library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)
library(tidyr)
library(stringr)
library(rlang)
library(RColorBrewer)
library(varhandle)
library(gtools)
library(R.utils)
library(hacksaw)
library(scales)
require(patchwork)
library(gridExtra)
library(qpdf)
library(taxonomizr)
library(MASS)
library(lsa)
library(readODS)
library(readxl)
library(pheatmap)
library(vegan)
library(nlme)
library(ggrepel)
library(ggpubr)
library(stringr)
library(lme4)
library(ggvenn)

source("functions/distributions.R")
library(fitdistrplus)


sdf<-read.table("helpfiles/summary_samples.tsv", sep = "\t", header = F) #sample treads breads vreads
cdf<-read.table("helpfiles/summary_contigs.tsv", sep = "\t", header = F)
#cdf<-subset(cdf, V3>=2000)
colnames(sdf)<-c("sample","total")
colnames(cdf)<-c("hepID","ID","length")
heps<-unique(cdf$hepID)
samples<-unique(sdf$sample)
### metadata
metadf<-read.table("helpfiles/metadata.tsv", sep = "\t", stringsAsFactors = F, header = T, fill = T, quote = '', na.strings = "")
metadf$sa<-sapply(metadf$Sampling.age,function(x){if(x=="12kk"){return("12_months")};if(x=="6kk"){return("6_months")};if(x=="1kk"){return("1_month")};if(x=="7kk"){return("7_months")};if(x=="2kk"){return("2_months")}})
metadf$sample<-paste0("HeP-",metadf$Perhe_ID,"_",metadf$sa)
meta_bysample_df<-read.table("helpfiles/metadata_bysample.tsv", sep = "\t", stringsAsFactors = F, header = T, fill = T, quote = '"', na.strings = "NA",numerals = "no.loss")
#meta_bysample_df<-read_ods("helpfiles/HELMI_plussa_samplelist_HK20230705.ods",strings_as_factors = F)
metamaster<-read_ods("helpfiles/HELMI_plussa_Master_list_010724.ods")
metamaster$sample<-metamaster$Sample_ID
meta_trust<-read.table("helpfiles/trusted_table.tsv", sep = "\t", header = T)
blacklist<-read.table("helpfiles/blacklist.txt", header = F)$V1
supptable<-read.table("helpfiles/supp_table1.tsv", header = F, sep = "\t")
colnames(supptable)<-c("hepID","sample","day","sampling_age","phage_enriched","delivery","X16S")

add_metadata<-function(df){
  if(!"samples" %in% colnames(df) && !"sample" %in% colnames(df)){
    df$delivery<-sapply(df$hepID, function(x){opt<-metadf[which(x==metadf$sample),"Delivery.mode"]; unique(opt[!is.na(opt)])})
    df$phage_enriched<-sapply(df$hepID, function(x){metadf[which(x==metadf$sample),"Phage.metagenomics"][1]})
    if(grepl("month",df$hepID[1])){
      df$age<-sapply(df$hepID, function(x){tmp<-strsplit(x,"_")[[1]];paste0(tmp[2]," ",tmp[3])})
    }
    return(df)
  }
  if(class(df$sample)=="factor"){
    df$sample<-unfactor(df$sample)  
  }
  if(!"hepID" %in% colnames(df)){
    df$hepID<-sapply(df$sample, function(x){tmp<-strsplit(x,"-")[[1]];paste0(tmp[1],"-",tmp[2])})
  }
  df$sample<-factor(df$sample, levels = unique(df[mixedorder(df$sample,decreasing = T),]$sample))
  df$sampling_age<-sapply(df$sample, function(s){meta_bysample_df[which(meta_bysample_df$sample==s),]$Sample_type})
  df$sa<-sapply(df$sampling_age,function(x)gsub(pattern = " ","_",x))
  df$hepID<-paste0(df$hepID,"_",df$sa)
  df$sa<-NULL
  df$delivery<-sapply(df$hepID, function(x){opt<-metadf[which(x==metadf$sample),"Delivery.mode"]; unique(opt[!is.na(opt)])})
  df<-subset(df, delivery %in% c("Vaginal","Cesarean"))
  df$phage_enriched<-sapply(df$hepID, function(x){metadf[which(x==metadf$sample),"Phage.metagenomics"][1]})
  df$day<-sapply(df$sample, function(x){as.numeric(strsplit(as.character(x),"-")[[1]][3])})
  mins<-df %>% group_by(hepID) %>% summarise(day=min(day))
  df$representative<-sapply(df$sample, function(x){
    tmp<-strsplit(as.character(x),"-")[[1]]
    hep<-paste(tmp[1:2], collapse = "-")
    d<-as.numeric(tmp[3])
    nrow(subset(mins, hepID==hep & day==d))==1
  })
  df
}

get_legend2 <- function(plot, legend = NULL) {
  
  gt <- ggplotGrob(plot)
  
  pattern <- "guide-box"
  if (!is.null(legend)) {
    pattern <- paste0(pattern, "-", legend)
  }
  
  indices <- grep(pattern, gt$layout$name)
  
  not_empty <- !vapply(
    gt$grobs[indices], 
    inherits, what = "zeroGrob", 
    FUN.VALUE = logical(1)
  )
  indices <- indices[not_empty]
  
  if (length(indices) > 0) {
    return(gt$grobs[[indices[1]]])
  }
  return(NULL)
}

print_stats<-function(df){
  print("########################")
  print(paste0("total: ",nrow(df)))
  ca<-round(nrow(subset(df, !is.na(family)))/nrow(df),2)
  print(paste0(ca," of contigs are family tax annotated | ",nrow(subset(df, !is.na(family)))))
  ca<-round(nrow(subset(df, !is.na(family) & !grepl(";",family)))/nrow(df),2)
  print(paste0(ca," of contigs are pure tax annotated | ",nrow(subset(df, !is.na(family) & !grepl(";",family)))))
  ca<-round(nrow(subset(df, !is.na(genus)))/nrow(df),2)
  print(paste0(ca," of contigs are genus tax annotated | ", nrow(subset(df, !is.na(genus)))))
  ca<-round(nrow(subset(df, !is.na(isphage)))/nrow(df),2)
  print(paste0(ca," of isphage annotated | ",nrow(subset(df, !is.na(isphage)))))
  ca<-round(nrow(subset(df, !is.na(host_g)))/nrow(df),2)
  print(paste0(ca," of host_g annotated | ",nrow(subset(df, !is.na(host_g)))))
  ca<-round(nrow(subset(df, !is.na(host_s)))/nrow(df),2)
  print(paste0(ca," of host_s annotated | ",nrow(subset(df, !is.na(host_s)))))
  ca<-round(nrow(subset(df, !is.na(phagetype)))/nrow(df),2)
  print(paste0(ca," of phagetype annotated | ",nrow(subset(df, !is.na(phagetype)))))
  print("########################")
  print(paste0("integral: ",round(nrow(subset(df, !is.na(family) | !is.na(host_g) | !is.na(phagetype)))/nrow(df),4)))
}

p2t<-function(p){
  cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1)
  symbols = c("****", "****", "***", "**", "*","")
  s<-symbols[min(which(p<=cutpoints))]
  if(is.na(s)){"ns"}
  s
}

library(CHNOSZ)
hgdbnames<-getnames("/data/tmp/HumGutDB/")
hgdbnodes<-getnodes("/data/tmp/HumGutDB/")
colnames(hgdbnodes)<-c("tid","parenttid","rank")
hgdbnames<-hgdbnames %>% filter(type=="scientific name")
colnames(hgdbnames)[1]<-"tid"
ncbi<-merge(hgdbnodes,hgdbnames[,c(1,2)], by="tid", all.x=T) %>% unique()
climbranks<-function(tid,lvl){
  idx<-which(ncbi$tid==tid)
  if(length(idx)==0){
    return(NA)
  }
  if(ncbi[idx,"tid"]==1){
    return(NA)
  }
  if(ncbi[idx,"rank"]==lvl){
    return(ncbi[idx,"name"])
  }else{
    return(climbranks(ncbi[idx,"parenttid"],lvl))
  }
}
name2rank<-function(n,rank="genus"){
  n<-gsub("'","",n)
  idx<-which(ncbi$name==n)
  climbranks(ncbi[idx,"tid"],rank)
}


pcolors<-c("Actinobacteria"="#1B9E77","Bacteroidetes"="#D95F02","Firmicutes"="#7570B3","Fusobacteria"="#E7298A","Proteobacteria"="#66A61E","Verrucomicrobia"="#E6AB02")
