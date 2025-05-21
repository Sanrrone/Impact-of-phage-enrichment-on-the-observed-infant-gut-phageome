setwd("~/Desktop/phage_enrichment_newassembly")
#prophage = temperate | use lysogenic cycle
#phage = virulent | use lytic stage

#lytic stage = replicate along with host genome
#lysogenic stage = replicate and assemble independent (using host resources)

l<-list.files(path = "taxprof_propagate", pattern = "_propagate.tsv", full.names = T)
df_prop<-lapply(l, function(tsv){
  sname<-gsub("_propagate.tsv","",basename(tsv))
  tmp<-read.table(tsv, header = T, sep = "\t")
  colnames(tmp)[1]<-"ID"
  colnames(tmp)[3]<-"stage"
  tmp$sample<-sname
  tmp$contig_len_noprophage<-tmp$host_len
  tmp$phage_host_lenratio<-tmp$prophage_len/tmp$host_len
  tmp[tmp$contig_len_noprophage<2000,"stage"]<-"ambiguous"
  tmp[tmp$phage_host_lenratio>0.5,"stage"]<-"ambiguous"
  
  #tmp<-tmp %>% filter(!stage%in%c("not present","ambiguous"))
  
  #tmp$fluct<-tmp$prophage_sd_cov/ifelse(tmp$prophage_mean_cov,tmp$prophage_mean_cov,1)
  #tmp[tmp$fluct>0.5,"stage"]<-"ambiguous"
  #tmp<-tmp %>% filter(!stage%in%c("not present","ambiguous"))
  
  #tmp[!is.na(tmp$prophage_cov_breadth)&tmp$prophage_cov_breadth<0.9,"stage"]<-"ambiguous"
  tmp<-tmp %>% filter(!stage%in%c("not present","ambiguous"))
  #tmp[tmp$host_mean_cov<5 & tmp$prophage_mean_cov<5,"stage"]<-"dormant"
  #tmp[tmp$prophage.host_ratio<2, "stage"]<-"dormant"
  
  tmp[,c("sample","ID","stage","contig_len_noprophage")]
  
}) %>% bind_rows()

