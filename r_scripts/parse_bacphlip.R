setwd("~/Desktop/phage_enrichment_newassembly")
#prophage = temperate | use lysogenic cycle
#phage = virulent | use lytic stage

#lytic stage = replicate along with host genome
#lysogenic stage = replicate and assemble independent (using host resources)


l<-list.files(path = "taxprof_contig_bacphlip/", pattern = ".tsv", full.names = T)
df_bacph<-lapply(l, function(tsv){
  sname<-gsub("_bacphlip.tsv","",basename(tsv))
  tmp<-read.table(tsv, header = T, sep = "\t")
  colnames(tmp)[1]<-"ID"
  tmp<-subset(tmp, Virulent>0.8 | Temperate>0.8)
  tmp$phagetype<-ifelse(tmp$Virulent>tmp$Temperate,"virulent","temperate")
  tmp$lf_prob_bacphlip<-ifelse(tmp$Virulent>tmp$Temperate,tmp$Virulent,tmp$Temperate)
  tmp$sepID<-sname
  tmp[,c("ID","sepID","phagetype","lf_prob_bacphlip")]
}) %>% bind_rows()
