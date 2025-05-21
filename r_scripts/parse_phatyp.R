setwd("~/Desktop/phage_enrichment_newassembly")
#prophage = temperate | use lysogenic cycle
#phage = virulent | use lytic stage

#lytic stage = replicate along with host genome
#lysogenic stage = replicate and assemble independent (using host resources)


l<-list.files(path = "taxprof_contig_phatyp", pattern = "_phatyp.tsv", full.names = T)
df_phatyp<-lapply(l, function(tsv){
  sname<-gsub("_phatyp.tsv","",basename(tsv))
  tmp<-read.table(tsv, header = T, sep = "\t")
  colnames(tmp)<-c("ID","Length","phagetype","lf_prob_phatyp")
  tmp<-subset(tmp,lf_prob_phatyp>0.7)
  tmp$sepID<-sname
  tmp[,c("sepID","ID","phagetype","lf_prob_phatyp")]
  
}) %>% bind_rows()
