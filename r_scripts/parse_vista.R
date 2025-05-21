setwd("~/Desktop/phage_enrichment_newassembly")
l<-list.files(path = "taxprof_contig_vista/", pattern = ".tsv", full.names = T)
df_vista<-lapply(l, function(tsv){
  sname<-gsub("_vista.tsv","",basename(tsv))
  tmp<-read.table(tsv,header = F, sep = "\t", na.strings = "")
  colnames(tmp)<-c("ID","family","genus","species","distance")
  #tmp<-subset(tmp, distance<=0.1) %>% group_by(ID) %>% slice(which.min(distance)) 
  if(nrow(tmp)==0){return(NULL)}
  tmp$sepID<-sname
  tmp$simi<- 1-tmp$distance
  tmp[,c("sepID","ID","family","genus","species","simi")]
}) %>% bind_rows()

