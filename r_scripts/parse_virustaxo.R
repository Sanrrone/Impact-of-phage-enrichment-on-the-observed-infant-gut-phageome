setwd("~/Desktop/phage_enrichment_newassembly")
 
## taxa abundances
l<-list.files(path = "taxprof_contig_virustaxo", pattern = "_virustaxo.tsv", full.names = T)
df_vtaxo<-lapply(l,function(tsv){
  sname<-gsub("_virustaxo.tsv","",basename(tsv))
  if(countLines(tsv)==0){return(NULL)}
  tmp<-read.table(tsv,header = T, sep = "\t",skipNul = T)
  colnames(tmp)<-c("ID","Length","genus","Entropy")
  tmp<-subset(tmp, Entropy<0.3) #1-0.3=0.7
  tmp$negentrophy<- 1-tmp$Entropy
  tmp$sepID<-sname
  
  ids<-getId(unique(tmp$genus),sqlFile = "/data/tmp/accessionTaxa.sql")
  taxdf<-taxonomizr::getTaxonomy(ids,sqlFile = "/data/tmp/accessionTaxa.sql")
  taxdf<-as.data.frame(taxdf)
  taxdf$tid<-rownames(taxdf)
  tmp<-merge(tmp,taxdf, by="genus", all.x = T)
  subset(tmp[,c("sepID","ID","family","genus","species","negentrophy")],!is.na(family))
  
}) %>% bind_rows()


