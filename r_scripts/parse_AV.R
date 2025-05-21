setwd("~/Desktop/phage_enrichment_newassembly")
avrc<-read.csv("/data/tmp/AVrC/AvRCv1.Merged_ViralDesc.csv", header = T)
## taxa abundances
l<-list.files(path = "taxprof_contig_AV/", pattern = ".massign", full.names = T)
df_av<-lapply(l,function(tsv){
  sname<-gsub(".massign","",basename(tsv))
  if(countLines(tsv)==0){return(NULL)}
  tmp<-read.table(tsv,header = T, sep = "\t")
  tmp<-subset(tmp, identity>0.7 & coverage>0.7)
  tmp$sepID<-sname
  
  taxdf<-avrc %>% filter(contig_id %in% unique(tmp$reference))
  taxdf$reference<-taxdf$contig_id
  
  tmp<-merge(tmp,taxdf,by="reference", all.x = T)
  tmp<-tmp %>% filter(Family!="Unclassified") %>% group_by(ID) %>% top_n(1,identity+coverage)
  tmp$vc_score<-tmp$identity+tmp$coverage
  tmp[,c("sepID","ID","Family","vc_score")]
  
}) %>% bind_rows()


