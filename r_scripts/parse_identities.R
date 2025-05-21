setwd("~/Desktop/phage_enrichment_newassembly")

#l<-list.files(path = "identities", pattern = ".idt", full.names = T)
#df_idn<-lapply(l, function(tsv){
#  sname<-gsub(".idt","",basename(tsv))
#  tmp<-read.table(tsv, header = F, sep = "\t")
#  colnames(tmp)<-c("ID1","ID2","identity")
#  tmp$hepID<-sname
#  tmp %>% filter(identity>0.2)
#  
#}) %>% bind_rows()
#df_idn$identity<-round(df_idn$identity,4)
#df_idn<-df_idn %>% filter(ID1!=ID2)

#l<-list.files(path = "bins_aai", pattern = "_ani_norank.tsv", full.names = T)
#df_idnani<-lapply(l, function(tsv){
#  sname<-gsub("_ani_norank.tsv","",basename(tsv))
#  sname<-gsub("aai_idn_","",sname)
#  if(countLines(tsv)==0){return(NULL)}
#  tmp<-read.table(tsv, header = F, sep = "\t")
#  colnames(tmp)<-c("ID1","ID2","identity")
#  tmp$hepID<-sname
#  tmp
#  
#}) %>% bind_rows()



#l<-list.files(path = "identities",pattern = "_family_edges.tsv", full.names = T)
#aaiidn_df<-lapply(l, function(tsv){
#  sname<-gsub("_family_edges.tsv","",basename(tsv))
#  sname<-gsub("aai_aaidn_","",sname)
#  if(countLines(tsv)==0){return(NULL)}
#  tmp<-read.table(tsv, header = F, sep = "\t")
#  colnames(tmp)<-c("ID1","ID2","identity")
#  tmp$hepID<-sname
#  tmp %>% filter(ID1!=ID2)
#  
#}) %>% bind_rows() %>% tibble()
#
#l<-list.files(path = "identities",pattern = "_maniac.csv", full.names = T)
#maniac_df<-lapply(l, function(csv){
#  sname<-gsub("_maniac.csv","",basename(csv))
#  if(countLines(csv)==0){return(NULL)}
#  tmp<-read.csv(csv, header = T)
#  colnames(tmp)[c(1:2)]<-c("ID1","ID2")
#  tmp$hepID<-sname
#  tmp<-tmp[,c(1,2,13,14)]
#  tmp %>% filter(ID1!=ID2)
#  
#  
#}) %>% bind_rows() %>% tibble()

l<-list.files(path = "atclust",pattern = "_norank.repatclust", full.names = T)
atclust_df<-lapply(l, function(csv){
  sname<-gsub("_norank.repatclust","",basename(csv))
  if(countLines(csv)==0){return(NULL)}
  tmp<-read.csv(csv, header = T)
  colnames(tmp)[c(3,5,6)]<-c("identity","ID1","ID2")
  tmp$hepID<-sname
  tmp<-tmp[,c(5,6,3,7)]
  tmp %>% filter(ID1!=ID2)
  
  
}) %>% bind_rows() %>% tibble()

#cluster by ATClust (kmer and shanon entropy based)
fid<-0.7
cls_df_atclust<-lapply(heps, function(h){
  #print(h)
  ftmp<-atclust_df%>%filter(identity>fid, hepID==h)
  if(nrow(ftmp)==0){return(NULL)}
  dtmp<-dcast(ftmp, ID1~ID2,value.var = "identity",fill=0)
  rownames(dtmp)<-dtmp$ID1
  dtmp$ID1<-NULL
  
  if(length(colnames(dtmp))<=1){return(data.frame(ID=c(ftmp$ID1,ftmp$ID2),bin=1, row.names = NULL))}
  hc<-hclust(as.dist(1-dtmp), method = "ward.D2")
  mhc<-max(hc$height, na.rm = TRUE)
  hc$height <- hc$height / ifelse(mhc==0,1,mhc)
  if (is.unsorted(hc$height)) hc$height <- sort(hc$height)
  #plot(hc)
  clusters <- cutree(hc, h = fid)
  
  gdf<-data.frame(ID=names(clusters),bin=clusters, row.names = NULL)
  gdf$hepID<-h
  gdf%>%arrange(bin)
  
}) %>% bind_rows() %>% tibble()

