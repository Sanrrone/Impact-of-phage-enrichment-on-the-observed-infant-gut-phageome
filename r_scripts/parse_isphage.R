setwd("~/Desktop/phage_enrichment_newassembly")

l<-list.files(path = "taxprof_phage_tag", pattern = "_mpp_exclusive.txt", full.names = T)
df_isp1<-lapply(l,function(txt){
  sname<-gsub("_mpp_exclusive.txt","",basename(txt))
  txtdf<-read.table(txt,header = T) %>% unique()
  if(nrow(txtdf)==0){return(NULL)}
  colnames(txtdf)[1]<-"ID"
  txtdf$isphage<-"yes"
  txtdf$sepID<-sname
  txtdf
  
}) %>% bind_rows()

l<-list.files(path = "taxprof_phage_tag", pattern = "_phamer_exclusive.txt", full.names = T)
df_isp2<-lapply(l,function(txt){
  sname<-gsub("_phamer_exclusive.txt","",basename(txt))
  txtdf<-read.table(txt,header = T) %>% unique()
  if(nrow(txtdf)==0){return(NULL)}
  colnames(txtdf)[1]<-"ID"
  txtdf$isphage<-"yes"
  txtdf$sepID<-sname
  txtdf
  
}) %>% bind_rows()

df_isphage<-merge(df_isp1,df_isp2,by=c("sepID","ID"))
df_isphage$isphage<-sapply(1:nrow(df_isphage),function(i){
  x<-df_isphage[i,] %>% pull(isphage.x)
  y<-df_isphage[i,] %>% pull(isphage.y)
  ifelse(x=="yes" || y=="yes", "yes",NA)
})
df_isphage<-df_isphage[,c("sepID","ID","isphage")]
