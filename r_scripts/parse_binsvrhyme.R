l<-list.files(path = "bins_vRhyme/",pattern = "_vrhyme", full.names = T)
vbins_df<-lapply(l,function(tsv){
  sname<-gsub(pattern = "_vrhyme.tsv", replacement = "", x=basename(tsv))
  tmp<-read.table(tsv, sep = "\t", header = T)
  colnames(tmp)<-c("ID","bin")
  tmp$sepID<-sname
  tmp
}) %>% bind_rows()
