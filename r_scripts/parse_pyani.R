l<-list.files(path = "pyani", pattern = "_percentage_identity.tab",full.names = T)
id_df<-lapply(l, function(i){
  sname<-gsub(pattern = "_percentage_identity.tab","",basename(l))
  tmp<-read.table(i, sep = "\t",header = T, check.names = F)
  tmp<-melt(tmp)
  colnames(tmp)<-c("ID_1","ID_2","identity")
  tmp
}) %>% bind_rows()


l<-list.files(path = "pyani", pattern = "_alignment_coverage.tab",full.names = T)
cov_df<-lapply(l, function(i){
  sname<-gsub(pattern = "_alignment_coverage.tab","",basename(l))
  tmp<-read.table(i, sep = "\t",header = T, check.names = F)
  tmp<-melt(tmp)
  colnames(tmp)<-c("ID_1","ID_2","coverage")
  tmp
}) %>% bind_rows()

pyani_df<-merge(id_df, cov_df, by=c("ID_1","ID_2"), all = T)

