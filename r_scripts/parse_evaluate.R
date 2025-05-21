#dvp
l<-list.files(path = "evaluate", pattern = "_dvp_vscores.tsv", full.names = T)
df_dvp<-lapply(l, function(tsv){
  #print(tsv)
  sname<-gsub("_dvp_vscores.tsv","",basename(tsv))
  if(countLines(tsv)<=1){return(NULL)}
  tmp<-read.table(tsv,header = T, sep = "\t",skipNul = T, stringsAsFactors = F)
  colnames(tmp)<-c("ID","len","dvp","pbb")
  tmp$hepID<-sname
  tmp[,c("ID","dvp","hepID")]
}) %>% bind_rows()

#vs2
l<-list.files(path = "evaluate", pattern = "_vs2_vscores.tsv", full.names = T)
df_vs2<-lapply(l, function(tsv){
  #print(tsv)
  sname<-gsub("_vs2_vscores.tsv","",basename(tsv))
  if(countLines(tsv)<=1){return(NULL)}
  tmp<-read.table(tsv,header = T, sep = "\t", stringsAsFactors = F)
  colnames(tmp)<-c("ID","vs2","max_score","max_score_group","group","length","hallmark","viral_gene","nonviral_gene")
  tmp$ID<-sapply(tmp$ID,function(i){strsplit(i,"\\|\\|")[[1]][1]})
  tmp$hepID<-sname
  tmp$vs2<-ifelse(is.nan(tmp$vs2),0,tmp$vs2)
  tmp[,c("ID","vs2","hepID")]
}) %>% bind_rows()

#mpp
l<-list.files(path = "evaluate", pattern = "_mpp_vscores.tsv", full.names = T)
df_mpp<-lapply(l, function(tsv){
  #print(tsv)
  sname<-gsub("_mpp_vscores.tsv","",basename(tsv))
  if(countLines(tsv)<=1){return(NULL)}
  tmp<-read.table(tsv,header = T, sep = "\t", stringsAsFactors = F)
  colnames(tmp)<-c("ID","length","mpp")
  tmp$hepID<-sname
  tmp[,c("ID","mpp","hepID")]
}) %>% bind_rows()

#phamer
l<-list.files(path = "evaluate", pattern = "_phamer_vscores.tsv", full.names = T)
df_phamer<-lapply(l, function(tsv){
  #print(tsv)
  sname<-gsub("_phamer_vscores.tsv","",basename(tsv))
  if(countLines(tsv)<=1){return(NULL)}
  tmp<-read.table(tsv,header = T, sep = "\t", stringsAsFactors = F)
  colnames(tmp)<-c("ID","length","pred","proportion","phamer","confidence")
  tmp$hepID<-sname
  tmp[,c("ID","phamer","hepID")]
}) %>% bind_rows()

df_eval<-Reduce(function(x, y) merge(x,y, all= TRUE), list(df_dvp,df_vs2,df_mpp,df_phamer)) %>% tibble()
df_eval$vs2<-ifelse(is.na(df_eval$vs2),0,df_eval$vs2)
df_eval$sum_score<-df_eval$dvp+df_eval$vs2+df_eval$mpp+df_eval$phamer
df_eval$mean_score<-(df_eval$dvp+df_eval$vs2+df_eval$mpp+df_eval$phamer)/4

write.table(df_eval,"df_eval.tsv",quote = F, sep = "\t",row.names = F, col.names = T)


p1<- ggplot(df_eval,aes(sum_score)) +
  geom_density(fill="green",alpha=0.2,linewidth = 0.2, color="darkgreen") +
  theme_minimal() + 
  geom_vline(xintercept = 1,linetype="dashed") +
  annotate(geom = "text",x = 1.2, y=0.2, label=paste0("Minimum cutoff (~20% filtered)"), angle=-90, size=3.5) +
  geom_vline(xintercept = 2,linetype="dashed") +
  annotate(geom = "text",x = 2.2, y=0.2, label=paste0("Medium cutoff (~60% filtered)"), angle=-90, size=3.5) +
  geom_vline(xintercept = 3,linetype="dashed") +
  annotate(geom = "text",x = 3.2, y=0.2, label=paste0("High cutoff (~80% filtered)"), angle=-90, size=3.5) +
  ggtitle("Distribution of viral score for HeP viruses") +
  xlab("Score") + ylab("Density")
#p1



