setwd("~/Desktop/phage_enrichment_newassembly")

## taxa abundances
l<-list.files(path = "taxprof_contig_phagcn", pattern = "_phagcn.tsv", full.names = T)
df_phagcn<-lapply(l, function(tsv){
  #print(tsv)
  sname<-gsub("_phagcn.tsv","",basename(tsv))
  if(countLines(tsv)==0){return(NULL)}
  tmp<-read.table(tsv,header = T, sep = "\t",skipNul = T)
  tmp<-tmp%>%filter(PhaGCNScore!="-1",grepl(";",Lineage))
  colnames(tmp)<-tolower(colnames(tmp))
  #tmp<-tmp%>%filter(prokaryotic.virus..bacteriophages.and.archaeal.virus.=="Y",genuscluster!="singleton")
  tmp$tmpcol<-sapply(1:nrow(tmp),function(i){
    lineage_split <- strsplit(tmp[i,3], ";")[[1]]
    score_split <- strsplit(tmp[i,4], ";")[[1]]
    
    fidx<-which(grepl("family:",lineage_split)==T)
    if(length(fidx)==0){
      return("no")
    }else{
      if(length(fidx)==2){fidx<-fidx[1]} # family -> subfamily
      if(grepl("subfamily:",lineage_split[fidx])){return("no")}
      fam<-lineage_split[fidx]
      fam<-strsplit(x = fam,split = ":")[[1]][2]
      return(paste0(fam,";",score_split[fidx]))
    }
  })
  tmp<-tmp%>%filter(tmpcol!="no")
  tmp$family<-sapply(tmp$tmpcol,function(x)strsplit(x,";")[[1]][1])
  tmp$pvalue<-sapply(tmp$tmpcol,function(x)strsplit(x,";")[[1]][2])
  tmp$pvalue<-as.numeric(tmp$pvalue)
  tmp<-tmp%>%filter(pvalue>0.7)
  tmp$sepID<-sname
  tmp$species<-NA
  colnames(tmp)[1]<-"ID"
  tmp[,c("sepID","ID","family","genus","species","pvalue")]
}) %>% bind_rows()
df_phagcn$genus<-ifelse(df_phagcn$genus=="-",NA,df_phagcn$genus)

