setwd("~/Desktop/phage_enrichment_newassembly")

get_mt4df<-function(lvl="Family"){
  lvl<-capitalize(lvl)
  availlvl<-c("Order","Family","Genus","Species","Strain")
  l<-list.files(path = "prof_mt4/", pattern = ".m4report", full.names = T)
  df_mt4<-lapply(l, function(x){
    
    sname<-gsub(".m4report","",basename(x))
    if(sname %in% blacklist){return(NULL)}
    if(countLines(x)==0){return(NULL)}
    tmp<-read.table(x,header = F, sep = "\t", stringsAsFactors = F)
    tmp$V1<-sapply(tmp$V1,function(x)(gsub(",.*","",x)))
    tmp<-separate(tmp, "V1", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species","Strain"), sep="\\|", fill="right")
    tmp<-separate(tmp, "V2", c("tid_k", "tid_p", "tid_c", "tid_o", "tid_f", "tid_g", "tid_s","tid_st"), sep="\\|", fill="right")
    
    idx<-which(lvl == availlvl)
    #letter<-c("tid_o","tid_f","tid_g","tid_s")[idx]
    tmp<-tmp %>% drop_na(all_of(lvl)) %>% keep_na(all_of(availlvl[idx+1]))
    tmp<-tmp[,c(lvl,"V3")]
    tmp[[lvl]]<-sapply(tmp[[lvl]],function(f){strsplit(f,"__")[[1]][2]})
    colnames(tmp)<-c(tolower(lvl),"abundance")
    #tmp$abundance<-tmp$abundance/100
    tmp$sample<-sname
    #tmp<-subset(tmp, tid %in% utids)
    tmp
  }) %>% bind_rows()
  
  
  p1<-ggplot(subset(df_mt4, sample %in% c("HeP-311-29")), aes(x=sample, y=abundance, fill=Family)) +
    geom_bar(stat="identity", position="stack")
  
  df_mt4
  
}
