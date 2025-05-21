setwd("~/Desktop/phage_enrichment_newassembly/")

globalscore<-0.3
humgutdf<-read.table("/data/tmp/HumGutDB/inspect.txt", sep = "\t", quote = "")
humgutdf$V6<-trimws(humgutdf$V6)
colnames(humgutdf)<-c("relabu","rawabu","rawabulvl","taxlvl","tid","name")
humgutdf<-humgutdf[,c("tid","taxlvl","name")]
tlvl<-c("D"=1,"P"=2,"C"=3,"O"=4,"F"=5,"G"=6,"S"=7)

get_idx<-function(x, lvl){
  tmpidx<-which(humgutdf[1:x,"taxlvl"]==lvl)
  ifelse(length(tmpidx)>0,max(tmpidx),1000000)
}

humgutdf_lineage<-function(tid){
  idx<-which(humgutdf$tid==tid)
  if(length(idx)==0){return(NULL)}
  ldf<-data.frame(tid=tid,D=NA,P=NA,C=NA,O=NA,F=NA,G=NA,S=NA)
  oriname<-humgutdf[idx,"name"]
  
  lowest_idx<-get_idx(idx,"D")
  for(lvl in tlvl){
    curr_idx<-get_idx(idx,names(tlvl[lvl]))
    if(lowest_idx<=curr_idx){
      ldf[1,names(tlvl[lvl])]<-humgutdf[curr_idx,"name"]
      lowest_idx<-curr_idx
    }
  }
  
  colnames(ldf)<-c("tid","domain","phylum","class","order","family","genus","species")
  ldf$species<-ifelse(is.na(ldf$species),oriname,ldf$species)
  return(ldf)
}

l<-list.files(path = "hostprof_contig_iphop/", pattern = ".csv", full.names = T)
df_iphop<-lapply(l, function(csv){
  sname<-gsub("_iphop.csv","",basename(csv))
  if(countLines(csv)<=1){return(NULL)}
  tmp<-read.table(csv,header = T, sep = ",",skipNul = T)
  tmp<-tmp%>%filter(Confidence.score>90)
  tmp$keep<-sapply(1:nrow(tmp),function(i){
    r<-tmp[i,]
    if(r['Confidence.score']<95&&!is.na(r['AAI.to.closest.RaFAH.reference'])){return("yes")}
    if(r['Confidence.score']>=95){return("yes")}
    return("no")
  })
  tmp<-tmp%>%filter(keep=="yes")
  tmp$host<-sapply(tmp$Host.genus,function(lineage){
    slin<-strsplit(lineage,";")[[1]]
    taxlvl<-slin[sapply(slin, function(x)grepl("g__",x))]
    taxlvl<-gsub("g__","",taxlvl)
    taxlvl
  })
  tmp<-subset(tmp) %>% group_by(Virus) %>% summarise(host=host[which.max(Confidence.score)], host_score=max(Confidence.score))
  tmp$sepID<-sname
  tmp$host_s<-NA
  tmp<-unique(tmp[,c("sepID","Virus","host","host_s","host_score")])
  colnames(tmp)<-c("sepID","ID","host_g","host_s","host_score")
  #tmp$host_score<-tmp$host_score/100
  tmp

}) %>% bind_rows()
head(df_iphop)


l<-list.files(path = "hostprof_contig_kr2", pattern = "_kr.tsv", full.names = T)
df_hostkr<-lapply(l,function(tsv){
  #print(tsv)
  sname<-gsub("_kr.tsv","",basename(tsv))
  if(countLines(tsv)==0){return(NULL)}
  tmp<-read.table(tsv,header = F, sep = "\t",skipNul = T)
  colnames(tmp)<-c("ID","tid","kmertax")
  tmp$tid<-as.numeric(tmp$tid)
  taxdf<-lapply(unique(tmp$tid),function(x){humgutdf_lineage(x)}) %>% bind_rows()
  
  #globalscore is 0.3(¿)
  tmp<-merge(tmp,taxdf, by="tid", all.x = T)
  tmp<-subset(tmp, kmertax>globalscore & domain=="Bacteria")[,c("ID","genus","species","kmertax")]
  tmp$kmertax<- -1*tmp$kmertax
  colnames(tmp)<-c("ID","host_g","host_s", "host_score")
  tmp<-tmp%>%filter(!is.na(host_g))
  tmp$sepID<-sname
  tmp
  
}) %>% bind_rows()

df_host<-lapply(unique(c(df_iphop$sepID, df_hostkr$sepID)),function(h){
  subtmpi<-subset(df_iphop, sepID==h)
  subtmpk<-subset(df_hostkr, sepID==h)
  tmp<-lapply(subset(cdf, sepID==h)$ID,function(id){
    idarr<-strsplit(id,"_")[[1]]
    if(length(idarr)>2){
      # >2 means phage is extracted from chromosome by checkv or phageboost.
      newid<-paste(idarr[1:2],collapse = "_")
      res<-subtmpk[which(newid == subtmpk$ID),]
      if(nrow(res)==0){res<-subtmpi[which(subtmpi$ID %in% id),]}
    }else{
      res<-subtmpi[which(subtmpi$ID %in% id),]
    }
    emptydf<-data.frame(sepID=h, ID=id, host_g=NA, host_s=NA, host_score=NA)
    if(nrow(res)==0){
      return(emptydf)
    }else{
      return(data.frame(sepID=h, ID=id, host_g=res$host_g, host_s=res$host_s, host_score=res$host_score))
    }
  }) %>% bind_rows()
  tmp
  
}) %>% bind_rows()


df_host$host_g<-sapply(df_host$host_g,function(g){
  strsplit(g,"_")[[1]][1]
})


#l<-list.files(path = "hostprof_contig_cherry", pattern = "_cherry.tsv", full.names = T)
#df_cherry<-lapply(l, function(tsv){
#  #print(tsv)
#  sname<-gsub("_cherry.tsv","",basename(tsv))
#  if(countLines(tsv)==0){return(NULL)}
#  tmp<-read.table(tsv,header = T, sep = "\t",skipNul = T, fill=T)
#  tmp<-tmp%>%filter(CHERRYScore>globalscore)
#  tmp$host_g<-sapply(tmp$Host,function(h){
#    if(grepl("species:",h)){
#      h<-gsub("species:","",h)
#      return(strsplit(h," ")[[1]][1])
#    }
#    if(grepl("genus:",h)){
#      return(gsub("genus:","",h))
#    }
#    return(NA)
#  })
#  tmp$host_s<-sapply(tmp$Host,function(h){
#    if(grepl("species:",h)){
#      return(gsub("species:","",h))
#    }
#    return(NA)
#  })
#  tmp$sepID<-sname
#  colnames(tmp)[1]<-"ID"
#  colnames(tmp)[4]<-"host_score"
#  tmp[,c("sepID","ID","host_g","host_s","host_score")]
#}) %>% bind_rows()
#
#
#df_host<-lapply(unique(c(df_cherry$sepID, df_hostkr$sepID)),function(h){
#  subtmpi<-subset(df_cherry, sepID==h)
#  subtmpk<-subset(df_hostkr, sepID==h)
#  m<-merge(subtmpi,subtmpk, by=c("sepID","ID"),all=T)
#  
#  tmp<-lapply(1:nrow(m),function(i){
#    sm<-m[i,]
#    sm$host_g<-NA
#    sm$host_s<-NA
#    sm$score<-NA
#    if(!is.na(sm[1,"host_g.x"])&&is.na(sm[1,"host_g.y"])){
#      sm[1,"host_g"]<-sm[1,"host_g.x"]
#      sm[1,"host_s"]<-sm[1,"host_s.x"]
#      sm[1,"score"]<-sm[1,"host_score.x"]
#      return(sm)
#    }
#    if(!is.na(sm[1,"host_g.y"])&&is.na(sm[1,"host_g.x"])){
#      sm[1,"host_g"]<-sm[1,"host_g.y"]
#      sm[1,"host_s"]<-sm[1,"host_s.y"]
#      sm[1,"score"]<-sm[1,"host_score.y"]
#      return(sm)
#    }
#
#    if(sm[1,"host_score.x"]>sm[1,"host_score.y"]){
#      sm[1,"host_g"]<-sm[1,"host_g.x"]
#      sm[1,"host_s"]<-sm[1,"host_s.x"]
#      sm[1,"score"]<-sm[1,"host_score.x"]
#    }else{
#      sm[1,"host_g"]<-sm[1,"host_g.y"]
#      sm[1,"host_s"]<-sm[1,"host_s.y"]
#      sm[1,"score"]<-sm[1,"host_score.y"]
#    }
#    return(sm)
#  }) %>% bind_rows()
#  tmp
#  
#}) %>% bind_rows()
#
#
#df_host$host_g<-sapply(df_host$host_g,function(g){
#  strsplit(g,"_")[[1]][1]
#})
