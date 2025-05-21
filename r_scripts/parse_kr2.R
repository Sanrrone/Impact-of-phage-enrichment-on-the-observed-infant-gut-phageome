library(CHNOSZ)
library(taxonomizr)
setwd("~/Desktop/phage_enrichment_newassembly")
uhgvdf<-read.table("/data/tmp/phantaDB/phantaDB.out", sep = "\t", quote = "")
uhgvdf$V6<-trimws(uhgvdf$V6)
colnames(uhgvdf)<-c("relabu","rawabu","rawabulvl","taxlvl","tid","name")
uhgvdf<-uhgvdf[,c("tid","taxlvl","name")]
phantanames<-getnames("/data/tmp/phantaDB")
phantanodes<-getnodes("/data/tmp/phantaDB")
colnames(phantanodes)<-c("tid","parenttid","rank")
phantanames<-phantanames %>% filter(type=="scientific name")
colnames(phantanames)[1]<-"tid"
ncbi<-merge(phantanodes,phantanames[,c(1,2)], by="tid", all.x=T) %>% unique()

climbranks<-function(tid,lvl){
  idx<-which(ncbi$tid==tid)
  if(length(idx)==0){
    return(NA)
  }
  if(ncbi[idx,"tid"]==1){
    return(NA)
  }
  if(ncbi[idx,"rank"]==lvl){
    return(ncbi[idx,"name"])
  }else{
    return(climbranks(ncbi[idx,"parenttid"],lvl))
  }
}


tlvl<-c("D"=1,"P"=2,"C"=3,"O"=4,"F"=5,"G"=6,"S"=7)

uhgv_metadata<-read.table("/data/tmp/phantaDB/votus_metadata_extended.tsv", sep = "\t", quote = "", header = T, na.strings = "NULL")
uhgv_metadata$kmers_host_lineage_gtdb_r207<-ifelse(uhgv_metadata$kmers_host_agreement>80,
                                                   uhgv_metadata$kmers_host_lineage_gtdb_r207,NA)
uhgv_metadata$host_lineage_gtdb_r207<-ifelse(is.na(uhgv_metadata$crispr_host_lineage_gtdb_r207)&
                                               is.na(uhgv_metadata$kmers_host_lineage_gtdb_r207),
                                               NA,uhgv_metadata$host_lineage_gtdb_r207)

uhgv_metadata<-uhgv_metadata%>%filter(checkv_quality %in% c("Complete","High-quality","Medium-quality"))                                             
uhgv_metadata<-uhgv_metadata[,c("uhgv_genome","uhgv_taxonomy","uhgv_votu","family_name",
                     "checkv_quality",
                     "class_name", "host_lineage_gtdb_r207",
                     "kmers_host_lineage_gtdb_r207",
                     "crispr_host_lineage_gtdb_r207")]
uhgv_metadata$vfam<-sapply(uhgv_metadata$uhgv_taxonomy,function(n){
  sn<-strsplit(x = n,split = ";")[[1]]
  idx<-which(grepl(pattern = "vFAM-",sn))
  if(length(idx)==0){
    return(NA)
  }else{
    return(sn[idx])
  }
})
uhgv_metadata$isphage<-sapply(1:nrow(uhgv_metadata),function(i){
  row<-uhgv_metadata[i,]
  ifelse("Caudoviricetes"==row[,"class_name"] || 
           grepl("d__Bacteria;",row[,"host_lineage_gtdb_r207"]) ||
           grepl("d__Bacteria;",row[,"kmers_host_lineage_gtdb_r207"]) ||
           grepl("d__Bacteria;",row[,"crispr_host_lineage_gtdb_r207"]),
         "yes",NA)
})
uhgv_metadata<-uhgv_metadata[,c("uhgv_genome","uhgv_votu","vfam","class_name",
                                    "family_name","uhgv_taxonomy","isphage")]

get_idx<-function(x, lvl){
  tmpidx<-which(uhgvdf[1:x,"taxlvl"]==lvl)
  ifelse(length(tmpidx)>0,max(tmpidx),1000000)
}

phanta_lineage<-function(tid){
  idx<-which(uhgvdf$tid==tid)
  ldf<-data.frame(tid=tid,D=NA,P=NA,C=NA,O=NA,F=NA,G=NA,S=NA)
  oriname<-uhgvdf[idx,"name"]

  lowest_idx<-get_idx(idx,"D")
  for(lvl in tlvl){
    curr_idx<-get_idx(idx,names(tlvl[lvl]))
    if(lowest_idx<=curr_idx){
      ldf[1,names(tlvl[lvl])]<-uhgvdf[curr_idx,"name"]
      lowest_idx<-curr_idx
    }
  }
  
  colnames(ldf)<-c("tid","domain","phylum","class","order","family","genus","species")
  ldf$oriname<-oriname
  return(ldf)
}

## taxa abundances
## contigs
l<-list.files(path = "taxprof_contig_kr2/", pattern = ".ktsv", full.names = T)
df_kr2<-lapply(l, function(tsv){
  #print(tsv)
  sname<-gsub(".ktsv","",basename(tsv))
  if(countLines(tsv)==0){return(NULL)}
  tmp<-read.table(tsv,header = F, sep = "\t",skipNul = T, stringsAsFactors = F)
  colnames(tmp)<-c("ID","tid","kmertax")
  tmp<-subset(tmp, round(kmertax,2)>=0.7)
  if(nrow(tmp)==0){
    return(NULL)
  }
  #taxdf<-lapply(unique(tmp$tid), function(x)phanta_lineage(x)) %>% bind_rows()
  
  taxdf<-lapply(unique(tmp$tid), function(tid){
    if(tid == 1 || tid == 131567 ){return(NULL)}
    lvec<-sapply(c("superkingdom","phylum","class","order","family","genus","species"), 
                 function(tl){
                   climbranks(tid,tl)
                 })
    lvec<-c("tid"=tid,lvec)
    ldf<-data.frame(as.list(lvec))
    ldf$tid<-as.numeric(ldf$tid)
    ldf$oriname<-ncbi[which(ncbi$tid==tid),"name"]
    return(ldf)
  }) %>% bind_rows()
  
  colnames(taxdf)[2]<-"domain"
  
  taxdf<-subset(taxdf,domain=="Viruses")
  
  tmp<-merge(tmp,taxdf, by="tid", all.x = T)
  tmp$hepID<-sname
  tmp$family<-sapply(tmp$family,function(x){
    #print(x)
    if(is.na(x)){return(NA)}
    idx<-which(uhgv_metadata$vfam==x)
    fn<-unique(uhgv_metadata[idx,]$family_name)
    fn<-fn[!is.na(fn)]
    if(length(fn)==0){return(NA)}
    if(length(fn)==1){return(fn)}
    #if(length(fn)>1){return(paste(fn,collapse = ";"))}
    if(length(fn)>1){return(NA)}
    
  })
  tmp$genus<-sapply(tmp$genus,function(x){
    #print(x)
    if(is.na(x)){return(NA)}
    idx<-which(uhgv_metadata$vfam==x)
    fn<-unique(uhgv_metadata[idx,]$family_name)
    fn<-fn[!is.na(fn)]
    if(length(fn)==0){return(NA)}
    if(length(fn)==1){return(fn)}
    if(length(fn)>1){return(paste(fn,collapse = ";"))}
  })
  tmp$isphage<-sapply(1:nrow(tmp),function(i){ #species works also as original name
    row<-tmp[i,]
    if(is.na(row[,"oriname"])){return(NA)}
    idx<-which(row[,"class"]==uhgv_metadata$class_name)
    if(length(idx)>0){
      ip<-unique(uhgv_metadata[idx,"isphage"])
      ip<-ip[!is.na(ip)]
      if(length(ip)>0)return(ip)
    }
    idx<-which(row[,"oriname"]==uhgv_metadata$vfam)
    if(length(idx)>0){
      ip<-unique(uhgv_metadata[idx,"isphage"]);return(ip[!is.na(ip)])
      ip<-ip[!is.na(ip)]
      if(length(ip)>0){return(ip)}else{return(NA)}
    }
    idx<-which(row[,"oriname"]==uhgv_metadata$uhgv_votu)
    if(length(idx)>0){
      ip<-unique(uhgv_metadata[idx,"isphage"])
      ip<-ip[!is.na(ip)]
      if(length(ip)>0){return(ip)}else{return(NA)}
    }
    idx<-which(row[,"oriname"]==uhgv_metadata$uhgv_genome)
    if(length(idx)>0){
      ip<-unique(uhgv_metadata[idx,"isphage"])
      ip<-ip[!is.na(ip)]
      if(length(ip)>0){return(ip)}else{return(NA)}
    }
    ip<-subset(uhgv_metadata, grepl(row[,"oriname"],uhgv_taxonomy))[,"isphage"]
    if(length(ip)>0){ip<-ip[!is.na(ip)];if(length(ip)>0)return(ip)}
    return(NA)
  })
  
  
  tmp<-subset(tmp[,c("hepID","ID","isphage","family","genus","species","kmertax")], !is.na(family))
  if(nrow(tmp)==0){
    return(NULL)
  }else{
    return(tmp)
  }
}) %>% bind_rows()

dfphagetag<-df_kr2[,c("hepID","ID","isphage")]
df_kr2$isphage<-NULL

#### general plot
p1<-ggplot(df_kr2, aes(x=hepID, y=kmertax)) +
  geom_boxplot() + ylim(0,1) +
  coord_flip() + theme_bw() +
  geom_hline(yintercept = 0.5, linetype="dashed", color="darkgreen") +
  ylab("K2 Confidence score") + xlab("HeP family") +
  ggtitle("Confidence score for Contigs from HeP samples")
####################################

