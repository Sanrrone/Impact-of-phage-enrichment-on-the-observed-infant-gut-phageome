source("functions/essentials.R")
options(scipen = 999)

if(file.exists("dfcounts.tsv")){
  df_rcounts<-read.table("dfcounts.tsv",sep="\t",header=T)
  print("**** File df_rcounts already exist (lodaded) *****")
}else{
  l<-list.files(path = "raw_counts_contigs/", pattern = ".cov", full.names = T)
  df_rcounts<-lapply(l, function(x){
    #print(x)
    sname<-gsub("_map.cov","",basename(x))
    if(sname %in% blacklist){return(NULL)}
    if(countLines(x)==0){return(NULL)}
    tmp<-read.table(x,sep = "\t",header = T)
    tmp$sample<-sname
    totalab<-sdf[which(sdf$sample==sname),"total"]
    if(length(totalab)==0){print(paste(sname, "not found in sdf"));return(NULL)}
    tmp$rel_abundance<-round(tmp$raw_abundance/totalab,8)
    
    tmp$cov<- round(tmp$bpcov1/tmp$length,7)
    tmp$cov3<- round(tmp$bpcov3/tmp$length,7)
    tmp$cov5<- round(tmp$bpcov5/tmp$length,7)
    tmp$cov10<- round(tmp$bpcov10/tmp$length,7)
    tmp$mdcov<- round(tmp$mdepth1*tmp$cov,7)
    tmp$mdcov3<- round(tmp$mdepth3*tmp$cov3,7)
    tmp$mdcov5<- round(tmp$mdepth5*tmp$cov5,7)
    tmp$mdcov10<- round(tmp$mdepth10*tmp$cov10,7)
    
    mdcsum<-sum(tmp$mdcov)
    mdcsum<-ifelse(mdcsum==0,1,mdcsum)
    mdcsum3<-sum(tmp$mdcov3)
    mdcsum3<-ifelse(mdcsum3==0,1,mdcsum3)
    mdcsum5<-sum(tmp$mdcov5)
    mdcsum5<-ifelse(mdcsum5==0,1,mdcsum5)
    mdcsum10<-sum(tmp$mdcov10)
    mdcsum10<-ifelse(mdcsum10==0,1,mdcsum10)
    
    tmp$rel_mdcov<-tmp$mdcov/mdcsum
    tmp$rel_mdcov3<-tmp$mdcov3/mdcsum3
    tmp$rel_mdcov5<-tmp$mdcov5/mdcsum5
    tmp$rel_mdcov10<-tmp$mdcov10/mdcsum10
    
    tmp
  }) %>% bind_rows
  df_rcounts<-add_metadata(df_rcounts)
  df_rcounts<-lapply(unique(df_rcounts$hepID),function(h){
    tmp0<-df_rcounts %>% filter(hepID==h)
    ddf<-dcast(tmp0 %>% group_by(sample,ID) %>% reframe(pa=ifelse(cov>0,1,0)),
               ID~sample,value.var = "pa", fill = 0, fun.aggregate = mean)
    ddf$prevalence<-rowSums(ddf[,-1])/ncol(ddf[,-1])
    ddf$hepID<-h
    ddf<-ddf[,c("hepID","ID","prevalence")]
    mprevs<-merge(tmp0,ddf, by=c("hepID","ID"), all=T)
    
    return(mprevs)
    
  }) %>% bind_rows()
  write.table(df_rcounts,"dfcounts.tsv", sep = "\t", quote = F, row.names = F, col.names = T)
  
  #generate propgage coordinates.
  
  smdf<-df_rcounts[,c("hepID","sample","ID")] %>% filter(grepl(" ",ID)|grepl("phage",ID))
  smdf$scaffold<-sapply(smdf$ID,function(x){paste(strsplit(x,"_")[[1]][1:2], collapse = "_")})
  smdf$fragment<-smdf$ID
  smdf$start<-sapply(smdf$ID,function(id){
    if(grepl(" ",id)){
      tmp<-strsplit(id," ")[[1]][2]
      tmp<-strsplit(tmp,"/")[[1]][1]
      tmp<-strsplit(tmp,"-")[[1]][1]
      return(tmp)
    }
    if(grepl("phage",id)){
      tmp<-strsplit(id,"_")[[1]][4]
      return(tmp)
    }
  })
  smdf$stop<-sapply(smdf$ID,function(id){
    if(grepl(" ",id)){
      tmp<-strsplit(id," ")[[1]][2]
      tmp<-strsplit(tmp,"/")[[1]][1]
      tmp<-strsplit(tmp,"-")[[1]][2]
      return(tmp)
    }
    if(grepl("phage",id)){
      tmp<-strsplit(id,"_")[[1]][5]
      return(tmp)
    }
  })
  smdf$ID<-NULL
  write.table(smdf,"helpfiles/prophages_coordinates.tsv",sep = "\t",quote = F,row.names = F)
}

add_abundances<-function(df){
  tmp<-merge(df_rcounts, df, by=c("ID","hepID"), all.x = T)
  tmp$family<-ifelse(is.na(tmp$family),"Viral",tmp$family)
  tmp
  
}

getfamily<-function(tidvec){
  #ids<-getId(namevec,sqlFile = "/data/tmp/accessionTaxa.sql")
  taxdf<-taxonomizr::getTaxonomy(tidvec,sqlFile = "/data/tmp/accessionTaxa.sql",desiredTaxa = c("family","genus","species"))
  taxdf<-as.data.frame(taxdf)
  taxdf
  
}

rename_cluster<-function(gdf){
  gs<-unique(gdf$group)
  ngs<-1:length(gs)
  names(ngs)<-gs
  gdf$group<-ngs[as.character(gdf$group)]
  gdf %>% arrange(group)
}
