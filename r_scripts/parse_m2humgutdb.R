setwd("~/Desktop/phage_enrichment_newassembly")

library(taxonomizr)
library(CHNOSZ)
library(dplyr)
library(data.table)
library(parallel)
#numCores <- detectCores()
numCores<-12

get_m2df<-function(){
  
  nodes<-getnodes("/data/tmp/HumGutDB")
  names<-getnames("/data/tmp/HumGutDB")
  colnames(nodes)<-c("tid","parenttid","rank")
  names<-names %>% filter(type=="scientific name")
  colnames(names)[1]<-c("tid")
  ncbi<-merge(nodes,names[,c(1,2)], by="tid", all.x=T)
  #prepareDatabase('/data/tmp/accessionTaxa.sql',getAccessions = T,tmpDir = "/data/tmp")
  #humtsv<-read.table("/data/tmp/HumGutDB/HumGut.tsv",header = T,sep = "\t")
  
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
  
  l<-list.files(path = "prof_m2p", pattern = ".bactstats", full.names = T)
  #totaltsv<-length(l)
  df_m2<-mclapply(1:length(l), function(i){
    tsv<-l[i]
    #print(paste0("Parsing file ",i))
    sname<-gsub(".bactstats","",basename(tsv))
    tmp<-read.table(tsv,header = T, sep = "\t",skipNul = T, stringsAsFactors = F)
    if(nrow(tmp)==0){stop(paste(sname, "wait, no rows in file :OOO"))}
    totalab<-sdf[which(sdf$sample==sname),"total"]
    if(length(totalab)==0){print(paste(sname, "not found in sdf"));return(NULL)}
    sum(tmp$raw_abundance)
    tmp$rel_abu<-round((tmp$raw_abundance/totalab),8)
    tmp$cov<- round(tmp$bpcov1/tmp$length,7)
    tmp$cov3<- round(tmp$bpcov3/tmp$length,7)
    tmp$cov5<- round(tmp$bpcov5/tmp$length,7)
    tmp$cov10<- round(tmp$bpcov10/tmp$length,7)
    tmp$mdcov<- round(tmp$mdepth1*tmp$cov,7)
    tmp$mdcov3<- round(tmp$mdepth3*tmp$cov3,7)
    tmp$mdcov5<- round(tmp$mdepth5*tmp$cov5,7)
    tmp$mdcov10<- round(tmp$mdepth10*tmp$cov10,7)

    
    tmp$tid<-sapply(tmp$ID,function(id)as.numeric(strsplit(id,split = "_")[[1]][1]))
    #tmp$species<-sapply(tmp$tid,function(tid)climbranks(tid,"species"))
    
    taxdf<-lapply(unique(tmp$tid), function(tid){
      lvec<-sapply(c("superkingdom","phylum","class","order","family","genus","species"), 
                   function(tl){
                     climbranks(tid,tl)
                   })
      lvec<-c("tid"=tid,lvec)
      ldf<-data.frame(as.list(lvec))
      ldf$tid<-as.numeric(ldf$tid)
      return(ldf)
    }) %>% bind_rows()
    
    colnames(taxdf)[2]<-"domain"
    
    taxdf$species<-sapply(taxdf$species,function(x)gsub(" ","_",x))
    taxdf$domain<-ifelse(!is.na(taxdf$domain),paste0("d__",taxdf$domain),NA)
    taxdf$phylum<-ifelse(!is.na(taxdf$phylum),paste0("p__",taxdf$phylum),NA)
    taxdf$class<-ifelse(!is.na(taxdf$class),paste0("c__",taxdf$class),NA)
    taxdf$order<-ifelse(!is.na(taxdf$order),paste0("o__",taxdf$order),NA)
    taxdf$family<-ifelse(!is.na(taxdf$family),paste0("f__",taxdf$family),NA)
    taxdf$genus<-ifelse(!is.na(taxdf$genus),paste0("g__",taxdf$genus),NA)
    taxdf$species<-ifelse(!is.na(taxdf$species),paste0("s__",taxdf$species),NA)
    taxdf$lineage<-apply(taxdf[,2:ncol(taxdf)], 1,function(x){
      tmp<-paste(x,collapse = "|")
      gsub("\\|NA","",tmp)
    })
    
    tmp<-merge(tmp,taxdf[,c("tid","lineage","species")], by="tid", all.x = T)
    tmp$species<-sapply(tmp$species,function(x){gsub("s__","",x)})
    tmp<-tmp %>% group_by(lineage,species) %>% reframe(raw_abu=sum(raw_abundance),
                                               mdcov=mean(mdcov),
                                               mdcov3=mean(mdcov3),
                                               mdcov5=mean(mdcov5),
                                               mdcov10=mean(mdcov10))
    mdcsum<-sum(tmp$mdcov)
    mdcsum3<-sum(tmp$mdcov3)
    mdcsum5<-sum(tmp$mdcov5)
    mdcsum10<-sum(tmp$mdcov10)
    tmp$rel_mdcov<-tmp$mdcov/ifelse(mdcsum==0,1,mdcsum)
    tmp$rel_mdcov3<-tmp$mdcov3/ifelse(mdcsum3==0,1,mdcsum3)
    tmp$rel_mdcov5<-tmp$mdcov5/ifelse(mdcsum5==0,1,mdcsum5)
    tmp$rel_mdcov10<-tmp$mdcov10/ifelse(mdcsum10==0,1,mdcsum10)

    tmp$sample<-sname
    tmp
  },mc.cores = numCores) %>% bind_rows()
  
}
