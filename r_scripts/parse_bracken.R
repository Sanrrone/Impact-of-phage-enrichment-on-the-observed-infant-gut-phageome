setwd("~/Desktop/phage_enrichment_newassembly/")

library(taxonomizr)
library(CHNOSZ)
library(parallel)
nc<-10
#prepareDatabase('/data/tmp/accessionTaxa.sql',getAccessions = T,tmpDir = "/data/tmp")
nodes<-getnodes("/data/tmp/HumGutDB")
names<-getnames("/data/tmp/HumGutDB")
colnames(nodes)<-c("tid","parenttid","rank")
names<-names %>% filter(type=="scientific name")
colnames(names)[1]<-c("tid")
ncbi<-merge(nodes,names[,c(1,2)], by="tid", all.x=T)

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

get_brkdf<-function(lvl="Family", mpa=TRUE, cs=0.5){
  lvl<-capitalize(lvl)
  availlvl<-c("Order","Family","Genus","Species","Strain")
  patt<-tolower(strsplit(lvl,"")[[1]][1])
  patt<-paste0("_",patt,".bracken")
  if(T){
    if(cs==0.5){
      l<-list.files(path = "prof_k2c05/", pattern = ".cs05.kreport", full.names = T) 
    }
    else{
      l<-list.files(path = "prof_k2c01", pattern = ".cs01.kreport", full.names = T) 
    }
    tlvl<-c("D"=1,"P"=2,"C"=3,"O"=4,"F"=5,"G"=6,"S"=7)
    
    get_idx<-function(x, lvl, df){
      tmpidx<-which(df[1:x,"taxlvl"]==lvl)
      ifelse(length(tmpidx)>0,max(tmpidx),1000000)
    }
    
    getlineage<-function(df){
      tdf<-lapply(unique(df$tid),function(tid){
        idx<-which(df$tid==tid)[1]
        ldf<-data.frame(tid=tid,D=NA,P=NA,C=NA,O=NA,F=NA,G=NA,S=NA)
        oriname<-df[idx,"name"]
        
        lowest_idx<-get_idx(idx,"D",df)
        for(l in tlvl){
          curr_idx<-get_idx(idx,names(tlvl[l]),df)
          if(lowest_idx<=curr_idx){
            ldf[1,names(tlvl[l])]<-df[curr_idx,"name"]
            lowest_idx<-curr_idx
          }
          if(l==5){
            break
          }
        }
        return(ldf)
      }) %>% bind_rows()
      
      colnames(tdf)<-c("tid","domain","phylum","class","order","family","genus","species")
      #ldf$species<-ifelse(is.na(ldf$family),oriname,ldf$species)
      return(tdf)
    }
    
    df_brck<-mclapply(l, function(tsv){
      #print(tsv)
      sname<-gsub(".kreport","",basename(tsv))
      tmp<-read.table(tsv,header = F, sep = "\t",skipNul = T, stringsAsFactors = F)
      colnames(tmp)<-c("relabu","rawabu","rawabulvl","taxlvl","tid","name")
      tmp$name<-trimws(tmp$name)
      maxr<-tmp[which(tmp$name=="Bacteria"),"rawabu"]
      tmp$relabu<-round((tmp$rawabu/maxr)*100,7)
      tmp<-tmp[,c("tid","taxlvl","name","rawabu","relabu")]
      tmp<-tmp %>% filter(taxlvl %in% names(tlvl))
      
      ## usar funcion propia con nodes y names y si el rank es missing usar taxonomizr.
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
  
      idxs<-sapply(1:nrow(taxdf), function(i){
        ldf<-taxdf[i,]
        idx<-which(is.na(ldf))
        if(length(idx)>0){
            expected_na<-length(idx)/length(min(idx):ncol(ldf))
            if(expected_na!=1){
              return(i)
            }
        }
      })
      idxs<-unlist(idxs)
      
      lastlvl<-sapply(idxs, function(i)taxdf[i,max(which(!is.na(taxdf[i,])))])
      lastlvl<-unlist(lastlvl)
      #sandro del futuro, ruminoccccc bacterium is NA, hacer lineage manual.      
      
      ids<-getId(lastlvl,'/data/tmp/accessionTaxa.sql')
      stillna<-idxs[which(is.na(ids))]
      nadf<-taxdf[stillna,]
      #nadf[is.na(nadf)]<-"Unknown"
      
      ids<-ids[!is.na(ids)]
      ldf<-getTaxonomy(ids,'/data/tmp/accessionTaxa.sql')
      ldf<-as.data.frame(ldf)
      ldf$tid<-rownames(ldf)
      colnames(ldf)[1]<-"domain"
      rownames(ldf)<-1:nrow(ldf)
      
      taxdf<-rbind(taxdf[-idxs,],ldf,nadf)
      
      #taxdf<-getlineage(tmp)
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
      taxdf$tid<-as.numeric(taxdf$tid)
      tmp<-merge(tmp,taxdf[,c("tid","lineage")], by="tid", all.x = T)
      tmp$sample<-sname
      
      tmp
    },mc.cores = nc) %>% bind_rows()
    
  }else{
    l<-list.files(path = "prof_bracken/", pattern = patt, full.names = T)
    df_brck<-lapply(l, function(x){
      sname<-gsub(patt,"",basename(x))
      if(sname %in% blacklist){return(NULL)}
      if(countLines(x)==0){return(NULL)}
      tmp<-read.table(x,header = T, sep = "\t", stringsAsFactors = F)
      tcount<-sum(tmp$new_est_reads)
      tmp$abundance<-round(tmp$new_est_reads/tcount,7)*100
      tmp<-tmp[,c("name","abundance")]
      colnames(tmp)<-c(tolower(lvl),"abundance")
      tmp$sample<-sname
      tmp
    }) %>% bind_rows()
  }

  
  df_brck
  
}

get_k2summary<-function(){
  dfs<-read.table("prof_bracken/summary.tsv", sep = "\t", header = F)
  dfs$V1<-sapply(dfs$V1,function(x)gsub(".kreport","",x))
  dfs
}

get_k2df_st<-function(lvl="Strain", cs=0.2){
  lvl<-capitalize(lvl)
  tlvl<-c("D"=1,"P"=2,"C"=3,"O"=4,"F"=5,"G"=6,"S"=7,"S1"=8)
  
  if(cs==0.5){
    l<-list.files(path = "prof_k2c05", pattern = ".cs05.kreport", full.names = T) 
  }
  else if(cs==0.2){
    l<-list.files(path = "prof_k2c02", pattern = ".cs02.kreport", full.names = T)
  }else if(cs==0.1){
    l<-list.files(path = "prof_k2c01", pattern = ".cs01.kreport", full.names = T)
  }else{
    l<-list.files(path = "prof_bracken", pattern = ".cs01.kreport", full.names = T) 
  }
  
  get_idx<-function(x, lvl, df){
    tmpidx<-which(df[1:x,"taxlvl"]==lvl)
    ifelse(length(tmpidx)>0,max(tmpidx),1000000)
  }
  
  getlineage<-function(df){
    tdf<-lapply(unique(df$tid),function(tid){
      idx<-which(df$tid==tid)[1]
      ldf<-data.frame(tid=tid,D=NA,P=NA,C=NA,O=NA,F=NA,G=NA,S=NA,S1=NA)
      oriname<-df[idx,"name"]
      
      lowest_idx<-get_idx(idx,"D",df)
      for(l in tlvl){
        curr_idx<-get_idx(idx,names(tlvl[l]),df)
        if(lowest_idx<=curr_idx){
          ldf[1,names(tlvl[l])]<-df[curr_idx,"name"]
          lowest_idx<-curr_idx
        }
        if(l==5){
          break
        }
      }
      return(ldf)
    }) %>% bind_rows()
    
    colnames(tdf)<-c("tid","domain","phylum","class","order","family","genus","species","strain")
    #ldf$species<-ifelse(is.na(ldf$family),oriname,ldf$species)
    return(tdf)
  }
  
  df_brck<-mclapply(l, function(tsv){
    #print(tsv)
    sname<-gsub(".kreport","",basename(tsv))
    tmp<-read.table(tsv,header = F, sep = "\t",skipNul = T, stringsAsFactors = F)
    colnames(tmp)<-c("relabu","rawabu","rawabulvl","taxlvl","tid","name")
    tmp$name<-trimws(tmp$name)
    maxr<-tmp[which(tmp$name=="Bacteria"),"rawabu"]
    tmp$relabu<-round((tmp$rawabu/maxr)*100,7)
    tmp<-tmp[,c("tid","taxlvl","name","rawabu","relabu")]
    tmp<-tmp %>% filter(taxlvl %in% names(tlvl))
    
    ## usar funcion propia con nodes y names y si el rank es missing usar taxonomizr.
    taxdf<-lapply(unique(tmp$tid), function(tid){
      lvec<-sapply(c("superkingdom","phylum","class","order","family","genus","species","strain"), 
                   function(tl){
                     climbranks(tid,tl)
                   })
      lvec<-c("tid"=tid,lvec)
      ldf<-data.frame(as.list(lvec))
      ldf$tid<-as.numeric(ldf$tid)
      return(ldf)
    }) %>% bind_rows()
    
    colnames(taxdf)[2]<-"domain"
    
    #taxdf<-getlineage(tmp)
    taxdf$species<-sapply(taxdf$species,function(x)gsub(" ","_",x))
    taxdf$domain<-ifelse(!is.na(taxdf$domain),paste0("d__",taxdf$domain),NA)
    taxdf$phylum<-ifelse(!is.na(taxdf$phylum),paste0("p__",taxdf$phylum),NA)
    taxdf$class<-ifelse(!is.na(taxdf$class),paste0("c__",taxdf$class),NA)
    taxdf$order<-ifelse(!is.na(taxdf$order),paste0("o__",taxdf$order),NA)
    taxdf$family<-ifelse(!is.na(taxdf$family),paste0("f__",taxdf$family),NA)
    taxdf$genus<-ifelse(!is.na(taxdf$genus),paste0("g__",taxdf$genus),NA)
    taxdf$species<-ifelse(!is.na(taxdf$species),paste0("s__",taxdf$species),NA)
    taxdf$strain<-ifelse(!is.na(taxdf$strain),paste0("st__",taxdf$strain),NA)
    taxdf$lineage<-apply(taxdf[,2:ncol(taxdf)], 1,function(x){
      tmp<-paste(x,collapse = "|")
      gsub("\\|NA","",tmp)
    })
    taxdf$tid<-as.numeric(taxdf$tid)
    tmp<-merge(tmp,taxdf[,c("tid","lineage")], by="tid", all.x = T)
    tmp$sample<-sname
    
    tmp
  },mc.cores = nc) %>% bind_rows()
  
  df_brck
  
}