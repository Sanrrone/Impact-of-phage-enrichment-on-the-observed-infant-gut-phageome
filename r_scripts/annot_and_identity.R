setwd("~/Desktop/phage_enrichment_newassembly")
source("functions/parse_hostpred.R")
source("functions/parse_taxannot.R")
source("functions/parse_binsvrhyme.R")
#source("functions/parse_binsaai.R")
source("functions/parse_phatyp.R")
source("functions/parse_identities.R")
source("functions/parse_evaluate.R")

#if(T) if we test

#cluster host
if(F){
  if(file.exists("hbins_p7_tc7_ac85_cv7.tsv")){
    hbins<-tibble(read.table("hbins_p7_tc7_ac85_cv7.tsv",sep="\t",header=T))
    print("**** File hbins.tsv already exist (lodaded) *****")
  }else{
    source("../viral_strain_tracking/refine_bins.R")
    hbins<-mclapply(unique(df_rcounts$sepID),function(h){
      #print(h)
      allowed<-df_eval%>%filter(sepID==h)
      tmp<-subset(df_rcounts, sepID==h & ID%in%allowed$ID)
      htmp<-binbypatterns(tmp,prevalence = 0.7,prune_tree = 0.7,agreement = 0.85,
                          magnitud_diff = 0.5,min_corr = c(0.7))
      htmp$sepID<-h
      htmp
    },mc.cores = 2) %>% bind_rows()
    write.table(hbins, "hbins_p7_tc7_ac85_cv7.tsv", sep = "\t", quote = F, row.names = F)
  }
}


## merge clusters methods
head(vbins_df)
head(cls_df_atclust)

#min identity cluster
if(F){
  df_pair<-lapply(seps, function(h){
    tdf<-tax_df %>% filter(sepID==h, !is.na(family))
    tdf$ID1<-tdf$ID
    tdf$ID<-NULL
    #tmp<-df_idn %>% filter(sepID==h)
    tmp<-atclust_df %>% filter(sepID==h)
    #tmp<-maniac_df %>% filter(sepID==h)
    #tmp<-atclust_df %>% filter(sepID==h)
    if("wgANI" %in% colnames(tmp)){
      ftmp<-tmp%>%filter(wgANI>0)
      icol<-"wgANI"
    }else{
      ftmp<-tmp%>%filter(identity>0)
      icol<-"identity"
    }
    ftmp<-merge(ftmp,tdf[,c("sepID","ID1","family")], by=c("sepID","ID1"),all.x = T)
    ftmp<-ftmp %>% filter(!is.na(family))
    colnames(ftmp)[which(colnames(ftmp)=="family")]<-"family_ID1"
    tdf$ID2<-tdf$ID1
    tdf$ID1<-NULL
    ftmp<-merge(ftmp,tdf[,c("sepID","ID2","family")], by=c("sepID","ID2"),all.x = T)
    ftmp<-ftmp %>% filter(!is.na(family_ID1),!is.na(family))
    colnames(ftmp)[which(colnames(ftmp)=="family")]<-"family_ID2"
    ftmp[,c("sepID","ID1","ID2",icol,"family_ID1","family_ID2")]
    
  }) %>% bind_rows() %>% tibble()
  df_pair$Annotation<-ifelse(df_pair$family_ID1==df_pair$family_ID2,"Same family","Different")
  #write.table(df_pair,"identities.tsv",sep = "\t",quote = F,row.names = F,col.names = T)
  #df_pair<-read.table("identities.tsv", sep = "\t",header = T)
  minid<-df_pair%>%filter(Annotation=="Different")%>%pull(identity)%>%max()
  
  ggplot(df_pair,aes(identity)) +
    geom_density(aes(fill=Annotation, color=Annotation), alpha=0.2,linewidth = 0.2) +
    theme_minimal() + 
    geom_vline(xintercept = minid,linetype="dashed") +
    annotate(geom = "text",x = minid+0.02, y=7, label=paste0("Agreement limit: ",minid),  angle=-90,size=3.5) +
    ggtitle("ALFATClust score between pairs")
}



cls_df<-lapply(seps, function(h){
  vtmp<-subset(vbins_df,sepID==h)
  atmp<-subset(cls_df_atclust, sepID==h)
  
  idv<-vtmp%>%pull(ID)%>%unique()
  id_notc<-atmp %>% filter(!ID%in%idv)
  colnames(id_notc)[2]<-"group"
  id_notc<-rename_cluster(id_notc)
  colnames(id_notc)[2]<-"bin"
  id_notc$bin<-id_notc$bin+max(vtmp$bin)
  rbind(vtmp,id_notc)
  
}) %>% bind_rows()
rownames(cls_df)<-1:nrow(cls_df)
#cls_df<-vbins_df # no combined methods
#get min Identity for family level



head(df_rcounts)
head(cdf) # all contig info


################################################

###
tax_host_df<-merge(cdf,tax_df%>%filter(!is.na(score)), by=c("sepID","ID"), all.x = T)
tax_host_df<-merge(tax_host_df,df_host, by=c("sepID","ID"), all.x = T)
tax_host_df<-merge(tax_host_df, df_phatyp, by=c("sepID","ID"), all.x=T)
tax_host_df$isphage<-NA

print_stats(tax_host_df)

#extend host by cluster
if(F){
  tax_host_df<-merge(tax_host_df,hbins,  by=c("sepID","ID"), all.x = T)
  tax_host_df<-lapply(unique(tax_host_df$sepID), function(h){
    #print(h)
    tmp<-subset(tax_host_df, sepID==h)
    ntmp<-subset(tmp,is.na(host_bin))
    nntmp<-tmp%>%filter(!is.na(host_bin)) %>% arrange(host_bin)

    etmp<-lapply(nntmp%>%pull(host_bin)%>%unique(), function(o){
      #print(o)
      stmp<-subset(nntmp, host_bin==o)
      uann<-unique(stmp$host_g)
      uann<-uann[!is.na(uann)]
      if(length(uann)==1){
        stmp$host_g<-uann
        uann<-unique(stmp$host_s)
        uann<-uann[!is.na(uann)]
        if(length(uann)==1){
          #print(paste0(h,": ",o))
          stmp$host_g<-stmp%>%drop_na(host_g)%>%pull(host_g)%>%unique()
          stmp$host_s<-uann
        }
      }
      stmp
    }) %>% bind_rows()
    rbind(etmp,ntmp)
  }) %>% bind_rows()
  
}


cdftax<-arrange(tibble(lapply(unique(tax_host_df$sepID),function(h){
  allowed<-df_eval%>%filter(sepID==h)
  stmp<-subset(tax_host_df, sepID==h & ID%in%allowed$ID)
  sbin<-subset(cls_df, sepID==h & ID%in%allowed$ID)
  tmp<-merge(stmp, sbin, by=c("sepID","ID"), all.x=T)
  seq<-max(tmp$bin[!is.na(tmp$bin)])+1:nrow(subset(tmp,is.na(bin)))
  tmp[is.na(tmp$bin),"bin"] <- seq
  tmp
}) %>% bind_rows()),sepID,bin)

##########################

#inherit annot by cluster
cdftax<-lapply(unique(cdftax$sepID), function(h){
  #print(h)
  tmp<-subset(cdftax, sepID==h)
  ntmp<-subset(tmp,is.na(bin))
  etmp<-subset(tmp,!is.na(bin))
  ## family
  etmp<-lapply(unique(etmp$bin), function(o){
    #print(o)
    stmp<-subset(etmp, bin==o)
    uann<-unique(stmp$family)
    uann<-uann[!is.na(uann)]
    uann<-unlist(sapply(uann,function(x)strsplit(x,";;")[[1]]))
    if(length(uann)==1){
      stmp$family<-uann
    }
    if(length(uann)>1){
      allf<-stmp$family[!is.na(stmp$family)]
      allf<-unlist(sapply(allf,function(x)strsplit(x,";;")[[1]]))
      allf<-unlist(sapply(allf,function(x)strsplit(x,";")[[1]]))
      tmpann<-sort(table(allf),decreasing=TRUE)
      stmp$family<-paste(names(tmpann[which(tmpann==max(tmpann))]),collapse = ";")
    }

    stmp
  }) %>% bind_rows()
  
  rbind(etmp,ntmp)
}) %>% bind_rows()

cdftax$isphage<-sapply(1:nrow(cdftax),function(i){
  if(!is.na(cdftax[i,"host_g"])){
    return("yes")
  }
  x<-cdftax[i,"family"]
  if(!is.na(x)){
    if(x %in% pflist){
      return("yes")
    }else{
      return(NA)
    }
  }
  return(NA)
})
cdftax$phagetype<-ifelse(cdftax$isphage=="yes",cdftax$phagetype,NA)
#cdftax$phagetype<-ifelse(!is.na(cdftax$host_g),"temperate",cdftax$phagetype)

#quizas filtrar >3000?
print_stats(cdftax)

#cdftax$length<-NULL

