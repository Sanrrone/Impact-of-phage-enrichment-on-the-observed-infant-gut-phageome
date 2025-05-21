setwd("~/Desktop/phage_enrichment_newassembly/")
source("functions/parse_phagcn.R")
#source("functions/parse_kr2.R")
source("functions/parse_virustaxo.R")
source("functions/parse_blast.R")

#source("functions/parse_AV.R")
#tag phage families
pf<-uhgv_metadata %>% 
  filter(isphage=="yes") %>% 
  drop_na(family_name) %>%
  unique()
pflist<-pf%>%pull(family_name)%>%unique()
pflist<-c(pflist,"Aliceevansviridae","Fiersviridae","Ahpuchviridae","Forsetiviridae", "Intestiviridae","Plectroviridae","Rudiviridae",
          "Steigviridae","Suoliviridae")

#head(df_kr2)
head(df_phagcn)
head(df_vtaxo)
#head(df_av)
head(df_blast)

#tmp<-merge(dfphagetag, df_isphage, by=c("sepID","ID"), all=T)
#tmp$isphage<-ifelse(tmp$isphage.x=="yes" | tmp$isphage.y=="yes", "yes", NA)
#tmp$isphage.x<-NULL
#tmp$isphage.y<-NULL
#df_isphage<-tmp


tax_df<-lapply(unique(cdf$sepID), function(h){
  ptmp<-subset(df_phagcn %>% filter(family%in%pflist), sepID==h)
  colnames(ptmp)[which(colnames(ptmp)=="family")]<-"ph_family"
  colnames(ptmp)[which(colnames(ptmp)=="genus")]<-"ph_genus"
  colnames(ptmp)[which(colnames(ptmp)=="species")]<-"ph_species"
  
  #ktmp<-subset(df_kr2 %>% filter(family%in%pflist), sepID==h)
  #colnames(ktmp)[which(colnames(ktmp)=="family")]<-"k2_family"
  #colnames(ktmp)[which(colnames(ktmp)=="genus")]<-"k2_genus"
  #colnames(ktmp)[which(colnames(ktmp)=="species")]<-"k2_species"
  
  vtmp<-subset(df_vtaxo %>% filter(family%in%pflist), sepID==h)
  colnames(vtmp)[which(colnames(vtmp)=="family")]<-"vt_family"
  colnames(vtmp)[which(colnames(vtmp)=="genus")]<-"vt_genus"
  colnames(vtmp)[which(colnames(vtmp)=="species")]<-"vt_species"
  
  btmp<-subset(df_blast %>% filter(family%in%pflist), sepID==h)
  colnames(btmp)[which(colnames(btmp)=="family")]<-"bl_family"
  colnames(btmp)[which(colnames(btmp)=="genus")]<-"bl_genus"
  colnames(btmp)[which(colnames(btmp)=="species")]<-"bl_species"
  
  alldf<-Reduce(function(x,y) merge(x,y,by=c("sepID","ID"), all=T), 
                list(ptmp,btmp,vtmp))#,vatmp))
  if(nrow(alldf)==0){return(NULL)}
  alldf$family<-sapply(1:nrow(alldf),function(i){
    tmp<-alldf[i,]
    tmpf<-unique(c(tmp$ph_family,tmp$bl_family,tmp$vt_family))#,tmp$vc_family))
    tmpf<-tmpf[!is.na(tmpf)]
    #tmp<-paste(tmpf[!is.na(tmpf)],collapse = ";;")
    if(length(tmpf)>1){
      return(NA)
    }else{
      if(is.na(tmpf)){return(NA)}
      if(tmpf==""){return(NA)}
    }
    tmpf
  })
  alldf$genus<-sapply(1:nrow(alldf),function(i){
    tmp<-alldf[i,]
    tmpf<-unique(c(tmp$ph_genus,tmp$bl_genus,tmp$vt_genus))#,tmp$vc_family))
    tmpf<-tmpf[!is.na(tmpf)]
    #tmp<-paste(tmpf[!is.na(tmpf)],collapse = ";;")
    if(length(tmpf)==0){return(NA)}
    if(length(tmpf)>1){
      return(NA)
    }else{
      if(is.na(tmpf)){return(NA)}
      if(tmpf==""){return(NA)}
    }
    tmpf
  })
  

  alldf$score<-sapply(1:nrow(alldf),function(i){
    tmp<-alldf[i,]
    if(is.na(tmp[,"family"])){return(NA)}
    pv<-round(tmp$pvalue,2)
    bv<-round(tmp$idcov,2)
    vv<-round(tmp$negentrophy,2)
    #sv<-round(tmp$vc_score,2)
    paste(collapse = ";;",c(paste0("pg",pv),paste0("bl",bv),paste0("vt",vv)))#,paste0("vc",sv)))
  })
  
  alldf[,c("sepID","ID","family","genus","score")]
  
}) %>% bind_rows

gc()

