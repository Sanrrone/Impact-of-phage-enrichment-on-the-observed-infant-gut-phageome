setwd("~/Desktop/phage_enrichment_newassembly")

uhgv_metadata<-read.table("/data/tmp/UHGV/votus_metadata_extended.tsv", sep = "\t", quote = "", header = T, na.strings = "NULL")
uhgv_metadata$kmers_host_lineage_gtdb_r207<-ifelse(uhgv_metadata$kmers_host_agreement>80,
                                                   uhgv_metadata$kmers_host_lineage_gtdb_r207,NA)
uhgv_metadata$host_lineage_gtdb_r207<-ifelse(is.na(uhgv_metadata$crispr_host_lineage_gtdb_r207)&
                                               is.na(uhgv_metadata$kmers_host_lineage_gtdb_r207),
                                             NA,uhgv_metadata$host_lineage_gtdb_r207)

uhgv_metadata<-uhgv_metadata%>%filter(checkv_quality %in% c("Complete","High-quality","Medium-quality"))                                             
uhgv_metadata<-uhgv_metadata[,c("uhgv_genome","uhgv_taxonomy","uhgv_votu","checkv_quality",
                                "family_name", "class_name", "genus_name","host_lineage_gtdb_r207",
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
                                "family_name","genus_name","uhgv_taxonomy","isphage")]




l<-list.files(path = "taxprof_contig_blast", pattern = "_blast.tsv", full.names = T)
df_blast<-lapply(l, function(tsv){
  #print(tsv)
  sname<-gsub("_blast.tsv","",basename(tsv))
  if(countLines(tsv)<=1){return(NULL)}
  tmp<-read.table(tsv,header = F, sep = "\t",skipNul = T, stringsAsFactors = F)
  colnames(tmp)<-c("ID","uhgv_genome","identity","qcov")
  tmp<-subset(tmp, identity>0.7 & qcov>0.7)
  if(nrow(tmp)==0){
    return(NULL)
  }
  tmp$idcov<-(tmp$identity*tmp$qcov)/10000
  tmp<-merge(tmp,uhgv_metadata, by="uhgv_genome",all.x=T)
  colnames(tmp)[which(colnames(tmp)=="family_name")]<-"family"
  colnames(tmp)[which(colnames(tmp)=="genus_name")]<-"genus"
  tmp$species<-NA
  tmp$sepID<-sname
  tmp
}) %>% bind_rows()
