#ICTV
ictv_df<-tibble(read.table("ICTV_all.tsv", header = T, sep = "\t", quote = "", comment.char = "", fill=T))
ictv_df<-ictv_df[,c("Sort","Kingdom","Phylum","Class","Order","Family","Genus","Species","Virus.GENBANK.accession", "Virus.REFSEQ.accession")]

#species
l<-list.files(path = "clusters/",pattern = "_s.clstr", full.names = T)
sclstr<-lapply(l, function(tsv){
  #print(tsv)
  sname<-strsplit(basename(tsv),"_s.clstr")[[1]][1]
  tmp<-read.table(tsv, header = F, sep = "\t", skip = 1)
  colnames(tmp)<-c("s__OTU","ID")
  tmp<-lapply(unique(tmp$s__OTU), function(x){
    stmp<-subset(tmp, s__OTU==x)
    stmp$s__center<-c("C",rep("M", nrow(stmp)-1))
    stmp
  }) %>% bind_rows
  tmp$sepID<-sname
  tmp
}) %>% bind_rows

#family
l<-list.files(path = "clusters/",pattern = "_f.clstr", full.names = T)
fclstr<-lapply(l, function(tsv){
  sname<-strsplit(basename(tsv),"_f.clstr")[[1]][1]
  tmp<-read.table(tsv, header = F, sep = "\t", skip = 1)
  colnames(tmp)<-c("f__OTU","ID")
  tmp<-lapply(unique(tmp$f__OTU), function(x){
    stmp<-subset(tmp, f__OTU==x)
    stmp$f__center<-c("C",rep("M", nrow(stmp)-1))
    stmp
  }) %>% bind_rows
  tmp$sepID<-sname
  tmp
}) %>% bind_rows()

clstr_df<-merge(sclstr, fclstr, by=c("sepID","ID"), all=T)
rm(fclstr, sclstr)
