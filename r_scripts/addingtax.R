setwd("~/Desktop/phage_enrichment_newassembly/")

source("functions/functions.R")
#write.table(df_rcounts, "dfcounts.tsv", sep="\t",quote = F,row.names = F, col.names = T)
#cdftax<-read.table("cdftax.tsv",sep = "\t",header = T)

source("functions/annot_and_identity.R")
#write.table(cdftax,"cdftax.tsv",quote = F,row.names = F,col.names = T,sep = "\t")
source("functions/parse_propagate.R")

m_df<-merge(df_rcounts, cdftax[,-which("length"==colnames(cdftax))], by=c("hepID","ID"), all.x = T)
m_df<-tibble(merge(m_df,df_prop, by=c("sample","ID"), all.x = T))

write.table(m_df,"alldata.tsv",sep = "\t", quote = F, row.names = F, col.names = T)


