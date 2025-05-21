setwd("~/Desktop/phage_enrichment_newassembly")
l<-list.files(path = "bins_aai/",pattern = "_family_clusters.txt", full.names = T)
aaibins_df<-lapply(l,function(clstr_file){
  sname<-gsub(pattern = "_family_clusters.txt", replacement = "", x=basename(clstr_file))
  
  clstr_content <- readLines(clstr_file)
  
  # Initialize variables
  current_cluster_id <- NULL
  sequence_ids <- list()
  
  # Process each line in the .clstr file
  cn<-1
  df<-data.frame(ID=NULL,bin=NULL)
  for (line in clstr_content) {
    contig_names <- strsplit(line, "\t")[[1]]
    
    # Create a dataframe
    df <- rbind(df,data.frame(ID=contig_names, bin = cn, stringsAsFactors = FALSE))
    cn<-cn+1
  }
  df$hepID<-sname
  df
  
}) %>% bind_rows()

l<-list.files(path = "bins_aai/",pattern = "_genus_clusters.txt", full.names = T)
g_aaibins_df<-lapply(l,function(clstr_file){
  sname<-gsub(pattern = "_genus_clusters.txt", replacement = "", x=basename(clstr_file))
  
  clstr_content <- readLines(clstr_file)
  
  # Initialize variables
  current_cluster_id <- NULL
  sequence_ids <- list()
  
  # Process each line in the .clstr file
  cn<-1
  df<-data.frame(ID=NULL,bin=NULL)
  for (line in clstr_content) {
    contig_names <- strsplit(line, "\t")[[1]]
    
    # Create a dataframe
    df <- rbind(df,data.frame(ID=contig_names, bin = cn, stringsAsFactors = FALSE))
    cn<-cn+1
  }
  df$hepID<-sname
  df
  
}) %>% bind_rows()




