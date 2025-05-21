setwd("~/Desktop/phage_enrichment_newassembly")
l<-list.files(path = "bins_cdhit/",pattern = "_cdhit95.clstr", full.names = T)

cdbins_df<-lapply(l,function(clstr_file){
# Read the content of the .clstr file
  sname<-gsub(pattern = "_cdhit95.clstr", replacement = "", x=basename(clstr_file))
  
clstr_content <- readLines(clstr_file)

# Initialize variables
current_cluster_id <- NULL
sequence_ids <- list()

# Process each line in the .clstr file
for (line in clstr_content) {
  if (startsWith(line, ">Cluster")) {
    # Extract the cluster ID from the line
    current_cluster_id <- as.numeric(strsplit(line, " ")[[1]][2])
  } else {
    # Extract the sequence ID from the line
    seq_id <- strsplit(line,"\\.\\.\\.")[[1]][1]
    seq_id<-strsplit(seq_id,">")[[1]][2]
    
    # Store the sequence ID and cluster ID
    sequence_ids[[seq_id]] <- current_cluster_id
  }
}

# Create a dataframe from the parsed data
clstr_df <- data.frame(
  ID = names(sequence_ids),
  bin = unlist(unname(sequence_ids)),
  stringsAsFactors = FALSE
)
clstr_df$bin<-clstr_df$bin+1
clstr_df$hepID<-sname
clstr_df
}) %>% bind_rows()
cdbins_df
