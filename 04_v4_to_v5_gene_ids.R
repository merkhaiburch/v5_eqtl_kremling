# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-04-11
# Updated... 2023-04-11
#
# Description:
# Go from v4 gene ids to v5
# ------------------------------------------------------------------------------

# Load in packages
library(dplyr)
library(data.table)

# Set variables
thread_n <- 60

# Set options
options(scipen = 999)


## Gather data -----------------------------------------------------------------

# Load in v4 --> v5 conversion table
v4_v5 <- data.table::fread("/workdir/mbb262/moaseq/mapping/B73v4_B73v5_liftoff_genemodel_CDS_xref_shortened.txt", header = TRUE)
colnames(v4_v5) <- c("v4", "v5")

# Remove any genes that did not map
v4_v5 <- v4_v5 %>% filter(!grepl('chr', v5)) 
v4_v5 <- v4_v5 %>% filter(!grepl('scaf', v5))

# only get unique v4 ids
v4_v5 <- v4_v5 %>% distinct(v4, .keep_all = TRUE)

# Find all result files in a directory
path_2_krem_results <- "/workdir/mbb262/moaseq/mapping/results/"
krem_results <- list.files(path = path_2_krem_results,
                           pattern = "*.txt",
                           full.names = FALSE)


## Function to go from v4 --> v5 gene Ids --------------------------------------

# Create loop that:
# - Loops though each file in a directory
# - Loads that result file in with fread
# - merges result file with the v42v5 conversion table
# - rearrages the merged results
# - exports results to file

# Function that does the above things
for(i in 1:length(krem_results)){
  # Load in an individual result
  result_df <- data.table::fread(paste0(path_2_krem_results, krem_results[i]), nThread = thread_n)
  message(paste0(krem_results[i]))

  # Merge result file with the v4_v5 conversion table
  setkey(v4_v5,v4)
  setkey(result_df,Trait)
  merge_result_df <- v4_v5[result_df, nomatch=0]
  
  # Rearrange merged results
  merge_result_df <- merge_result_df %>% select(-"v4")
  
  # Change column name back to Trait (is now in v5)
  colnames(merge_result_df)[1] <- "Trait"
  
  # Export result back to file
  data.table::fwrite(x = merge_result_df, 
                     file = paste0("/workdir/mbb262/moaseq/mapping/results/v4_2_v5/v5_", krem_results[i]),
                     nThread = thread_n)
}









