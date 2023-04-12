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

# Load in v4 --> v5 conversion table
v42v5 <- read.delim()

# Find all files in a directory
path_2_krem_results <- ""
krem_results <- list.files(path = path_2_krem_results,
                           pattern = "*.txt",
                           full.names = FALSE)

# Create loop that:
# - Loops though each file in a directory
# - Loads that result file in with fread
# - merges result file with the v42v5 conversion table
# - rearrages the merged results
# - exports results to file

# Function that does the above things
for(i in 1:length(krem_results)){
  # Load in an individual result
  result_df <- data.table::fread(paste0(path_2_krem_results, krem_results[i]))
  
  # Merge result file with the v42v5 conversion table
  merge_result_df <- merge(v42v5, result_df, by.x = "v4", by.y = "Trait", all = FALSE)
  
  # Rearrange merged results
  merge_result_df <- merge_result_df %>% 
    select("v4", "Marker", "Chr", "Pos", "df", "r2", "p")
  
  # Change column name back to Trait (is now in v5)
  colnames(merge_result_df)[1] <- "Trait"
  
  # Export result back to file
  data.table::fwrite(x = paste0("v5_", krem_results), file = "path/to/object/out")
}