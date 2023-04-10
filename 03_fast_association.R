# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-03-24
# Updated... 2023-03-27
#
# Description:
# Preform eQTL mapping using uplifted (v4-->v5) Hapmap 321 SNPs with the 
# expression data from Kremling et al., 2018
# ------------------------------------------------------------------------------

# Set directory, load packages, format -----------------------------------------

# Run garbage collection
gc(full = TRUE, verbose = TRUE)

# Setting memory and global parameters
options(java.parameters = c("-Xmx500g"))
numThreads <- 63
p_value <- 0.00001

setwd("/workdir/mbb262/moaseq/mapping/results/")

# Load in packages
library(dplyr)
library(magrittr)
library(data.table)
library(rTASSEL)

# Start logging
rTASSEL::startLogger(fullPath = NULL, fileName = NULL)


# ------------------------------------------------------------------------------
# Load PCs and PEERs
# ------------------------------------------------------------------------------

# Get all file names
peer_dir = list.files("/workdir/mbb262/moaseq/mapping/peers", pattern="*.txt", full.names = TRUE)

# Load pcs + peers into r as sepearte dfs within a single list
pcs_peers <- lapply(peer_dir,
                    data.table::fread, 
                    header=TRUE,
                    nThread = numThreads)

# Shorten names on pcs + peers
peer_short_names <- gsub("/workdir/mbb262/moaseq/mapping/peers/TASSEL_HEADER_df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_", "", peer_dir)
peer_short_names <- gsub("_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs_avgDups.txt", "", peer_short_names)
names(pcs_peers) <- peer_short_names # rename each element in the list to the tissue

# Check pcs + peers
lapply(pcs_peers, function(x) x[1:5,1:8])
lapply(pcs_peers, function(x) dim(x))

# Remove phenotype columns from karls data
pcs_peers <- lapply(pcs_peers, function(x) x[,1:31])
lapply(pcs_peers, function(x) dim(x))

# ------------------------------------------------------------------------------
# Load in Kremling expression phenotypes (in v3)
# ------------------------------------------------------------------------------

# Get all file names
phenotype_dir <- list.files("/workdir/mbb262/moaseq/mapping/phenotypes", pattern="*.csv", full.names = TRUE)

# Load data into r as sepearte dfs within a single list
expression_phenos <- lapply(phenotype_dir,
                     data.table::fread, 
                     header=TRUE,
                     nThread = numThreads)

# Shorten names
short_names <- gsub("_kremling_formatted_v4_hapmapids.csv", "", phenotype_dir)
short_names <- gsub("/workdir/mbb262/moaseq/mapping/phenotypes/", "", short_names)
names(expression_phenos) <- short_names # rename each element in the list to the tissue

# Check
lapply(expression_phenos, function(x) x[1:5,1:5])
lapply(expression_phenos, function(x) dim(x))


# ------------------------------------------------------------------------------
# Join expression pehnotypes with PCs + PEERs
# ------------------------------------------------------------------------------

# Check order of both lists, must be the same
names(expression_phenos)
names(pcs_peers)

# Innerjoin two lists by index and by taxa name
expression_pcs_peers <- purrr::map2(pcs_peers, expression_phenos, inner_join, by = "Taxa")

# Check
str(expression_pcs_peers)
names(expression_pcs_peers)
expression_pcs_peers$GRoot[1:5,1:33]
lapply(expression_pcs_peers, function(x) dim(x))

# Tasselize merged phenotypes + PCs + PEERS
tasPhenoDF <- lapply(expression_pcs_peers, function(x) {rTASSEL::readPhenotypeFromDataFrame(
  phenotypeDF = x,
  taxaID = "Taxa",
  attributeTypes = c(rep("covariate", 30), rep("data", ncol(x)-31)))
})

# Check
tasPhenoDF
names(tasPhenoDF[1])


# ------------------------------------------------------------------------------
#       Map phenotypes
# ------------------------------------------------------------------------------

# Local variables
vcf <- "/workdir/mbb262/moaseq/genotypes/"

# Set directory for results
setwd("/workdir/mbb262/moaseq/mapping/results/")

# Run fast association
lapply(seq_len(10), function(i) {
  message("I am on chromosome ", i)
  
  # Load in genotype table
  vcf <-  rTASSEL::readGenotypeTableFromPath(path = paste(vcf,"hmp321_282_agpv4_merged_chr", i, "_imputed_goodman282.vcf.gz", sep = ""),
                                             keepDepth = FALSE)
  print(vcf)
  
  # Iterate through phenotypes
  lapply(seq_len(length(tasPhenoDF)), function(j) {
    message("I am on phenotype file ", names(tasPhenoDF)[j])
    
    # Join genotypes with (phenotypes + PCs + PEERs)
    tasPhenoDFGenotype <- rTASSEL::readGenotypePhenotype(
      genoPathOrObj = vcf,
      phenoPathDFOrObj = tasPhenoDF[[j]])
    print(tasPhenoDFGenotype)
    
    # Do a light MAF filter to remove invariant sites
    tasGenoPhenoFilt <- rTASSEL::filterGenotypeTableSites(
      tasObj = tasPhenoDFGenotype,
      siteRangeFilterType = "none",
      siteMinAlleleFreq = 0.05,
      siteMaxAlleleFreq = 1.0,
      siteMinCount = 100)
    
    # Run fast association, write files to disk.
    rTASSEL::assocModelFitter(
      tasObj = tasGenoPhenoFilt,
      formula = . ~ ., 
      fitMarkers = TRUE,
      kinship = NULL,
      fastAssociation = TRUE,
      maxP = p_value,
      maxThreads = numThreads,
      outputFile = paste("chrom_", i, "_fast_assoc_results_", names(tasPhenoDF)[j], "_Kremling_2018", sep = ""))
  })
})


## Go from v3 to v5 gene ids ---------------------------------------------------

# Load in cross reference table, rename columns
v4_v5 <- read.delim("/workdir/mbb262/moaseq/mapping/B73v4_B73v5_liftoff_genemodel_CDS_xref_shortened.txt", header = TRUE)
colnames(v4_v5) <- c("v4", "v5")

# Remove any genes that did not map
v4_v5 <- v4_v5 %>% filter(!grepl('chr', v5)) 
v4_v5 <- v4_v5 %>% filter(!grepl('scaf', v5))

# only get unique v4 ids
v4_v5 <- v4_v5 %>% distinct(v4, .keep_all = TRUE)

# Change names by:
# 1) Loop though each file in the results directory, load 1 in at a time
# 2) merge v3 gene id names with v5 names
# 3) drop v3 names, rearrange columns so v5 names are first
# 4) export in new directory

# Find all files to iterate through
results_dir <- "/workdir/mbb262/moaseq/mapping/results/"
fast_assoc_results <- list.files(results_dir, pattern="*.txt")

# Iterate through files





