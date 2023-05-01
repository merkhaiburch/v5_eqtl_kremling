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
rTASSEL::startLogger(fullPath = "/workdir/mbb262/moaseq/mapping/results", fileName = "debug_fast_association.txt")


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
lapply(expression_phenos, function(x) dim(x)-30)


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


# ------------------------------------------------------------------------------
#       Map phenotypes
# ------------------------------------------------------------------------------

# Local variables
vcf <- "/workdir/mbb262/moaseq/genotypes/v5/"

# Set directory for results
setwd("/workdir/mbb262/moaseq/mapping/results/")

# # Run fast association
# lapply(seq_len(10), function(i) {
#   message("I am on chromosome ", i)
#   
#   # Load in genotype table
#   vcf <-  rTASSEL::readGenotypeTableFromPath(path = paste(vcf,"hmp321_282_agpv5_merged_chr", i, "_sorted.vcf.gz", sep = ""),
#                                              keepDepth = FALSE)
#   print(vcf)
#   
#   # Iterate through phenotypes
#   lapply(seq_len(length(tasPhenoDF)), function(j) {
#     message("I am on phenotype file ", names(tasPhenoDF)[j])
#     
#     # Join genotypes with (phenotypes + PCs + PEERs)
#     tasPhenoDFGenotype <- rTASSEL::readGenotypePhenotype(
#       genoPathOrObj = vcf,
#       phenoPathDFOrObj = tasPhenoDF[[j]])
#     print(tasPhenoDFGenotype)
#     
#     # Do a light MAF filter to remove invariant sites
#     tasGenoPhenoFilt <- rTASSEL::filterGenotypeTableSites(
#       tasObj = tasPhenoDFGenotype,
#       siteRangeFilterType = "none",
#       siteMinAlleleFreq = 0.05,
#       siteMaxAlleleFreq = 1.0,
#       siteMinCount = 100)
#     
#     # Run fast association, write files to disk.
#     rTASSEL::assocModelFitter(
#       tasObj = tasGenoPhenoFilt,
#       formula = . ~ ., 
#       fitMarkers = TRUE,
#       kinship = NULL,
#       fastAssociation = TRUE,
#       maxP = p_value,
#       maxThreads = numThreads,
#       outputFile = paste("chrom_", i, "_fast_assoc_results_", names(tasPhenoDF)[j], "_Kremling_2018", sep = ""))
#   })
# })


## Count the number of SNPs post merge -----------------------------------------

# Create a subsetted list
tasPhenoDF_sub <- tasPhenoDF$GRoot

# Run fast association
lapply(seq_len(10), function(i) {
  message("I am on chromosome ", i)
  
  # Load in genotype table
  tas_vcf <-  rTASSEL::readGenotypeTableFromPath(path = paste(vcf,"hmp321_282_agpv5_merged_chr", i, "_sorted.vcf.gz", sep = ""),
                                             keepDepth = FALSE)
  
  # Join genotypes with (phenotypes + PCs + PEERs)
  tasPhenoDF_subGenotype <- rTASSEL::readGenotypePhenotype(
    genoPathOrObj = tas_vcf,
    phenoPathDFOrObj = tasPhenoDF_sub)
  
  # Do a light MAF filter to remove invariant sites
  tasGenoPhenoFilt <- rTASSEL::filterGenotypeTableSites(
    tasObj = tasPhenoDF_subGenotype,
    siteRangeFilterType = "none",
    siteMinAlleleFreq = 0.05,
    siteMaxAlleleFreq = 1.0,
    siteMinCount = 100)
  
  # Return count of SNPs
  print(tasGenoPhenoFilt)
  
})

