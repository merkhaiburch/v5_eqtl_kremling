# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-03-24
# Updated... 2023-03-24
#
# Description:
# Preform eQTL mapping using uplifted (v4-->v5) Hapmap 321 SNPs with the 
# expression data from Kremling et al., 2018
# ------------------------------------------------------------------------------

# Make directories on cbsu -------
# mkdir /workdir/mbb262/
# mkdir /workdir/mbb262/peers_copy
# mkdir /workdir/mbb262/phenotypes/
# mkdir /workdir/mbb262/phenotypes/kremling
# mkdir /workdir/mbb262/results
# mkdir /workdir/mbb262/results/kremling_all
# mkdir /workdir/mbb262/genotypes
# mkdir /workdir/mbb262/genotypes/goodman

# Download data ------------------

# Get PEERS from my blfs1 account --> do not have to remove first lines of file
# scp mbb262@cbsublfs1.tc.cornell.edu:/data1/users/mbb262/phenotypes/1_individual_datasets/kremling_2018_naturegenetics_expression/peers/*.txt /workdir/mbb262/peers_copy

# Code used to make these files in phenotype merging file, these are the cleaned up files on blfs1
# scp mbb262@cbsublfs1.tc.cornell.edu:/data1/users/mbb262/phenotypes/1_individual_datasets/kremling_2018_naturegenetics_expression/normalized_counts/kremling_raw_count_v4_hapmap321taxaid/*.csv /workdir/mbb262/phenotypes/kremling/

# mkdir /workdir/mbb262/genotypes/
# scp mbb262@cbsublfs1.tc.cornell.edu:/data1/users/mbb262/genotypes/goodman282/*_imputed_goodman282.vcf.gz /workdir/mbb262/genotypes/goodman


# Run garbage collection
gc(full = TRUE, verbose = TRUE)

# Setting memory and global parameters
options(java.parameters = c("-Xmx500g"))
numThreads <- 63
p_value <- 0.00001


# Set directory ------------------
setwd("/workdir/mbb262/results/kremling_all/")

# Load in packages
library(dplyr)
library(magrittr)
library(data.table)
library(rTASSEL)

# Start logging
rTASSEL::startLogger(fullPath = NULL, fileName = NULL)


# ----------------------------------
#   Load in pre-analyzed Covariates
# ----------------------------------

# Load in global PCs (Ames to Goodman method)
ames2goodman_gPCs_3gPCs <- read.csv("/home/mbb262/git_projects/pleiotropy/data/Q_all_3gPCs_allGoodman.csv", header = TRUE)
ames2goodman_gPCs_3gPCs <- ames2goodman_gPCs_3gPCs[,-c(5:6)]
colnames(ames2goodman_gPCs_3gPCs) <- c("Taxa", "PC1", "PC2", "PC3")

# Download in PEERs from irods (do above methotd instead)
# icd /iplant/home/shared/commons_repo/curated/Kremling_Nature3RNASeq282_March2018/Expression_matrix/TASSEL_fmt_expression_w_covariates
# iget -r ./
# Remove first two lines with <Phenotype> from tassel header (I would load these into R with Rtassel but I keep getting an error)
# sed -i '1d' ./*.txt # ran twice because I'm not very clever
# sed -i '1d' ./*.txt

# Load in peers
setwd("/workdir/mbb262/peers_copy/")
temp = list.files(pattern="*.txt")
myfiles = lapply(temp, 
                 data.table::fread, 
                 header=TRUE, 
                 sep = "\t", 
                 select = c(1, 7:31), 
                 nThread = numThreads)

# Create shorter name for lists and then rename myfiles
short_names <- gsub("TASSEL_HEADER_df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_", "", temp)# Merge the files together
short_names <- gsub("_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs_avgDups.txt", "", short_names)
names(myfiles) <- short_names

# Join PEERS with my PCs
myfiles <- lapply(myfiles, inner_join, x = ames2goodman_gPCs_3gPCs)


# ------------------------
# Format Kremling data
# ------------------------

# Get all file names
setwd("/workdir/mbb262/phenotypes/kremling/")
temp = list.files("/workdir/mbb262/phenotypes/kremling/", pattern="*.csv")

# Load data into r as sepearte dfs within a single list
formatted_kremling <- lapply(temp, 
                             data.table::fread, 
                             header=TRUE, 
                             nThread = numThreads)

# Shorten names
short_names <- gsub("_kremling_formatted_v4_hapmapids.csv", "", temp)
names(formatted_kremling) <- short_names # rename each element in the list to the tissue

# Remove L3 from the list
# formatted_kremling <- formatted_kremling[-4]

# Add on tissue name to gene names
for (i in seq_along(formatted_kremling)){
  colnames(formatted_kremling[[i]]) <- paste0(colnames(formatted_kremling[[i]]), "_", names(formatted_kremling[i]))
  colnames(formatted_kremling[[i]])[1] <- "Taxa"
}

# Check
lapply(formatted_kremling, function(x) x[1:5,1:5])


# -------------------------------------
# Join expression data with PCs + PEERs
# -------------------------------------

# Check order of both lists, must be the same
names(formatted_kremling)
names(myfiles)

# Innerjoin two lists by index and by taxa name
expression_pcs_peers <- purrr::map2(myfiles, formatted_kremling, inner_join, by = "Taxa")

# Check
# str(expression_pcs_peers)
# names(expression_pcs_peers)
# expression_pcs_peers$GRoot[1:5,1:30]
# expression_pcs_peers$LMAN[1:5,1:30]
# str(expression_pcs_peers$LMAN[1:5,1:30])

# Tasselize merged phenotypes + PCs + PEERS
tasPhenoDF <- lapply(expression_pcs_peers, function(x) {rTASSEL::readPhenotypeFromDataFrame(
  phenotypeDF = x,
  taxaID = "Taxa",
  attributeTypes = c(rep("covariate", 28), rep("data", ncol(x)-29)))
})
tasPhenoDF
names(tasPhenoDF[1])

# Subset list to re-run problematic files (L3Tip, L3Base, GShoot, Kern)
# tasPhenoDF <- list(tasPhenoDF$GShoot, tasPhenoDF$Kern, tasPhenoDF$L3Base, tasPhenoDF$L3Tip)
# names(tasPhenoDF) <- c("GShoot", "Kern", "L3Base", "L3Tip")


# ---------------------------
#       Map phenotypes
# ---------------------------

# Local variables
vcf <- "/workdir/mbb262/genotypes/goodman/"

# Set directory for results
setwd("/workdir/mbb262/results/kremling_all/")

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

# Save on blfs1
# scp /workdir/mbb262/results/kremling_all/* mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/results/pleiotropy/gwa_results_goodman_panel_kremling_maffilter/filtered_goodman
# scp /workdir/mbb262/results/kremling_all/*gz mbb262@cbsumm32.biohpc.cornell.edu:/workdir/mbb262/results/kremling_filtered


