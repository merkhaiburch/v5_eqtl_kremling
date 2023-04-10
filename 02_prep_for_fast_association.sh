# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-03-24
# Updated... 2023-03-27
#
# Description:
# Prep to preform eQTL mapping using
# ------------------------------------------------------------------------------

# Make directories on cbsu -------
mkdir -p /workdir/mbb262/moaseq/mapping

HOME_DIR=/workdir/mbb262/moaseq/mapping

mkdir $HOME_DIR/peers
mkdir $HOME_DIR/phenotypes
mkdir $HOME_DIR/results

## Download data ---------------------------------------------------------------

# Download in PEERs from Cyverse
wget https://de.cyverse.org/anon-files//iplant/home/shared/commons_repo/curated/Kremling_maizeDiversityPanelRNAseq_2018/Expression_matrix/TASSEL_fmt_expression_w_covariates_corrected_header/TASSEL_HEADER_df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_GRoot_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs_avgDups.txt.gz --directory $HOME_DIR/peers
wget https://de.cyverse.org/anon-files//iplant/home/shared/commons_repo/curated/Kremling_maizeDiversityPanelRNAseq_2018/Expression_matrix/TASSEL_fmt_expression_w_covariates_corrected_header/TASSEL_HEADER_df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_Kern_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs_avgDups.txt.gz --directory $HOME_DIR/peers
wget https://de.cyverse.org/anon-files//iplant/home/shared/commons_repo/curated/Kremling_maizeDiversityPanelRNAseq_2018/Expression_matrix/TASSEL_fmt_expression_w_covariates_corrected_header/TASSEL_HEADER_df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_L3Base_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs_avgDups.txt.gz --directory $HOME_DIR/peers
wget https://de.cyverse.org/anon-files//iplant/home/shared/commons_repo/curated/Kremling_maizeDiversityPanelRNAseq_2018/Expression_matrix/TASSEL_fmt_expression_w_covariates_corrected_header/TASSEL_HEADER_df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_L3Tip_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs_avgDups.txt.gz --directory $HOME_DIR/peers
wget https://de.cyverse.org/anon-files//iplant/home/shared/commons_repo/curated/Kremling_maizeDiversityPanelRNAseq_2018/Expression_matrix/TASSEL_fmt_expression_w_covariates_corrected_header/TASSEL_HEADER_df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_LMAD_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs_avgDups.txt.gz --directory $HOME_DIR/peers
wget https://de.cyverse.org/anon-files//iplant/home/shared/commons_repo/curated/Kremling_maizeDiversityPanelRNAseq_2018/Expression_matrix/TASSEL_fmt_expression_w_covariates_corrected_header/TASSEL_HEADER_df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_LMAN_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs_avgDups.txt.gz --directory $HOME_DIR/peers
wget https://de.cyverse.org/anon-files//iplant/home/shared/commons_repo/curated/Kremling_maizeDiversityPanelRNAseq_2018/Expression_matrix/TASSEL_fmt_expression_w_covariates_corrected_header/TASSEL_HEADER_df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_GShoot_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs_avgDups.txt.gz --directory $HOME_DIR/peers
gunzip $HOME_DIR/peers/*

# Remove first two lines with <Phenotype> from tassel header
for FILE in $HOME_DIR/peers/TASSEL_H*
do
  sed -i '1d' $FILE # ran twice because I'm not very clever
  sed -i '1d' $FILE
done

# Grab phenotypes
# Code used to make these files in phenotype merging file, these are the cleaned up files on blfs1
scp mbb262@cbsublfs1.tc.cornell.edu:/data1/users/mbb262/phenotypes/1_individual_datasets/kremling_2018_naturegenetics_expression/normalized_counts/kremling_raw_count_v4_hapmap321taxaid/*.csv $HOME_DIR/phenotypes

# remove L3
rm $HOME_DIR/phenotypes/L3_kremling_formatted_v4_hapmapids.csv

# Download the v4 to v5 cross reference file
wget https://de.cyverse.org/anon-files//iplant/home/maizegdb/maizegdb/B73v5_JBROWSE_AND_ANALYSES/B73v3-B73v5_and_B73v4-B73v5_gene_model_associations/B73v4_B73v5_liftoff_genemodel_CDS_xref_shortened.txt.gz -P $HOME_DIR
gunzip $HOME_DIR/B73v4_B73v5_liftoff_genemodel_CDS_xref_shortened.txt.gz

# Remove first line of cross reference file
sed -i '1d' $HOME_DIR/B73v4_B73v5_liftoff_genemodel_CDS_xref_shortened.txt


