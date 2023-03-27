# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-03-24
# Updated... 2023-03-24
#
# Description:
# Prep to preform eQTL mapping using
# ------------------------------------------------------------------------------

# Make directories on cbsu -------
mkdir -p /workdir/mbb262/mapping

HOME_DIR=/workdir/mbb262/mapping

mkdir $HOME_DIR/peers
mkdir $HOME_DIR/phenotypes
mkdir $HOME_DIR/results
mkdir $HOME_DIR/genotypes

# Download data ------------------

# Get PEERS from my blfs1 account --> do not have to remove first lines of file
scp mbb262@cbsublfs1.tc.cornell.edu:/data1/users/mbb262/phenotypes/1_individual_datasets/kremling_2018_naturegenetics_expression/peers/*.txt /workdir/mbb262/peers_copy

# Code used to make these files in phenotype merging file, these are the cleaned up files on blfs1
scp mbb262@cbsublfs1.tc.cornell.edu:/data1/users/mbb262/phenotypes/1_individual_datasets/kremling_2018_naturegenetics_expression/normalized_counts/kremling_raw_count_v4_hapmap321taxaid/*.csv /workdir/mbb262/phenotypes/kremling/

# scp mbb262@cbsublfs1.tc.cornell.edu:/data1/users/mbb262/genotypes/goodman282/*_imputed_goodman282.vcf.gz /workdir/mbb262/genotypes/goodman
