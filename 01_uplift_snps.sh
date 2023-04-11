#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-03-23
# Updated... 2023-04-10
#
# Description:
# Uplift Hapmap 3.2.1 SNPs from Ranstein et al 2021 to v5
# https://doi.org/10.1534/genetics.120.303025
# ------------------------------------------------------------------------------

# Create variables
HOME_DIR_GENO=/workdir/mbb262/moaseq/genotypes
V4_DIR=$HOME_DIR_GENO/v4
V5_DIR=$HOME_DIR_GENO/v5

# Make directories
mkdir -p $HOME_DIR_GENO
mkdir -p $V4_DIR
mkdir -p $V5_DIR


## Gather data -----------------------------------------------------------------

# Get SNPs from cbsu
scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/gr226/Hmp321/imputed/AGPv4/*vcf.gz $V4_DIR

# Get cross-map file from maizegdb and v5 genome
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/chain_files/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain -P $HOME_DIR_GENO
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz -P $HOME_DIR_GENO

# Unzip genome
gunzip $HOME_DIR_GENO/Zm-B73-REFERENCE-NAM-5.0.fa.gz


## Use crossmap to uplift files ------------------------------------------------
# Set path
export PYTHONPATH=/programs/CrossMap-0.6.1/lib64/python3.9/site-packages:/programs/CrossMap-0.6.1/lib/python3.9/site-packages
export PATH=/programs/CrossMap-0.6.1/bin:$PATH

# Iterate through each chromosome & uplifts 282/GAP coordinates
cd $V5_DIR
for CHROM in {1..10}
do
  echo "I am on chr ${CHROM}"

  CrossMap.py vcf \
    $HOME_DIR_GENO/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain \
    $V4_DIR/hmp321_282_agpv4_merged_chr${CHROM}.imputed.vcf.gz \
    $HOME_DIR_GENO/Zm-B73-REFERENCE-NAM-5.0.fa \
    $V5_DIR/hmp321_282_agpv5_merged_chr${CHROM}.vcf

  bgzip --threads 20 $V5_DIR/hmp321_282_agpv5_merged_chr${CHROM}.vcf

  echo "I just finished chr ${CHROM}"
done


## Sort genotype file ----------------------------------------------------------
for CHROM in {1..10}
do
  echo "I am on chr ${CHROM}"

  /home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl \
    -debug $V5_DIR/debug2.txt \
    -Xmx420g \
    -maxThreads 60 \
    -SortGenotypeFilePlugin \
    -inputFile $V5_DIR/hmp321_282_agpv5_merged_chr${CHROM}.vcf.gz \
    -outputFile $V5_DIR/hmp321_282_agpv5_merged_chr${CHROM}_sorted.vcf.gz

  echo "I just finished chr ${CHROM}"
done

# Alternatively
bgzip -c hmp321_282_agpv5_merged_chr10.vcf > hmp321_282_agpv5_merged_chr10.vcf.gz
tabix -p vcf hmp321_282_agpv5_merged_chr10.vcf.gz
/programs/bcftools-1.15.1-r/bin/bcftools sort -Oz hmp321_282_agpv5_merged_chr10.vcf.gz -o hmp321_282_agpv5_merged_chr10_sorted_bcf.vcf.gz
