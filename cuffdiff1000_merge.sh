#!/bin/bash

## merge results with perl scripts, remove runs

## modify:
output_dir="/home/prak27/cuffdiff_test"


## input
## get species and tissue from input parameters
## ./cuffdiff_1000.sh species tissue timepoint
species=$1
tissue=$2
comparisons=$3


## create merged output for genes and isoforms
perl avg_row_genes.pl "${output_dir}/${species}/${tissue}/cuffdiff_output_${comparisons}/" > "${output_dir}/${species}/${tissue}/cuffdiff_output_${comparisons}/genes.fpkm_merged.txt"

perl avg_row_isoforms.pl "${output_dir}/${species}/${tissue}/cuffdiff_output_${comparisons}/" > "${output_dir}/${species}/${tissue}/cuffdiff_output_${comparisons}/isoforms.fpkm_merged.txt"

## remove all other files
rm -rf ${output_dir}/${species}/${tissue}/cuffdiff_output_${comparisons}/run*	
