#!/bin/bash

## modify:
output_dir="/home/cuffdiff_test"
numThreads=60

## input
## get species and tissue from input parameters
## ./cuffdiff_1000.sh species tissue timepoint1,timepoint2 runID
species=$1
tissue=$2
comparison=$3
runID=$4 

t1="$(cut -d',' -f1 <<<"$comparison")"
t2="$(cut -d',' -f2 <<<"$comparison")"

path1="/mnt/fass2/hacken_2017/"
path2="/mappings/polyA/tophat2/normal_aging/"

## read time points and files
	
case "$species" in
danio_rerio) 	fasta_file="/mnt/aging_splicing/danio_rerio/genome/Danio_rerio.GRCz10.genome.fa"
		gtf_file="/mnt/aging_splicing/danio_rerio/annotation/Danio_rerio.GRCz10.release83.annotation.gtf"
		;;
homo_sapiens) 	fasta_file="/mnt/aging_splicing/homo_sapiens/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
		gtf_file="/mnt/aging_splicing/homo_sapiens/annotation/Homo_sapiens.GRCh38.85.chr.sorted.gtf"
		;;
mus_musculus) 	fasta_file="/mnt/aging_splicing/mus_musculus/genome/Mus_musculus.GRCm38.fa"
		gtf_file="/mnt/aging_splicing/mus_musculus/annotation/Mus_musculus.GRCm38.87.chr.sorted.chr.gtf"
		;;
nothobranchius_furzeri) fasta_file="/mnt/aging_splicing/nothobranchius_furzeri/genome/nothobranchius_furzeri_platzer.fa"
		gtf_file="/mnt/aging_splicing/nothobranchius_furzeri/annotation/nothobranchius_furzeri_annotation.sorted.gff"
		;;
esac
	

file="$path1/$species/$path2/$tissue/"
echo $file
mkdir -p "$output_dir/$species/$tissue/"

## compare each time point with its successor, the last one against the first

cond1="$file/$t1"
shopt -s nullglob
files_cond1=($cond1/*.bam)

cond2="$file/$t2"
comparisons="${t1}vs${t2}"
shopt -s nullglob
files_cond2=($cond2/*.bam)

all_files_cond1="${files_cond1[0]},"
cond1Len=${#files_cond1[@]}
for (( j=1; j<$((${cond1Len}-1)); j++ ));
do    		
	all_files_cond1="$all_files_cond1${files_cond1[j]},"
done
all_files_cond1="$all_files_cond1${files_cond1[$((${cond1Len}-1))]}"

all_files_cond2="${files_cond2[0]},"
cond2Len=${#files_cond2[@]}
for (( j=1; j<$((${cond2Len}-1)); j++ ));
do    
	all_files_cond2="$all_files_cond2${files_cond2[j]},"
done
all_files_cond2="$all_files_cond2${files_cond2[$((${cond2Len}-1))]}"

## run cuffdiff
cuffdiff -o "${output_dir}/${species}/${tissue}/cuffdiff_output_${comparisons}/run$runID" -b ${fasta_file} -p $numThreads -L "${t1},${t2}" -u ${gtf_file} ${all_files_cond1} ${all_files_cond2} &	






