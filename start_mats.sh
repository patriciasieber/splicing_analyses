#!/bin/bash

output_dir="/mnt/dessertlocal/emanuel/hacken_alternative_splicing/differential_splicing/"


## input
## get species and tissue from input parameters
## ./start_mats.sh species tissue timepoint1,timepoint2

species=$1
tissue=$2
comparison=$3

t1="$(cut -d',' -f1 <<<"$comparison")"
t2="$(cut -d',' -f2 <<<"$comparison")"

path1="/mnt/dessertlocal/emanuel/hacken_alternative_splicing/mappings/"


## read time points and files
	
case "$species" in
dre) 	gtf_file="/home/prak27/Danio_rerio.GRCz10.release83.annotation.modified.gtf"
		;;
hsa) 	gtf_file="/home/prak27/Homo_sapiens.GRCh38.85.chr.sorted.nochr.modified.gtf"
		;;
mmu) 	gtf_file="/home/prak27/Mus_musculus.GRCm38.87.chr.sorted.chr.modified.gtf"
		;;
nfu) 	gtf_file="/home/prak27/nothobranchius_furzeri_annotation.modified.gtf"
		;;
esac

file="$path1/$species/$tissue/"
echo $file

cond1="$file/$t1"
shopt -s nullglob
files_cond1=($cond1/*.bam)

cond2="$file/$t2"
shopt -s nullglob
files_cond2=($cond2/*.bam)

all_files_cond1="${files_cond1[0]}"
cond1Len=${#files_cond1[@]}
for (( j=1; j<$((${cond1Len})); j++ ));
do    		
	all_files_cond1=",$all_files_cond1${files_cond1[j]}"
done

all_files_cond2="${files_cond2[0]}"
cond2Len=${#files_cond2[@]}
for (( j=1; j<$((${cond2Len})); j++ ));
do    
	all_files_cond2=",$all_files_cond2${files_cond2[j]}"
done

t1short="$(cut -d'_' -f1 <<<"$t1")"
t2short="$(cut -d'_' -f1 <<<"$t2")"
comparisons="${t1short}vs${t2short}"

## run mats
python /mnt/prostlocal/programs/rmats/3.0.8/RNASeq-MATS.py -b1 ${all_files_cond1} -b2 ${all_files_cond2} -gtf ${gtf_file} -o "${output_dir}/${species}/${tissue}/${comparisons}" -t single -len 50

rm -rf "${output_dir}/${species}/${tissue}/${comparisons}/SAMPLE_1"
rm -rf "${output_dir}/${species}/${tissue}/${comparisons}/SAMPLE_2"







