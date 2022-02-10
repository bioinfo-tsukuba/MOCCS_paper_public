#!/bin/bash
#$ -S /bin/bash

PROJECT_PATH=/home/s-tahara/MOCCS-DB
RESOURCE_PATH_FAI=/home/s-tahara/MOCCS-DB/RESOURCE/UCSC/hg38/hg38.chrom.sizes
RESOURCE_PATH=/home/s-tahara/MOCCS-DB/RESOURCE/UCSC/hg38/hg38.fa
ChIP_ATLAS_BED_PARH=/home/s-tahara/DROMPA/code_hg38_1/ftp.biosciencedbc.jp/archive/chip-atlas/data/hg38/eachData/bed05/

#srxfile=${PROJECT_PATH}/WORK/MOCCS/MOCCS-SCRIPT/check_ID_list_hg19_TF.txt
srxfile=${PROJECT_PATH}/WORK/MOCCS/MOCCS-SCRIPT/Antigen_list_hg19.txt


declare -a srx_array=()
row_N=`cat ${srxfile} | wc -l`

for row_num in `seq 1 ${row_N}`
do
 line=`head -n ${row_num} ${srxfile} | tail -n 1`
 SRX=$(echo ${line} | awk '{print $1}')
 srx_array+=(${SRX})
done

declare -a srx_array_uniq=()
srx_array_uniq=(`echo ${srx_array[@]} | tr ' ' '\n' | uniq | xargs` )

#for output
mkdir -p ${PROJECT_PATH}/WORK/MOCCS/MOCCS-OUTPUT-HG38_v2/${srx_array_uniq[((${SGE_TASK_ID}-1))]}/
mkdir -p ${PROJECT_PATH}/WORK/LIST_DUMP_HG38_v2/${srx_array_uniq[((${SGE_TASK_ID}-1))]}/
echo ${srx_array_uniq[((${SGE_TASK_ID}-1))]}

bed_file_path=/home/s-tahara/DROMPA/code_hg38_1/ftp.biosciencedbc.jp/archive/chip-atlas/data/hg38/eachData/bed05/${srx_array_uniq[((${SGE_TASK_ID}-1))]}.05.bed
#bed file 1bp
awk '{OFS="\t"; c=int(($3+$2)/2); $2=c; $3=c+1; print $0}' ${bed_file_path} > ${PROJECT_PATH}/WORK/LIST_DUMP_HG38_v2/${srx_array_uniq[((${SGE_TASK_ID}-1))]}/${srx_array_uniq[((${SGE_TASK_ID}-1))]}_new_6mer.bed

#bed file 701bp
bedtools slop -i ${PROJECT_PATH}/WORK/LIST_DUMP_HG38_v2/${srx_array_uniq[((${SGE_TASK_ID}-1))]}/${srx_array_uniq[((${SGE_TASK_ID}-1))]}_new_6mer.bed -g ${RESOURCE_PATH_FAI} -b 350 > ${PROJECT_PATH}/WORK/LIST_DUMP_HG38_v2/${srx_array_uniq[((${SGE_TASK_ID}-1))]}/${srx_array_uniq[((${SGE_TASK_ID}-1))]}_ext_6mer.bed

#change to fasta format
bedtools getfasta -fi ${RESOURCE_PATH} -bed ${PROJECT_PATH}/WORK/LIST_DUMP_HG38_v2/${srx_array_uniq[((${SGE_TASK_ID}-1))]}/${srx_array_uniq[((${SGE_TASK_ID}-1))]}_ext_6mer.bed > ${PROJECT_PATH}/WORK/LIST_DUMP_HG38_v2/${srx_array_uniq[((${SGE_TASK_ID}-1))]}/${srx_array_uniq[((${SGE_TASK_ID}-1))]}_ext_6mer.fasta

#MOCCS
perl /home/s-tahara/MOCCS-DB/WORK/MOCCS/moccs/MOCCS.pl -i ${PROJECT_PATH}/WORK/LIST_DUMP_HG38_v2/${srx_array_uniq[((${SGE_TASK_ID}-1))]}/${srx_array_uniq[((${SGE_TASK_ID}-1))]}_ext_6mer.fasta -k 6 --label ${srx_array_uniq[((${SGE_TASK_ID}-1))]}_6mer_v2 --mask --low-count-threshold -1

mv ${srx_array_uniq[((${SGE_TASK_ID}-1))]}_6mer_v2.* ${PROJECT_PATH}/WORK/MOCCS/MOCCS-OUTPUT-HG38_v2/${srx_array_uniq[((${SGE_TASK_ID}-1))]}/


