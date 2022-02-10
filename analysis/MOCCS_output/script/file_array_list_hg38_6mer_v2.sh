#!/bin/bash
#$ -S /bin/bash

PROJECT_PATH=/home/s-tahara/MOCCS-DB/
mkdir -p ${PROJECT_PATH}WORK/LIST_DUMP_HG38_v2
LIST_DUMP=${PROJECT_PATH}WORK/LIST_DUMP_HG38_v2
MEMORY_REQ=24G

declare -a bed_id_array=()
declare -a bed_file_array=()

#srxfile=${PROJECT_PATH}WORK/MOCCS/MOCCS-SCRIPT/check_ID_list_hg19_TF.txt
srxfile=${PROJECT_PATH}/WORK/MOCCS/MOCCS-SCRIPT/Antigen_list_hg19.txt

declare -a srx_array=()
row_N=`cat ${srxfile} | wc -l`

for row_num in `seq 1 ${row_N}`
do
 line=`head -n ${row_num} ${srxfile} | tail -n 1`
 SRX=$(echo ${line} | awk '{print $1}')
 srx_array+=(${SRX})
done
#echo ${srx_array[@]}

declare -a srx_array_uniq=()
srx_array_uniq=(`echo ${srx_array[@]} | tr ' ' '\n' | uniq | xargs` )
echo ${srx_array_uniq[1]}
echo ${srx_array_uniq[2]}

JOB_NUM=${#srx_array_uniq[*]}
CMD=/home/s-tahara/MOCCS-DB/WORK/MOCCS/MOCCS-SCRIPT/getfast_array_list_hg38_6mer_v2.sh

qsub -cwd -o ${LIST_DUMP}/file_array_20220121_hg38_o.log -e ${LIST_DUMP}/file_array_20210121_hg38_e.log -tc 50 -t 1:${JOB_NUM} -V -l s_vmem=${MEMORY_REQ} -l mem_req=${MEMORY_REQ} ${CMD}


