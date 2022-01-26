#!/bin/bash
#$ -S /bin/bash

srxfile=/home/s-tahara/DROMPA/code_hg38_2/SRX_TO_SRR_hg19.txt
MEMORY_REQ=24G

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
CMD=/home/s-tahara/DROMPA/code_hg38_2/drompa_ssp_hg38_ddbj_v2.sh
#JOB_NUM=`echo ${srx_array[@]} | tr ' ' '\n' | uniq | wc -l | awk '{print $1}'`
JOB_NUM=${#srx_array_uniq[@]}
echo after sed JOB_NUM ${JOB_NUM}

mkdir -p /home/s-tahara/DROMPA/code_hg38_2/array_log/
qsub -cwd -o /home/s-tahara/DROMPA/code_hg38_2/array_log/ -e /home/s-tahara/DROMPA/code_hg38_2/array_log/ -t 1:${JOB_NUM} -V -l s_vmem=${MEMORY_REQ} -l mem_req=${MEMORY_REQ} ${CMD}
