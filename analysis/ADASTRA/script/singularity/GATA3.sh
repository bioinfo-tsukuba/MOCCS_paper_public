#!/bin/sh
#$ -S /bin/sh

target_tf=GATA3
srxfile=/home/s-tahara/singularity_ADASTRA/ID_list/${target_tf}.txt
MEMORY_REQ=32G

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

singularity run lolcow.sif 
singularity exec --bind tmp:/tmp --bind singularity_r_package:/usr/local/lib/R/site-library singularity-rstudio-tidyverse.simg Rscript --no-save --slave --vanilla /home/s-tahara/allele_binding_SNP_hg38/ADASTRA/script/job_ADASTRA.R ${srx_array_uniq[((${SGE_TASK_ID}-1))]} ${target_tf} > /home/s-tahara/singularity_ADASTRA/log/${target_tf}.R.log.${srx_array_uniq[((${SGE_TASK_ID}-1))]}
