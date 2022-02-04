#!/bin/sh
#$ -S /bin/sh
#$ -l medium
#$ -cwd
#$ -t 1-1
#$ -tc 20
#$ -l s_vmem=4G -l mem_req=4G


srxfile=/home/s-tahara/singularity_ADASTRA/ID_list/GATA3.txt
MEMORY_REQ=4G

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


echo HOME: $HOME
echo USER: $USER
echo JOB_ID: $JOB_ID
echo JOB_NAME: $JOB_NAME
echo HOST_NAME: $HOSTNAME
echo SGE_TASK_ID: ${srx_array_uniq[((${SGE_TASK_ID}-1))]} $SGE_TASK_ID
echo SGE_TASK_FIRST: $SGE_TASK_FIRST
echo SGE_TASK_LAST: $SGE_TASK_LAST
echo SGE_TASK_STEPSIZE: $SGE_TASK_STEPSIZE
pwd 

#ARGS="'--args x=${SGE_TASK_ID}'"
ARGS="'--args ${srx_array_uniq[((${SGE_TASK_ID}-1))]}'"
echo ARGS:$ARGS

singularity run lolcow.sif 
#R_CMD="R CMD BATCH --slave --vanilla ${ARGS} /home/s-tahara/allele_binding_SNP_hg38/ADASTRA/script/test.R log/test.R.log.${srx_array_uniq[((${SGE_TASK_ID}-1))]}"
#R_CMD="singularity exec --bind tmp:/tmp --bind singularity_r_package:/usr/local/lib/R/site-library singularity-rstudio-tidyverse.simg Rscript --no-save --slave --vanilla --args ${srx_array_uniq[((${SGE_TASK_ID}-1))]} /home/s-tahara/allele_binding_SNP_hg38/ADASTRA/script/test.R > /home/s-tahara/singularity_ADASTRA/log/test_singularity.R.log.${srx_array_uniq[((${SGE_TASK_ID}-1))]}"


# --args以降をとるとうまく言ったやつ
#singularity exec --bind tmp:/tmp --bind singularity_r_package:/usr/local/lib/R/site-library singularity-rstudio-tidyverse.simg Rscript --no-save --vanilla --default-packages=utils --args ${srx_array_uniq[((${SGE_TASK_ID}-1))]} /home/s-tahara/allele_binding_SNP_hg38/ADASTRA/script/test.R > /home/s-tahara/singularity_ADASTRA/log/test_singularity.R.log.${srx_array_uniq[((${SGE_TASK_ID}-1))]}
singularity exec --bind tmp:/tmp --bind singularity_r_package:/usr/local/lib/R/site-library singularity-rstudio-tidyverse.simg Rscript --no-save --slave --vanilla /home/s-tahara/allele_binding_SNP_hg38/ADASTRA/script/test.R ${srx_array_uniq[((${SGE_TASK_ID}-1))]} > /home/s-tahara/singularity_ADASTRA/log/test_singularity.R.log.${srx_array_uniq[((${SGE_TASK_ID}-1))]}

#echo $R_CMD
#eval $R_CMD

echo "Finished."
