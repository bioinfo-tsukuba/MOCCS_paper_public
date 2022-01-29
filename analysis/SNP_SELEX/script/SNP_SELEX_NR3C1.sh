#!/bin/bash
#$ -S /bin/bash

target_TF=NR3C1
singularity run lolcow.sif 
singularity exec --bind tmp:/tmp --bind singularity_r_package:/usr/local/lib/R/site-library singularity-rstudio-tidyverse.simg Rscript --no-save /home/s-tahara/allele_binding_SNP_hg38/SNP_SELEX/script/exec_get_dMOCCS2score_SNP_SELEX_${target_TF}.R  > /home/s-tahara/singularity_SNP_SELEX/job_log/job_${target_TF}.log
