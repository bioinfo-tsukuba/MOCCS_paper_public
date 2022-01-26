#!/bin/bash
#$ -S /bin/bash

singularity run lolcow.sif 
singularity exec --bind tmp:/tmp --bind singularity_r_package:/usr/local/lib/R/site-library singularity-rstudio-tidyverse.simg Rscript --no-save /home/s-tahara/allele_binding_SNP_hg38/GWAS/script/job_SLE.R  > /home/s-tahara/singularity_GWAS_ver4/job_log/job_SLE.log
