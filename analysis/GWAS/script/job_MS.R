library(tidyverse)

# 1. load ID list
ID_hard <- readRDS("/home/s-tahara/DROMPA/code_hg38_2/hg38_hard_filter_ID.rds")

# 2. prepare progress bar
library("progress")
n <- length(ID_hard)
pb <- progress_bar$new(total = n,
                       format = "[:bar] :percent : :eta",
                       clear = TRUE)


# 3. prepare function
target_phenotype <- "MS"
source(paste0("/home/s-tahara/allele_binding_SNP_hg38/GWAS/script/get_dMOCCS2score_GWAS_peaks_rand_phenotype_", target_phenotype,".R"))

input_bed_file_path <- "/home/s-tahara/DROMPA/code_hg38_1/ftp.biosciencedbc.jp/archive/chip-atlas/data/hg38/eachData/bed05/"
bed_output_dir = paste0("/home/s-tahara/allele_binding_SNP_hg38/GWAS/data/bed_output_", target_phenotype)
fasta_file_name =  "/home/s-tahara/MOCCS-DB/RESOURCE/UCSC/hg38/hg38.fa"
input_MOCCS_file_dir = "/home/s-tahara/MOCCS-DB/WORK/MOCCS/"
result_output_dir = paste0("/home/s-tahara/allele_binding_SNP_hg38/GWAS/data/result_output_",target_phenotype,"/")
k_mer_num = 6
simu_N = 20
ref_genome = "HG38"

# 4. apply function to all hard filter samples
for (i in 1:n) {
  target_ID <- ID_hard[i]
  print(paste0(i, " ", target_ID))
  df_phenotype <- get_dMOCCS2score_GWAS_peaks_rand_phenotype(target_phenotype, 
                                                             target_ID,
                                                             input_bed_file_path, 
                                                             bed_output_dir, 
                                                             fasta_file_name, 
                                                             input_MOCCS_file_dir, 
                                                             result_output_dir,
                                                             k_mer_num, 
                                                             simu_N, 
                                                             ref_genome)
  
  if(is.null(df_phenotype) == FALSE){
    saveRDS(df_phenotype, paste0("/home/s-tahara/allele_binding_SNP_hg38/GWAS/data/result_output_", target_phenotype,"/",target_ID, "_", target_phenotype,"_peak_rand.rds"))
  }
}
