rm(list=ls())
library(tidyverse)
source("/home/s-tahara/allele_binding_SNP_hg38/SNP_SELEX/script/get_dMOCCS2score_SNP_SELEX.R")

input_bed_file_path = "/home/s-tahara/DROMPA/code_hg38_1/ftp.biosciencedbc.jp/archive/chip-atlas/data/hg38/eachData/bed05/"
input_allele_binding_SNP_file_path = "/home/s-tahara/allele_binding_SNP_hg38/SNP_SELEX/data/GSE118725_pbs_novel_batch_hg38.tsv"
bed_output_dir = "/home/s-tahara/allele_binding_SNP_hg38/SNP_SELEX/data/bed_output"
k_mer_num = 6
ref_genome = "HG38"
input_MOCCS_file_dir = "/home/s-tahara/MOCCS-DB/WORK/MOCCS/"
fasta_file_name = "/home/s-tahara/MOCCS-DB/RESOURCE/UCSC/hg38/hg38.fa"

target_tf <- "GATA3"
totalization <- readRDS("/home/s-tahara/MOCCS-DB/WORK/MOCCS/MOCCS-ANALYSIS/output_hg38_v2/MOCCSout_hg38_all_qval_annotated.rds")
ID_hard <- readRDS("/home/s-tahara/DROMPA/code_hg38_2/hg38_hard_filter_ID.rds")
Antigen_list <- totalization %>% filter(ID %in% ID_hard) %>%filter(Antigen == target_tf) %>%  .$ID %>% unique() %>% as.character() 

dMOCCS2score_df <- list()
for (x in 1:length(Antigen_list)) {
  input_SRX_ID <- Antigen_list[x]
  dMOCCS2score <- get_dMOCCS2score(input_bed_file_path,
                                   input_allele_binding_SNP_file_path,
                                   bed_output_dir,
                                   k_mer_num,
                                   ref_genome, 
                                   target_tf, 
                                   input_SRX_ID,
                                   input_MOCCS_file_dir,
                                   fasta_file_name 
  )
  dMOCCS2score_df[[input_SRX_ID]] <-dMOCCS2score 
}
saveRDS(dMOCCS2score_df, paste0("/home/s-tahara/allele_binding_SNP_hg38/SNP_SELEX/data/dMOCCS2score/dMOCCS2score_hg38_", target_tf, ".rds"))


df_all <-c()
for(i in 1:length(dMOCCS2score_df)){
  target_ID  <- Antigen_list[[i]]
  target_df <- dMOCCS2score_df[[target_ID]]
  if(target_df == "no common region" | target_df == "no target TF"){
    print("skip")
  }else{
    target_df <- target_df %>% mutate(ID = target_ID)
    df_all <- rbind(df_all, target_df)
  }
}
write_tsv(df_all, paste0("/home/s-tahara/allele_binding_SNP_hg38/SNP_SELEX/data/dMOCCS2score/dMOCCS2score_hg38_", target_tf, "_all.tsv"))