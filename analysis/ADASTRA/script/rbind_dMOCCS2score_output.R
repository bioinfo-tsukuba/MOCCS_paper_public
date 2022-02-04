rm(list=ls())
target_tf <- "GATA3"
ID_list <- read_tsv(paste0("/home/s-tahara/singularity_ADASTRA/ID_list/", target_tf, ".txt"), col_names = F)
length(ID_list$X1)

file_path_list <- list.files("/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/data/dMOCCS2score_output", pattern=paste0("ADASTRA_dMOCCS2score_", target_tf))
length(file_path_list)

for (i in 1:length(file_path_list)) {
  target_file_path <- file_path_list[i]
  target_ID_pre <- gsub("ADASTRA_dMOCCS2score_GATA3_", "", target_file_path)
  target_ID <- gsub("_qval.rds", "",target_ID_pre)
  target_df1 <- readRDS(paste0("/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/data/dMOCCS2score_output/", target_file_path))
  target_df2 <- target_df1 %>% mutate(ID = rep(target_ID, nrow(target_df1)))
  if(i == 1){
    df_all <- target_df2
  }else{
    df_all <- df_all %>% add_row(target_df2)
  }
}
df_all <- df_all %>% drop_na(dMOCCS2score) 

saveRDS(df_all, paste0("/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/data/dMOCCS2score_output_all/ADASTRA_dMOCCS2score_", target_tf, "_all.rds"))
