Read_df <- function(){

  library(dplyr)
  
  # Get MOCCSout_hg38_hard_filter_annotated.rds
  system("wget https://figshare.com/ndownloader/files/34065695 -P ~/MOCCS-DB_paper/data/Fig2/out")
  system("mv ~/MOCCS-DB_paper/data/Fig2/out/34065695 ~/MOCCS-DB_paper/data/Fig2/out/MOCCSout_hg38_hard_filter_annotated.rds")
  
  # Get MOCCSout_hg38_all_qval.rds
  system("wget https://figshare.com/ndownloader/files/34065689 -P ~/MOCCS-DB_paper/data/Fig2/qval")
  system("mv ~/MOCCS-DB_paper/data/Fig2/qval/34065689 ~/MOCCS-DB_paper/data/Fig2/qval/MOCCSout_hg38_all_qval.rds")
  
  df_raw_pre <- readRDS("~/MOCCS-DB_paper/data/Fig2/out/MOCCSout_hg38_hard_filter_annotated.rds")
  df_raw_pre %>% mutate(ID_kmer = paste0(ID, kmer)) -> df_raw_pre_2
  df_qval <- readRDS("~/MOCCS-DB_paper/data/Fig2/qval/MOCCSout_hg38_all_qval.rds")
  df_qval %>% mutate(ID_kmer = paste0(ID, kmer)) -> df_qval_2
  df_qval_2$q_value[is.na(df_qval_2$q_value)] <- 1
  df_raw <- inner_join(df_raw_pre_2, df_qval_2[,c("q_value", "ID_kmer")], by = "ID_kmer")
  return(df_raw)  
  
}