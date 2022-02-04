library(tidyverse)
df <- readRDS("/home/s-tahara/MOCCS-DB/WORK/MOCCS/MOCCS-ANALYSIS/output_hg38/MOCCSout_hg38_hard_filter_annotated.rds")
df2 <- df %>% select(ID, Antigen) %>% distinct()
target_tf_list <- c("CTCF", "FOXA1", "BRD4", "ESR1", "SPI1", "EP300", "JUND", "MYCN", "MAX", "EZH2", "SUMO2", "CEBPB", 
                    "GATA2", "GATA3", "HDAC2", "FOS")

for (i in 1:length(target_tf_list)) {
  
  target_tf <- target_tf_list[i]
  target_ID_list <- df2 %>% filter(Antigen == target_tf) %>% .$ID %>% unique() %>% as.character()
  ID_df <- tibble(ID = target_ID_list)
  write_tsv(ID_df, paste0("/home/s-tahara/singularity_ADASTRA/ID_list/", target_tf, ".txt"),col_names = FALSE)
  
}







