FigS8_plot <- function(path, target_phenotype, target_tf){
  
  library(tidyverse)
  target_df <- readRDS(paste0(path, "result_output_binded_all/", target_phenotype, "_peak_rand_binded_all.rds"))
  
  totalization_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds"
  totalization <- readRDS(totalization_path)
  annotation <- totalization %>% select(ID, Antigen, Cell_type_class, Cell_type) %>% distinct()
  df_phenotype_binded_all_selected_annotated <- target_df %>% left_join(annotation, by = "ID")
  
  df1 <- df_phenotype_binded_all_selected_annotated %>% filter(Antigen == target_tf)
}