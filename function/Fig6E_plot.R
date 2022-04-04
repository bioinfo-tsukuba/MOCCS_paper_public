Fig6E_plot <- function(target_phenotype, target_rs, target_position, annotation_path, threshold){
  
  library(tidyverse)
  target_df <- readRDS(paste0(annotation_path, "result_output_binded_all/", target_phenotype, "_peak_rand_binded_all.rds"))
  
  totalization_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds"
  totalization <- readRDS(totalization_path)
  annotation <- totalization %>% select(ID, Antigen, Cell_type_class, Cell_type) %>% distinct()
  df_phenotype_binded_all_selected_annotated <- target_df %>% left_join(annotation, by = "ID")
  df1 <- df_phenotype_binded_all_selected_annotated %>% mutate(abs_dMOCCS2score = abs(dMOCCS2score)) %>% select(ID, p_value, Cell_type_class, Cell_type, Antigen, abs_dMOCCS2score) %>% distinct()
  
  ID_snp_for_join <- df_phenotype_binded_all_selected_annotated %>% 
    filter(snp == target_position) %>% 
    unite("ID_position", c(ID, snp)) %>% 
    .$ID_position %>% 
    as.character()
  df_phenotype_binded_all_selected_annotated2 <- df_phenotype_binded_all_selected_annotated %>% 
    filter(snp == target_position) %>%
    mutate(ID_position = ID_snp_for_join) %>% 
    drop_na(dMOCCS2score) 
  df2 <- df_phenotype_binded_all_selected_annotated2 %>% mutate(abs_dMOCCS2score = abs(dMOCCS2score)) %>% 
    select(ID_position, snp, abs_dMOCCS2score, ID, Antigen, Cell_type_class, Cell_type) %>%
    distinct() %>%
    group_by(ID, Antigen, Cell_type_class, snp) %>% 
    summarise(max_abs_dMOCCS2score = max(abs_dMOCCS2score))
  
  x_lab_ID <- df2  %>% filter(max_abs_dMOCCS2score > threshold) %>% arrange(desc(max_abs_dMOCCS2score)) %>% .$ID %>% as.character()
  Fig6D_CL <- df2  %>% unite("ID_position", c(ID, snp)) %>% filter(max_abs_dMOCCS2score > threshold) %>%
    ggplot(aes(x = reorder(ID_position, max_abs_dMOCCS2score), y = max_abs_dMOCCS2score, fill = Cell_type_class))+
    geom_col() +
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_line(colour="gray"),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=6,face="bold", angle = 45, hjust = 1),
          axis.text.y =element_text(size=10,face="bold"),
          axis.title=element_text(size=14,face="bold"),
          #legend.position = 'none',
          legend.title = element_blank()
    )+
    scale_x_discrete("ID", labels = x_lab_ID) +
    ylab("max(|dMOCCS2score|)") +
    ggtitle(paste0(target_rs, " ", target_phenotype)) +
    coord_flip()
  
  Fig6D_TF <- df2  %>% unite("ID_position", c(ID, snp)) %>% filter(max_abs_dMOCCS2score > threshold) %>%
    ggplot(aes(x = reorder(ID_position, max_abs_dMOCCS2score), y = max_abs_dMOCCS2score, fill = Antigen))+
    geom_col() +
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_line(colour="gray"),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=6,face="bold", angle = 45, hjust = 1),
          axis.text.y =element_text(size=10,face="bold"),
          axis.title=element_text(size=14,face="bold"),
          #legend.position = 'none',
          legend.title = element_blank()
    )+
    scale_x_discrete("ID", labels = x_lab_ID) +
    ylab("max(|dMOCCS2score|)") +
    ggtitle(paste0(target_rs, " ", target_phenotype)) +
    coord_flip()
  
  return(list(Fig6D_CL, Fig6D_TF))
  
}