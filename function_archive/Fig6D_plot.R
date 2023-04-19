Fig6D_plot <- function(target_phenotype, annotation_path, threshold){
  
  library(tidyverse)
  if(target_phenotype == "CD"){
    target_df <- readRDS(url("https://figshare.com/ndownloader/files/34065827", "rb"))
  }else if(target_phenotype == "MS"){
    target_df <- readRDS(url("https://figshare.com/ndownloader/files/34065830", "rb"))
  }else if(target_phenotype == "SLE"){
    target_df <- readRDS(url("https://figshare.com/ndownloader/files/34660126", "rb"))
  }else if(target_phenotype == "IBD"){
    target_df <- readRDS(url("https://figshare.com/ndownloader/files/34660132", "rb"))
  }
  
  totalization <- readRDS(url("https://figshare.com/ndownloader/files/34065686","rb")) #MOCCSout_hg38_all_qval_annotated.rds
  annotation <- totalization %>% select(ID, Antigen, Cell_type_class, Cell_type) %>% distinct()
  df_phenotype_binded_all_selected_annotated <- target_df %>% left_join(annotation, by = "ID") %>% filter(q_value < 0.05)
  df1 <- df_phenotype_binded_all_selected_annotated %>% mutate(abs_dMOCCS2score = abs(dMOCCS2score)) %>% select(ID, p_value, Cell_type_class, Cell_type, Antigen, abs_dMOCCS2score) %>% distinct()
  
  ID_snp_for_join <- df_phenotype_binded_all_selected_annotated %>% unite("ID_snp_for_join", c(ID, position)) %>% .$ID_snp_for_join %>% as.character()
  df_phenotype_binded_all_selected_annotated2 <- df_phenotype_binded_all_selected_annotated %>% mutate(ID_snp_join = ID_snp_for_join)
  df2 <- df_phenotype_binded_all_selected_annotated2 %>% mutate(abs_dMOCCS2score = abs(dMOCCS2score)) %>% group_by(ID, position, Antigen, Cell_type_class, Cell_type) %>% summarise(max_abs_dMOCCS2score = max(abs_dMOCCS2score))
  
  x_lab_ID <- df2  %>% unite("ID_position", c(ID, position)) %>% filter(max_abs_dMOCCS2score > threshold) %>% arrange(desc(max_abs_dMOCCS2score)) %>% select(ID_position)
  x_lab_ID2 <- x_lab_ID %>% separate(ID_position, into = c("ID", "chr", "position"), sep = "_") %>% .$ID %>% as.character() %>% rev()
  Fig6C_CL <- df2  %>% unite("ID_position", c(ID, position)) %>% filter(max_abs_dMOCCS2score > threshold) %>%
    ggplot(aes(x = reorder(ID_position, max_abs_dMOCCS2score), y = max_abs_dMOCCS2score, fill = Cell_type_class))+
    geom_col() +
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=6,face="bold", angle = 45, hjust = 1),
          axis.text.y =element_text(size=10,face="bold"),
          axis.title=element_text(size=14,face="bold"),
          #legend.position = 'none',
          legend.title = element_blank()
    )+
    scale_x_discrete("ID", labels = x_lab_ID2) +
    ylab("max(|dMOCCS2score|)") +
    ggtitle(target_phenotype) +
    coord_flip()
  
  Fig6C_TF <- df2  %>% unite("ID_position", c(ID, position)) %>% filter(max_abs_dMOCCS2score > threshold) %>%
    ggplot(aes(x = reorder(ID_position, max_abs_dMOCCS2score), y = max_abs_dMOCCS2score, fill = Antigen))+
    geom_col() +
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=6,face="bold", angle = 45, hjust = 1),
          axis.text.y =element_text(size=10,face="bold"),
          axis.title=element_text(size=14,face="bold"),
          #legend.position = 'none',
          legend.title = element_blank()
    )+
    scale_x_discrete("ID", labels = x_lab_ID2) +
    ylab("max(|dMOCCS2score|)") +
    ggtitle(target_phenotype) +
    coord_flip()

  return(list(Fig6C_CL, Fig6C_TF))
  
}
