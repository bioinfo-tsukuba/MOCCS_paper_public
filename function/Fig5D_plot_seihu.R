Fig5D_plot_seihu <- function(target_phenotype, target_rs, target_position,  threshold){
  
  library(tidyverse)
  options(timeout=1000)
  if(target_phenotype == "CD"){
    #target_df <- readRDS(url("https://figshare.com/ndownloader/files/34065827", "rb"))
    target_df <- readRDS("/Users/saeko/MOCCS_paper_public/data/Fig6/CD_peak_rand_binded_all_qval.rds")
  }else if(target_phenotype == "MS"){
    target_df <- readRDS(url("https://figshare.com/ndownloader/files/34065830", "rb"))
  }else if(target_phenotype == "SLE"){
    target_df <- readRDS(url("https://figshare.com/ndownloader/files/34660126", "rb"))
  }else if(target_phenotype == "IBD"){
    target_df <- readRDS(url("https://figshare.com/ndownloader/files/34660132", "rb"))
  }
  
  #totalization <- readRDS(url("https://figshare.com/ndownloader/files/34065686","rb")) #MOCCSout_hg38_all_qval_annotated.rds
  totalization <- readRDS("~/MOCCS_paper_public/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds")
  annotation <- totalization %>% select(ID, Antigen, Cell_type_class, Cell_type) %>% distinct()
  df_phenotype_binded_all_selected_annotated <- target_df %>% left_join(annotation, by = "ID")
  
  ID_snp_for_join <- df_phenotype_binded_all_selected_annotated %>% 
    #filter(snp == target_position) %>% 
    filter(position == target_position) %>% 
    #unite("ID_position", c(ID, snp)) %>% 
    unite("ID_position", c(ID, position)) %>% 
    .$ID_position %>% 
    as.character()
  df_phenotype_binded_all_selected_annotated2 <- df_phenotype_binded_all_selected_annotated %>% 
    #filter(snp == target_position) %>%
    filter(position == target_position) %>% 
    mutate(ID_position = ID_snp_for_join) %>% 
    drop_na(dMOCCS2score) 
  df2 <- df_phenotype_binded_all_selected_annotated2 %>% mutate(abs_dMOCCS2score = abs(dMOCCS2score)) %>% 
    #select(ID_position, snp, abs_dMOCCS2score, ID, Antigen, Cell_type_class, Cell_type) %>%
    select(ID_position, position, abs_dMOCCS2score, ID, Antigen, Cell_type_class, Cell_type) %>%
    distinct() %>%
    #group_by(ID, Antigen, Cell_type_class, snp) %>% 
    group_by(ID, Antigen, Cell_type_class, position) %>% 
    summarise(max_abs_dMOCCS2score = max(abs_dMOCCS2score))
  
  #x_lab_ID <- df2  %>% filter(max_abs_dMOCCS2score > threshold) %>% arrange(desc(max_abs_dMOCCS2score)) %>% .$ID %>% as.character()
  #x_lab_ID <- df2  %>% unite("ID_position", c(ID, snp)) %>% filter(max_abs_dMOCCS2score > threshold) %>% arrange(desc(max_abs_dMOCCS2score)) %>% select(ID_position)
  x_lab_ID <- df2  %>% unite("ID_position", c(ID, position)) %>% filter(max_abs_dMOCCS2score > threshold) %>% arrange(desc(max_abs_dMOCCS2score)) %>% select(ID_position)
  x_lab_ID2 <- x_lab_ID %>% separate(ID_position, into = c("ID", "chr", "position"), sep = "_") %>% .$ID %>% as.character() %>% rev()
  
  Fig6E_CL <- df2  %>% 
    #unite("ID_position", c(ID, snp)) %>%
    unite("ID_position", c(ID, position)) %>%
    filter(max_abs_dMOCCS2score > threshold) %>%
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
    scale_x_discrete("ID", labels = x_lab_ID2) +
    ylab("max(|dMOCCS2score|)") +
    ggtitle(paste0(target_rs, " ", target_phenotype)) +
    coord_flip()
  
  Fig6E_TF <- df2  %>% 
    #unite("ID_position", c(ID, snp)) %>% 
    unite("ID_position", c(ID, position)) %>%
    filter(max_abs_dMOCCS2score > threshold) %>%
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
    scale_x_discrete("ID", labels = x_lab_ID2) +
    ylab("max(|dMOCCS2score|)") +
    ggtitle(paste0(target_rs, " ", target_phenotype)) +
    coord_flip()
  
  # dMOCCS2score >0 -----
  df_phenotype_binded_all_selected_annotated_sei <- df_phenotype_binded_all_selected_annotated%>% filter(dMOCCS2score > 0 & q_value < 0.05)
  ID_snp_for_join <- df_phenotype_binded_all_selected_annotated_sei %>% 
    #filter(snp == target_position) %>% 
    filter(position == target_position) %>% 
    #unite("ID_position", c(ID, snp)) %>% 
    unite("ID_position", c(ID, position)) %>% 
    .$ID_position %>% 
    as.character()
  df_phenotype_binded_all_selected_annotated_sei2 <- df_phenotype_binded_all_selected_annotated_sei %>% 
    #filter(snp == target_position) %>%
    filter(position == target_position) %>% 
    mutate(ID_position = ID_snp_for_join) %>% 
    drop_na(dMOCCS2score) 
  df2 <- df_phenotype_binded_all_selected_annotated_sei2 %>% mutate(abs_dMOCCS2score = abs(dMOCCS2score)) %>% 
    #select(ID_position, snp, abs_dMOCCS2score, ID, Antigen, Cell_type_class, Cell_type) %>%
    select(ID_position, position, dMOCCS2score,abs_dMOCCS2score, ID, Antigen, Cell_type_class, Cell_type) %>%
    distinct() %>%
    #group_by(ID, Antigen, Cell_type_class, snp) %>% 
    group_by(ID, Antigen, Cell_type_class, position) %>% 
    summarise(max_dMOCCS2score= max(dMOCCS2score),min_MOCCS2score = min(dMOCCS2score),  max_abs_dMOCCS2score = max(abs_dMOCCS2score))
  
  #x_lab_ID <- df2  %>% filter(max_abs_dMOCCS2score > threshold) %>% arrange(desc(max_abs_dMOCCS2score)) %>% .$ID %>% as.character()
  #x_lab_ID <- df2  %>% unite("ID_position", c(ID, snp)) %>% filter(max_abs_dMOCCS2score > threshold) %>% arrange(desc(max_abs_dMOCCS2score)) %>% select(ID_position)
  x_lab_ID <- df2  %>% unite("ID_position", c(ID, position)) %>% filter(max_abs_dMOCCS2score > threshold) %>% arrange(desc(max_abs_dMOCCS2score)) %>% select(ID_position)
  x_lab_ID2 <- x_lab_ID %>% separate(ID_position, into = c("ID", "chr", "position"), sep = "_") %>% .$ID %>% as.character() %>% rev()
  
  Fig6E_TF_sei <- df2  %>% 
    #unite("ID_position", c(ID, snp)) %>% 
    unite("ID_position", c(ID, position)) %>%
    filter(max_abs_dMOCCS2score > threshold) %>%
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
    scale_x_discrete("ID", labels = x_lab_ID2) +
    ylab("max(|dMOCCS2score|) (dMOCCS2score > 0)") +
    ggtitle(paste0(target_rs, " ", target_phenotype)) +
    coord_flip()
  
  # dMOCCS2score <0 ------
  threshold <- 5
  df_phenotype_binded_all_selected_annotated_minus <- df_phenotype_binded_all_selected_annotated%>% filter(dMOCCS2score < 0 & q_value < 0.05)
  ID_snp_for_join <- df_phenotype_binded_all_selected_annotated_minus %>% 
    #filter(snp == target_position) %>% 
    filter(position == target_position) %>% 
    #unite("ID_position", c(ID, snp)) %>% 
    unite("ID_position", c(ID, position)) %>% 
    .$ID_position %>% 
    as.character()
  df_phenotype_binded_all_selected_annotated_minus2 <- df_phenotype_binded_all_selected_annotated_minus %>% 
    #filter(snp == target_position) %>%
    filter(position == target_position) %>% 
    mutate(ID_position = ID_snp_for_join) %>% 
    drop_na(dMOCCS2score) 
  df2 <- df_phenotype_binded_all_selected_annotated_minus2 %>% mutate(abs_dMOCCS2score = abs(dMOCCS2score)) %>% 
    #select(ID_position, snp, abs_dMOCCS2score, ID, Antigen, Cell_type_class, Cell_type) %>%
    select(ID_position, position, dMOCCS2score,abs_dMOCCS2score, ID, Antigen, Cell_type_class, Cell_type) %>%
    distinct() %>%
    #group_by(ID, Antigen, Cell_type_class, snp) %>% 
    group_by(ID, Antigen, Cell_type_class, position) %>% 
    summarise(max_dMOCCS2score= max(dMOCCS2score),min_MOCCS2score = min(dMOCCS2score),  max_abs_dMOCCS2score = max(abs_dMOCCS2score))
  
  x_lab_ID <- df2  %>% unite("ID_position", c(ID, position)) %>% filter(max_abs_dMOCCS2score > threshold) %>% arrange(desc(max_abs_dMOCCS2score)) %>% select(ID_position)
  x_lab_ID2 <- x_lab_ID %>% separate(ID_position, into = c("ID", "chr", "position"), sep = "_") %>% .$ID %>% as.character() %>% rev()
  
  Fig6E_TF_minus <- df2  %>% 
    #unite("ID_position", c(ID, snp)) %>% 
    unite("ID_position", c(ID, position)) %>%
    filter(max_abs_dMOCCS2score > threshold) %>%
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
    scale_x_discrete("ID", labels = x_lab_ID2) +
    ylab("max(|dMOCCS2score|) (dMOCCS2score < 0)") +
    ggtitle(paste0(target_rs, " ", target_phenotype)) +
    coord_flip()
  
  return(list(Fig6E_TF_sei, Fig6E_TF_minus))
  
}