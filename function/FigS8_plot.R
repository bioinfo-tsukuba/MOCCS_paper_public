FigS8_plot <- function(path, target_phenotype, target_tf){
  
  library(tidyverse)
  target_df <- readRDS(paste0(path, "result_output_binded_all/", target_phenotype, "_peak_rand_binded_all_qval.rds"))
  
  totalization_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds"
  totalization <- readRDS(totalization_path)
  annotation <- totalization %>% select(ID, Antigen, Cell_type_class, Cell_type) %>% distinct()
  df_phenotype_binded_all_selected_annotated <- target_df %>% left_join(annotation, by = "ID")
  
  # 1.number of peak-overlapping snps (target tf)
  df1 <- df_phenotype_binded_all_selected_annotated %>% filter(Antigen == target_tf) %>% drop_na(dMOCCS2score)
  df2 <- df1 %>% mutate(significant = ifelse(q_value < 0.05, "significant", "non-significant"))
  df3 <- df2 %>% group_by(ID, significant) %>% summarise(snp_num = n()) 
  
  df4 <- df3 %>% pivot_wider(names_from = "significant", values_from = "snp_num")
  colnames(df4) <- c("ID", "non_significant", "significant")
  df5 <- df4 %>% pivot_longer(-ID, names_to = "significant", values_to = "snp_num")
  
  p1 <- df5 %>%  ggplot(aes(x = reorder(ID, -snp_num),  y = snp_num, fill = significant)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("#696969", "#DC143C"), labels = c(ratio_nonsig = "non significant", ratio_sig ="significant")) +
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
          legend.title = element_blank(),
          aspect.ratio = 1
    )+
    ggtitle(target_tf) +
    ylab("Number of peak-overlapping snps") +
    labs(fill = "") 
  
  
  # 2. all TFs number of snps (mean per sample)
  tf_list <- c("CTCF", "FOXA1", "BRD4", "ESR1", "SPI1", "EP300", "JUND", "MYCN", "MAX", "EZH2", "SUMO2", "CEBPB", "GATA2", "GATA3", "HDAC2", "FOS")
  df_all <- tibble()
  for(i in 1:length(tf_list)){
    target_tf <- tf_list[i]
    df11 <- df_phenotype_binded_all_selected_annotated %>% filter(Antigen == target_tf) %>% drop_na(dMOCCS2score)
    df12 <- df11 %>% mutate(significant = ifelse(q_value < 0.05, "significant", "non-significant"))
    df13 <- df12 %>% group_by(ID, significant) %>% summarise(snp_num = n()) 
    df14 <- df13 %>% pivot_wider(names_from = "significant", values_from = "snp_num")
    colnames(df14) <- c("ID", "non_significant", "significant")
    df14 <- df14 %>% drop_na(non_significant)  %>% drop_na(significant)
    
    nonsig_mean_num <- mean(as.numeric(df14$non_significant))
    sig_mean_num <- mean(as.numeric(df14$significant))
    
    df_part <- tibble(tf = target_tf, nonsig_mean_num = nonsig_mean_num, sig_mean_num = sig_mean_num)
    
    if(nrow(df_all) == 0){
      df_all <- df_part
    }else{
      df_all <- df_all %>% add_row(df_part)
    }
  }
  
  
  df_all2 <- df_all %>% pivot_longer(-tf, names_to = "significant", values_to = "snp_num_mean")
  p2 <- df_all2 %>% ggplot(aes(x = reorder(tf, -snp_num_mean),  y = snp_num_mean, fill = significant))+
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("#696969", "#DC143C"), labels = c(nonsig_mean_num = "non significant", sig_mean_num ="significant")) +
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
          legend.title = element_blank(),
          aspect.ratio = 0.7
    )+
    ggtitle("Number of peak-overlapping snps (mean)") +
    ylab("Number of peak-overlapping snps") +
    labs(fill = "") 
  
  
  
  return(list(p1, p2))
  
}