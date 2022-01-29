Fig1G_plot <- function(target_TF, annotation_path){
  
  library(tidyverse)
  library(RColorBrewer)
  
  totalization_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds"
  totalization <- readRDS(totalization_path)
  ID_hard <- readRDS(paste0(annotation_path, "hg38_hard_filter_ID.rds"))
  
  # filter target TF PWM table
  PWM_table_all <- readRDS("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/PWM_likelihood_HOMER.rds")
  target_PWM <- PWM_table_all[[target_TF]]
  
  # join MOCCS output and PWM table
  target_MOCCS <- totalization %>% filter(ID %in% ID_hard) %>% filter(Antigen == target_TF)
  df_join1 <- target_MOCCS %>% left_join(target_PWM, by = "kmer")
  df_join1[is.na(df_join1)] <- 0 #NAを0に置換
    
  # PWM scoreをsignificant k-merとそれ以外で分けてplot
  df_join2 <- df_join1 %>% mutate(sig_kmer = ifelse(q_value < 0.05, "significant k-mer", "non significant k-mer")) 
  df_join3 <- transform(df_join2, sig_kmer = factor(sig_kmer, levels = c("non significant k-mer", "significant k-mer")))
  p <- df_join3 %>% ggplot(aes(x=sig_kmer, y = PWMscore, color = sig_kmer)) +
    geom_violin()  +
    xlab("k-mer")+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_line(colour="gray"),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
          axis.text.y =element_text(size=10,face="bold"),
          axis.title=element_text(size=14,face="bold")
    )
  
  
  
  
}