FigS7_plot <- function(path, target_tf, annotation_path){
  
  
  library(tidyverse)
  
  # difference among position of SNPs in k-mer
  target_df <- read_tsv(paste0(annotation_path, "SNP_SELEX/dMOCCS2score_hg38_", target_tf, "_all.tsv"))
  totalization <- readRDS("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds")
  annotation <- totalization %>% select(ID, Antigen, Cell_type_class, Cell_type) %>% filter(Antigen == target_tf) %>% distinct()
  target_df <- target_df %>% left_join(annotation, by = "ID")
  
  p_list <- list()
  for (i in 1:6) {
    target_df2 <- target_df %>% filter(posi == i)
    cor <- cor(target_df2$pbs, target_df2$dMOCCS2score, method = "spearman", use = "complete.obs")
    cor <- format(cor, digits = 3)
    
    p_list[[i]] <- target_df2 %>% ggplot(aes(x = dMOCCS2score, y = pbs))+
      geom_point(size = 1)+
      #scale_fill_viridis_c(breaks= c(100, 500, 1500, 2500, 4000))+
      geom_smooth(se = FALSE, method = lm) +
      #ggtitle(target_tf)+
      xlab("") +
      ylab("") +
      annotate("text", x=-Inf,y=Inf, label= paste0("cor = ", cor),hjust=-.2,vjust=3,fontface="bold", colour="blue",size=7) +
      theme(plot.title = element_text(face="bold",hjust = 0.5), 
            panel.grid.major = element_line(colour = "gray"),
            panel.grid.minor =  element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour="black"),
            axis.text=element_text(size=12,face="bold"),
            axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
            axis.text.y =element_text(size=10,face="bold"),
            axis.title=element_text(size=14,face="bold"),
            aspect.ratio = 1
      )
    
  }
  
  
  original_df <- read_tsv(paste0(path, "ADASTRA_data/" , target_tf, "_HUMAN.tsv"))
  pwm_df <- original_df %>% select('#chr', pos, motif_log_pref, motif_log_palt, motif_fc, motif_pos) %>% unite("snp", c('#chr', pos))
  
  # dMOCCS2scoreのデータと結合
  df <- readRDS(paste0(path, "dMOCCS2score/ADASTRA_dMOCCS2score_", target_tf, "_all.rds"))
  df2 <- df %>% mutate(pos = end + posi - 6) %>% unite("snp", c(chr, pos))
  df3 <- left_join(df2, pwm_df, by = "snp")
  df4 <- df3 %>% #mutate(color = ifelse((q_value < 0.05 & (motif_fc > 2 | motif_fc < -2)), "differential", "non-differential"))
    mutate(color = ifelse((q_value < 0.05 & motif_fc > 0 & dMOCCS2score < 0) | (q_value < 0.05 & motif_fc < 0 & dMOCCS2score > 0), "concordant", "disconcordant"))
  
  p1 <- df4 %>%
    ggplot(aes(x = dMOCCS2score, y = motif_fc)) +
    geom_bin_2d(bins = 100) +
    scale_fill_viridis_c(breaks= c(100, 500, 1500, 2500, 4000))+
    #geom_point(size = 0.1, alpha = 0.2) +
    #scale_color_manual(values = c("#DC143C", "#696969")) +
    geom_smooth(se = FALSE, method = lm) +
    ylab("motif fc") +
    xlab("dMOCCS2score") +
    ggtitle(target_tf)+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
          axis.text.y =element_text(size=10,face="bold"),
          axis.title=element_text(size=14,face="bold")
    )
  
  
  
  return(list(p_list, p1))
  
}