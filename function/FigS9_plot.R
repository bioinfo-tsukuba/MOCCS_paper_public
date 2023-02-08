FigS9_plot <- function(path, target_TF){
  
  
  library(tidyverse)
  
  # Fig. S9A
  df_pbs <- read_csv("~/MOCCS_paper_public/data/supplement/FigS9_pbs_table.csv")
  tf_order <- df_pbs %>% filter(label == "original") %>% arrange(desc(sp_cor)) %>% .$Antigen 
  p_allTF <- df_pbs %>% ggplot(aes(x = reorder(Antigen, sp_cor), y = sp_cor, fill = label, color = label))+
    geom_point(size = 5)+
    scale_x_discrete(limit = tf_order)+
    xlab("") +
    ylab("Spearman correlation coefficient")+
    scale_color_manual(values = c("red", "blue"))+
    scale_fill_manual(values = c("red", "blue"))+
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
  p_allTF 
  
  tf_order2 <- df_pbs %>% filter(label == "original") %>% arrange(desc(pea_cor)) %>% .$Antigen 
  p_allTF2 <- df_pbs %>% ggplot(aes(x = reorder(Antigen, pea_cor), y = pea_cor, fill = label, color = label))+
    geom_point(size = 5)+
    scale_x_discrete(limit = tf_order2)+
    xlab("") +
    ylab("Pearson correlation coefficient")+
    scale_color_manual(values = c("red", "blue"))+
    scale_fill_manual(values = c("red", "blue"))+
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
  p_allTF2
  
  # difference among position of SNPs in k-mer (Fig. S9B)
  if(target_TF == "GATA3"){
    target_df <- read_tsv(url("https://figshare.com/ndownloader/files/34336355", "rb"))
  }else if(target_TF == "FOXA1"){
    target_df <- read_tsv(url("https://figshare.com/ndownloader/files/34336352", "rb"))
  }
  totalization <- readRDS(url("https://figshare.com/ndownloader/files/34065686","rb")) #MOCCSout_hg38_all_qval_annotated.rds
  annotation <- totalization %>% select(ID, Antigen, Cell_type_class, Cell_type) %>% filter(Antigen == target_TF) %>% distinct()
  target_df <- target_df %>% left_join(annotation, by = "ID")
  
  p_list <- list()
  for (i in 1:6) {
    target_df2 <- target_df %>% filter(posi == i)
    cor <- cor(target_df2$pbs, target_df2$dMOCCS2score, method = "spearman", use = "complete.obs")
    cor <- format(cor, digits = 3)
    
    p_list[[i]] <- target_df2 %>% ggplot(aes(x = dMOCCS2score, y = pbs))+
      geom_point(size = 1)+
      geom_smooth(se = FALSE, method = lm) +
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
  
  if(target_TF == "GATA3"){
    original_df <- read_tsv(url("https://figshare.com/ndownloader/files/34337135", "rb"))
  }else if(target_TF == "FOXA1"){
    original_df <- read_tsv(url("https://figshare.com/ndownloader/files/34337141", "rb"))
  }
  pwm_df <- original_df %>% select('#chr', pos, motif_log_pref, motif_log_palt, motif_fc, motif_pos) %>% unite("snp", c('#chr', pos))
  
  # dMOCCS2scoreのデータと結合
  #df <- readRDS(paste0(path, "dMOCCS2score/ADASTRA_dMOCCS2score_", target_tf, "_all.rds"))
  url_list <- read_csv("~/MOCCS_paper_public/data/Fig5/Fig5E_url_list.csv", col_names = FALSE)
  colnames(url_list) <- c("TF", "url")
  target_url <- url_list %>% filter(TF == target_TF) %>% .$url %>% as.character()
  df <- readRDS(url(target_url, "rb"))
  df2 <- df %>% mutate(pos = end + posi - 6) %>% unite("snp", c(chr, pos))
  df3 <- left_join(df2, pwm_df, by = "snp")
  df4 <- df3 %>% #mutate(color = ifelse((q_value < 0.05 & (motif_fc > 2 | motif_fc < -2)), "differential", "non-differential"))
    mutate(color = ifelse((q_value < 0.05 & motif_fc > 0 & dMOCCS2score < 0) | (q_value < 0.05 & motif_fc < 0 & dMOCCS2score > 0), "concordant", "disconcordant"))
  
  # Fig. S9C
  p1 <- df4 %>%
    ggplot(aes(x = dMOCCS2score, y = motif_fc)) +
    geom_bin_2d(bins = 100) +
    scale_fill_viridis_c(breaks= c(100, 500, 1500, 2500, 4000))+
    geom_smooth(se = FALSE, method = lm) +
    ylab("motif fc") +
    xlab("dMOCCS2score") +
    ggtitle(target_TF)+
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
  
  
  
  return(list(p_list, p1, p_allTF, p_allTF2))
  
}