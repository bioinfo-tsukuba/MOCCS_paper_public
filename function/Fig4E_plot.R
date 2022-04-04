Fig4E_plot <- function(target_ID1 = "", 
                       target_ID2 = "", 
                       target_CT = ""){
  
  library(tidyverse)
  library(ggrepel)
  totalization <- readRDS("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/MOCCSout_hg38_hard_filter_annotated.rds")
  
  target_MOCCS1 <- totalization %>% filter(ID == target_ID1) %>% select(kmer, auc, count, MOCCS2score, Antigen, Cell_type_class, Cell_type)
  colnames(target_MOCCS1) <- c("kmer", "auc1", "count1", "MOCCS2score1", "Antigen1", "Cell_type_class1", "Cell_type1")
  TF1 <- target_MOCCS1$Antigen1 %>% unique() %>% as.character()
  CTC1 <- target_MOCCS1$Cell_type_class1 %>% unique() %>% as.character()
  CT1 <- target_MOCCS1$Cell_type1 %>% unique() %>% as.character()
  
  target_MOCCS2 <- totalization %>% filter(ID == target_ID2) %>% select(kmer, auc, count, MOCCS2score, Antigen, Cell_type_class, Cell_type)
  colnames(target_MOCCS2) <- c("kmer", "auc2", "count2", "MOCCS2score2", "Antigen2", "Cell_type_class2", "Cell_type2")
  TF2 <- target_MOCCS2$Antigen2 %>% unique() %>% as.character()
  CTC2 <- target_MOCCS2$Cell_type_class2 %>% unique() %>% as.character()
  CT2 <- target_MOCCS2$Cell_type2 %>% unique() %>% as.character()
  
  df_join <- target_MOCCS1 %>% left_join(target_MOCCS2, by = "kmer") %>% drop_na(auc1, auc2)
  
  # calculate dMOCCS2score and p value
  W <- 350
  p_list_2 <- c()
  p_list_sig_1 <- c()
  p_list_sig_2 <- c()
  if(nrow(df_join) != 0){
    for(z in 1:nrow(df_join)){
      
      count_i <- df_join[z,"count1"] %>% as.numeric()
      count_j <- df_join[z,"count2"] %>% as.numeric()
      
      #AUC
      auc_i <- df_join[z,"auc1"] %>% as.numeric()
      auc_j <- df_join[z,"auc2"] %>% as.numeric()
      
      var_i <- W^2/12/count_i
      var_j <- W^2/12/count_j
      
      # MOCCS2score
      moccs2score_i <- auc_i / sqrt(var_i)
      moccs2score_j <- auc_j / sqrt(var_j)
      
      # sample内のpvalue (significant k-mer)
      target_p_i <- 1-pnorm(auc_i, mean = 0, sd = sqrt(W^2/12/count_i))
      target_p_j <- 1-pnorm(auc_j, mean = 0, sd = sqrt(W^2/12/count_j))
      
      p_list_sig_1 <- c(p_list_sig_1, target_p_i)
      p_list_sig_2 <- c(p_list_sig_2, target_p_j)
      
      # 2sample間のpvalue (differential k-mer)
      if (auc_i >= auc_j){
        p <- 1 - pnorm((auc_i-auc_j)/sqrt(var_i + var_j) , 0, 1) 
        p_list_2 <- c(p_list_2, p)
      }else{
        p <- 1 - pnorm((auc_j-auc_i)/sqrt(var_j + var_i) , 0, 1) 
        p_list_2 <- c(p_list_2, p) 
      }
    } 
  }
  
  df_join2 <- df_join %>% mutate(p_value = p_list_2, p_list_sig_1 = p_list_sig_1 , p_list_sig_2 = p_list_sig_2)
  q_list_2 <- p.adjust(p_list_2)
  q_list_sig_1 <- p.adjust(p_list_sig_1)
  q_list_sig_2 <- p.adjust(p_list_sig_2)
  df_join3 <- df_join2 %>% mutate(q_value = q_list_2, q_list_sig_1 = q_list_sig_1, q_list_sig_2 = q_list_sig_2) 
  df_join4 <- df_join3 %>% mutate(color = ifelse(q_value < 0.05, "differential", "non differential")) %>% distinct()
  if(target_ID1 == "SRX150546" & target_ID2 == "SRX190276"){
    df_join5 <- df_join4 %>% mutate(color4 = ifelse((str_detect(kmer, "CCCCC") | str_detect(kmer, "GGGGG")), "CTCF PWM motif", "non significant"))
    df_join6 <- df_join5 %>% mutate(color5 = ifelse((kmer %in%  "TGACTC" | kmer %in%  "GACTCA" | kmer %in% "TGAGTC" | kmer %in% "GAGTCA"), "JUN PWM motif", "non significant"))
    
  }else{
    df_join5 <- df_join4 %>% mutate(color4 = "")
    df_join6 <- df_join5 %>% mutate(color5 = "")
  }
  df_join7 <- df_join6 %>% unite("color6", c(color4, color5), sep = "_")
  df_join7$color6 <- gsub("non significant", "", df_join7$color6)
  df_join7$color6 <- gsub("_", "", df_join7$color6)
  df_join8 <- df_join7 %>% mutate(kmer_lab = ifelse((color6 == "significantbothCTCF PWM motif" | color6 == "significantbothJUN PWM motif"), kmer, "")) %>% distinct()
  
  selected_kmer <- df_join8 %>% filter(color6 == "CTCF PWM motif" | color6 == "JUN PWM motif") %>% .$kmer %>% unique()
  df_join9 <- df_join8 %>% mutate(kmer2 = df_join8$kmer)
  df_join9$kmer2[!df_join9$kmer2 %in% selected_kmer] <- ""
  
  p2 <- df_join9 %>% ggplot(aes(x = MOCCS2score1, y = MOCCS2score2, color = color, label = kmer2)) +
    geom_point(size = 0.8, alpha = 0.8) +
    #geom_text(aes(x =  MOCCS2score1 + 5, y = MOCCS2score2 + 5), size = 3) +
    ggtitle(target_CT) +
    geom_text_repel(size = 7) +
    xlab(paste0(target_ID1, "  (", TF1, ")")) +
    ylab(paste0(target_ID2, "  (", TF2, ")")) +
    scale_color_manual(values = c("#ff0000", "#000080")) +
    #xlim(c(-10, 70)) +
    #ylim(c(-10, 70)) +
    labs(color="") +
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_line(colour="gray"),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
          axis.text.y =element_text(size=10,face="bold"),
          axis.title=element_text(size=14,face="bold"),
          legend.position = 'none',
          aspect.ratio = 1
    )
  
  
  return(p2)
}