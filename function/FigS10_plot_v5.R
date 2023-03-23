FigS10_plot_v5 <- function(target_phenotype, target_tf_list, FigS10B_plot_need){
  library(tidyverse)
  target_phenotype_list <- c("SLE", "MS", "IBD", "CD")
  print(target_phenotype)
  
  
  if(target_phenotype == "CD"){
    #target_df <- readRDS(url("https://figshare.com/ndownloader/files/39276440", "rb"))
    target_df <- readRDS("~/MOCCS_paper_public/data/Fig6/CD_binded_all.rds")
  }else if(target_phenotype == "MS"){
    #target_df <- readRDS(url("https://figshare.com/ndownloader/files/39276446", "rb"))
    target_df <- readRDS("~/MOCCS_paper_public/data/Fig6/MS_binded_all.rds")
  }else if(target_phenotype == "SLE"){
    #target_df <- readRDS(url("https://figshare.com/ndownloader/files/39276452", "rb"))
    target_df <- readRDS("~/MOCCS_paper_public/data/Fig6/SLE_binded_all.rds")
  }else if(target_phenotype == "IBD"){
    #target_df <- readRDS(url("https://figshare.com/ndownloader/files/39276443", "rb"))
    target_df <- readRDS("~/MOCCS_paper_public/data/Fig6/IBD_binded_all.rds")
  }
  
  #totalization <- readRDS(url("https://figshare.com/ndownloader/files/34065686","rb")) #MOCCSout_hg38_all_qval_annotated.rds
  totalization <- readRDS("~/MOCCS_paper_public/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds")
  annotation <- totalization %>% select(ID, Antigen, Cell_type_class, Cell_type) %>% distinct()
  
  df_phenotype_binded_all_selected_annotated <- target_df %>% left_join(annotation, by = "ID")
  ID_hard <- readRDS("~/MOCCS_paper_public/data/Fig1/hg38_hard_filter_ID.rds")
  TF_list  <-  annotation %>% filter(ID %in% ID_hard) %>% group_by(Antigen) %>% summarise(n= n()) %>% arrange(desc(n)) %>% filter(n >= 15) %>% distinct() %>% .$Antigen %>%  as.character()
  
  p_list_within  <- c()
  p_list_without <- c()
  for (target_tf in TF_list) {
    print(target_tf)
    
    # within peaks ---
    df1 <- df_phenotype_binded_all_selected_annotated %>% filter(Antigen == target_tf & peak == "within peaks") %>% drop_na(dMOCCS2score)
    df2 <- df1 %>% mutate(significant = ifelse(q_val < 0.05, "significant", "non-significant"))
    df3 <- df2 %>% unite("chr_posi", c(chr, chr_position), sep = "_") %>% group_by(ID, chr_posi) %>% summarise(max_abs_dMO = max(abs(dMOCCS2score)))
    join_key <- df3 %>% unite("key", c(ID, chr_posi), sep = "_") %>% .$key
    df3_2 <- df2 %>% unite("chr_posi", c(chr, chr_position), sep = "_") %>% unite("key", c(ID, chr_posi), sep = "_") %>% select(key, q_val) %>%
      group_by(key) %>% summarise(min_qval = min(q_val))
    df3_3 <- df3 %>% ungroup() %>% mutate(key = join_key) %>% left_join(df3_2, by = "key")
    df3_4 <- df3_3 %>% mutate(significant = ifelse(min_qval < 0.05, "significant", "non-significant")) %>% 
      group_by(ID, significant) %>% summarise(snp_num = n()) 
    
    df4 <- df3_4 %>% pivot_wider(names_from = "significant", values_from = "snp_num")
    
    if(ncol(df4)  == 3){
      colnames(df4) <- c("ID","significant", "non_significant")
      df5 <- df4 %>% pivot_longer(-ID, names_to = "significant", values_to = "snp_num")
      
      p1 <- df5 %>% drop_na(snp_num) %>% ggplot(aes(x = reorder(ID, -snp_num),  y = snp_num, fill = significant)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c("#696969", "#DC143C"), labels = c(ratio_nonsig = "non significant", ratio_sig ="significant")) +
        ggtitle(target_tf) +
        ylab("Number of peak-overlapping snps") +
        theme(plot.title = element_text(face="bold",hjust = 0.5), 
              panel.grid.major = element_line(colour = "gray"),
              panel.grid.minor = element_line(colour="gray"),
              panel.background = element_blank(), 
              axis.line = element_line(colour="black"),
              axis.text=element_text(size=5,face="bold"),
              #axis.text.x =element_text(size=5,face="bold", angle = 45, hjust = 1),
              axis.title.x = element_blank(),
              axis.text.y =element_text(size=3,face="bold"),
              axis.title=element_text(size=5,face="bold"),
              legend.position = 'none',
              legend.title = element_blank(),
              aspect.ratio = 1
        )+
        labs(fill = "") 
      
      #plot(p1)
    }
    
    
    #  without peaks ----
    df6 <- df_phenotype_binded_all_selected_annotated %>% filter(Antigen == target_tf & peak == "without peaks") %>% drop_na(dMOCCS2score)
    df7 <- df6 %>% mutate(significant = ifelse(q_val < 0.05, "significant", "non-significant"))
    df8 <- df7 %>% unite("chr_posi", c(chr, chr_position), sep = "_") %>% group_by(ID, chr_posi) %>% summarise(max_abs_dMO = max(abs(dMOCCS2score)))
    join_key <- df8 %>% unite("key", c(ID, chr_posi), sep = "_") %>% .$key
    df8_2 <- df7 %>% unite("chr_posi", c(chr, chr_position), sep = "_") %>% unite("key", c(ID, chr_posi), sep = "_") %>% select(key, q_val) %>%
      group_by(key) %>% summarise(min_qval = min(q_val))
    df8_3 <- df8 %>% ungroup() %>% mutate(key = join_key) %>% left_join(df8_2, by = "key")
    df8_4 <- df8_3 %>% mutate(significant = ifelse(min_qval < 0.05, "significant", "non-significant")) %>% 
      group_by(ID, significant) %>% summarise(snp_num = n()) 
    
    df9 <- df8_4 %>% pivot_wider(names_from = "significant", values_from = "snp_num")
    
    
    if(ncol(df9)  == 3){
      colnames(df9) <- c("ID","significant", "non_significant")
      df10 <- df9 %>% pivot_longer(-ID, names_to = "significant", values_to = "snp_num")
      
      p2<- df10 %>% drop_na(snp_num) %>% ggplot(aes(x = reorder(ID, -snp_num),  y = snp_num, fill = significant)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c("#696969", "#DC143C"), labels = c(ratio_nonsig = "non significant", ratio_sig ="significant")) +
        ggtitle(target_tf) +
        ylab("Number of no-peak-overlapping snps") +
        theme(plot.title = element_text(face="bold",hjust = 0.5), 
              panel.grid.major = element_line(colour = "gray"),
              panel.grid.minor = element_line(colour="gray"),
              panel.background = element_blank(), 
              axis.line = element_line(colour="black"),
              axis.text=element_text(size=5,face="bold"),
              #axis.text.x =element_text(size=5,face="bold", angle = 45, hjust = 1),
              axis.title.x = element_blank(),
              axis.text.y =element_text(size=3,face="bold"),
              axis.title=element_text(size=5,face="bold"),
              legend.position = 'none',
              legend.title = element_blank(),
              aspect.ratio = 1
        )+
        labs(fill = "") 
      
      #plot(p2)
    }
    
    
    if(target_tf == TF_list[1]){
      p_list_within <- p1
      p_list_without <- p2
    }else{
      p_list_within <- p_list_within + p1
      p_list_without <- p_list_without + p2
    }
  }
  
  p_list_within
  p_list_without
}