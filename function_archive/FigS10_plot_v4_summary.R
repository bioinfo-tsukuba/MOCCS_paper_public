library(tidyverse)
library(patchwork)

if(target_phenotype == "CD"){
  #target_df <- readRDS(url("https://figshare.com/ndownloader/files/34065827", "rb"))
  target_df <- readRDS(url("https://figshare.com/ndownloader/files/39276440", "rb"))
  #target_df <- readRDS("~/Downloads/CD_binded_all.rds")
}else if(target_phenotype == "MS"){
  #target_df <- readRDS(url("https://figshare.com/ndownloader/files/34065830", "rb"))
  target_df <- readRDS(url("https://figshare.com/ndownloader/files/39276446", "rb"))
}else if(target_phenotype == "SLE"){
  #target_df <- readRDS(url("https://figshare.com/ndownloader/files/34660126", "rb"))
  target_df <- readRDS(url("https://figshare.com/ndownloader/files/39276452", "rb"))
}else if(target_phenotype == "IBD"){
  #target_df <- readRDS(url("https://figshare.com/ndownloader/files/34660132", "rb"))
  target_df <- readRDS(url("https://figshare.com/ndownloader/files/39276443", "rb"))
}

totalization <- readRDS(url("https://figshare.com/ndownloader/files/34065686","rb")) #MOCCSout_hg38_all_qval_annotated.rds
#totalization <- readRDS("~/MOCCS_paper_public/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds")
annotation <- totalization %>% select(ID, Antigen, Cell_type_class, Cell_type) %>% distinct()
df_phenotype_binded_all_selected_annotated <- target_df %>% left_join(annotation, by = "ID")

ID_hard <- readRDS("~/MOCCS_paper_public/data/Fig1/hg38_hard_filter_ID.rds")
TF_list  <-  annotation %>% filter(ID %in% ID_hard) %>% group_by(Antigen) %>% summarise(n= n()) %>% arrange(desc(n)) %>% filter(n >= 15) %>% distinct() %>% .$Antigen %>%  as.character()
length(TF_list)


# plot ---
p_list_within  <- c()
p_list_without <- c()
for (target_tf in TF_list) {
  print(target_tf)
  
  # within peaks ---
  df1 <- df_phenotype_binded_all_selected_annotated %>% filter(Antigen == target_tf & peak == "within peaks") %>% drop_na(dMOCCS2score)
  if(nrow(df1) != 0){
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
  
}

p_list_within
p_list_without
#plot(p_list_within + plot_layout(nrow = 6, ncol = 6), height = 50)
#plot(p_list_without + plot_layout(nrow = 6, ncol = 6), height = 50)

# make summary table ----

summary_tib <- tibble()
for (target_tf in TF_list) {
  print(target_tf)
  # within peaks ---
  df1 <- df_phenotype_binded_all_selected_annotated %>% filter(Antigen == target_tf & peak == "within peaks") %>% drop_na(dMOCCS2score)
  if(nrow(df1) != 0){
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
    }
    
    df11 <- df5 %>% mutate(label = "within")
    df12 <- df10 %>% mutate(label = "without")
    df13 <- rbind(df11, df12)
    df13 <- df13 %>% as_tibble() %>% mutate(TF = target_tf)
    
    if(nrow(summary_tib) == 0){
      summary_tib <- df13
    }else{
      summary_tib <- rbind(summary_tib, df13)
    }
  }
  
}

p <- summary_tib %>% ggplot(aes(x = label, y = snp_num, fill = significant)) +
  geom_violin()
p

summary_tib2 <- tibble()
for (tgt_ID in unique(summary_tib$ID)) {
  tmp_within <- tibble()
  tmp_without <- tibble()
  tgt_within <- summary_tib  %>% filter(ID == tgt_ID &label == "within") %>% drop_na(snp_num)
  tgt_without <- summary_tib  %>% filter(ID == tgt_ID &label == "without") %>% drop_na(snp_num)
  ratio_within <- tgt_within$snp_num/sum(tgt_within$snp_num)
  ratio_without <- tgt_without$snp_num/sum(tgt_without$snp_num)
  tmp_within <- tgt_within %>% mutate(ratio = ratio_within)
  tmp_without <- tgt_without %>% mutate(ratio = ratio_without)
  if(nrow(tmp_within) != 0 & nrow(tmp_without) != 0){
    tmp_df <- rbind(tmp_within, tmp_within)
  }else if(nrow(tmp_within) == 0){
    tmp_df <- tmp_without
  }else{
    tmp_df <- tmp_within
  }
  if(tgt_ID == unique(summary_tib$ID)[1]){
    summary_tib2 <- tmp_df
  }else{
    summary_tib2 <- rbind(summary_tib2, tmp_df)
  }
}
write_tsv(summary_tib2, paste0("~/MOCCS_paper_public/data/Fig6/sig_snp_ratio_", target_phenotype, ".tsv"))

summary_tib3 <- summary_tib2 %>% mutate(label2 = ifelse(label == "within", "within peaks", "out of peaks"))
summary_tib3 %>% filter(significant == "significant") %>%
  ggplot(aes(x = TF, y = ratio, fill = label2)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#DC143C", "#696969")) +
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_line(colour="gray"),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        #legend.position = 'none',
        legend.title = element_blank(),
        aspect.ratio = 0.3
  )+
  ggtitle(target_phenotype) +
  ylab("Ratio of SNPs with significant dMOCCS2score") +
  labs(fill = c("within peaks", "out of peaks")) 



# 
patch <- c()
for (tgt_tf in as.character(unique(summary_tib3$Antigen))) {
  tmp <- summary_tib3 %>% filter(Antigen == tgt_tf)
  if(length(unique(tmp$label)) != 1){

    p <- summary_tib3 %>% filter(Antigen == tgt_tf) %>%
      ggplot(aes(x = label, y = ratio, fill = label)) +
      #geom_violin() +
      geom_boxplot()+
      geom_point() +
      ggtitle(tgt_tf) +
      ylab("Ratio of significant SNPs") +
      scale_fill_manual(values = c("#DC143C", "#696969")) +
      theme(plot.title = element_text(face="bold",hjust = 0.5),
            panel.grid.major = element_line(colour = "gray"),
            panel.grid.minor = element_line(colour="gray"),
            panel.background = element_blank(),
            axis.line = element_line(colour="black"),
            axis.text=element_text(size=5,face="bold"),
            axis.text.x =element_text(size=5,face="bold", angle = 45, hjust = 1),
            axis.text.y =element_text(size=5,face="bold"),
            axis.title=element_text(size=5,face="bold"),
            legend.position = 'none',
            legend.title = element_blank(),
            aspect.ratio = 1.2
      )+
      labs(fill = "")
    plot(p)
    if(length(patch) == 0){
      patch <- p
    }else{
      patch <- patch + p
    }
  }
}
patch
