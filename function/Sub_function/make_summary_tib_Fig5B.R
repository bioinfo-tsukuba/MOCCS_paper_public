library(tidyverse)
library(exactRankTests)
target_phenotype_list <- c("SLE", "MS", "IBD", "CD")
#totalization <- readRDS(url("https://figshare.com/ndownloader/files/34065686","rb")) #MOCCSout_hg38_all_qval_annotated.rds
totalization <- readRDS("~/MOCCS_paper_public/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds")
annotation <- totalization %>% select(ID, Antigen, Cell_type_class, Cell_type) %>% distinct()


for (target_phenotype in target_phenotype_list) {
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
  
  df_phenotype_binded_all_selected_annotated <- target_df %>% left_join(annotation, by = "ID")
  ID_hard <- readRDS("~/MOCCS_paper_public/data/Fig1/hg38_hard_filter_ID.rds")
  #TF_list  <-  annotation %>% filter(ID %in% ID_hard) %>% group_by(Antigen) %>% summarise(n= n()) %>% arrange(desc(n)) %>% filter(n >= 15) %>% distinct() %>% .$Antigen %>%  as.character()
  TF_list  <-  annotation %>% filter(ID %in% ID_hard) %>% group_by(Antigen) %>% summarise(n= n()) %>% arrange(desc(n)) %>% filter(n >= 12) %>% distinct() %>% .$Antigen %>%  as.character()
  
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
        df4[is.na(df4)] <- 0
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
        df9[is.na(df9)] <- 0
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
      tmp_df <- rbind(tmp_within, tmp_without)
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
  #write_tsv(summary_tib2, paste0("~/MOCCS_paper_public/data/Fig6/sig_snp_ratio_", target_phenotype, ".tsv"))
  write_tsv(summary_tib2, paste0("~/MOCCS_paper_public/data/Fig6/sig_snp_ratio_", target_phenotype, "_TF42.tsv"))
}

# wilcoxon test per TFs
q_list_all <- list()
p_list <- list()
for (target_phenotype in target_phenotype_list) {
  print(target_phenotype)
  #summary_tib2 <- read_tsv(paste0("~/MOCCS_paper_public/data/Fig6/sig_snp_ratio_", target_phenotype, ".tsv"))
  summary_tib2 <- read_tsv(paste0("~/MOCCS_paper_public/data/Fig6/sig_snp_ratio_", target_phenotype, "_TF42.tsv"))
  for (target_tf  in TF_list) {
    print(target_tf)
    tgt_df <- summary_tib2 %>% filter(TF == target_tf)
    within_df <- tgt_df %>% filter(label == "within") %>% filter(significant == "significant") 
    out_of_df <- tgt_df %>% filter(label == "without") %>% filter(significant == "significant")
    share_ID <- intersect(unique(within_df$ID), unique(out_of_df$ID))
    if(length(share_ID) != 0){
      within_ratio <- within_df %>% filter(ID %in% share_ID) %>% .$ratio %>% as.numeric()
      out_of_ratio <- out_of_df %>% filter(ID %in% share_ID) %>% .$ratio %>% as.numeric()
      
      res <- wilcox.exact(within_ratio, out_of_ratio, paired = T, alternative = "greater")
      p <- res$p.value %>% as.numeric()
      p_list[[target_tf]] <- p
    }
  }
  q_list_all[[target_phenotype]] <- p.adjust(p_list)
}

df_qval <- tibble()
for (target_phenotype in target_phenotype_list) {
  print(target_phenotype)
  df1 <- tibble(TF = names(q_list_all[[target_phenotype]]), qval = q_list_all[[target_phenotype]])
  df2 <- df1 %>% mutate(phenotype = target_phenotype)
  if(nrow(df_qval) == 0){
    df_qval <- df2
  }else{
    df_qval <- rbind(df_qval, df2)
  }
}

df_qval2 <- df_qval %>% mutate(significance_wilcoxon = ifelse(qval < 0.05, "TRUE", "FALSE"))
#write_tsv(df_qval2, "~/MOCCS_paper_public/data/Fig6/sig_ratio_wilcoxon_exact.tsv")
write_tsv(df_qval2, "~/MOCCS_paper_public/data/Fig6/sig_ratio_wilcoxon_exact_TF42.tsv")


for (target_phenotype in target_phenotype_list) {
  print(target_phenotype)
  summary_tib2 <- read_tsv(paste0("~/MOCCS_paper_public/data/Fig6/sig_snp_ratio_", target_phenotype, ".tsv"))
  summary_tib3 <- summary_tib2 %>% mutate(label2 = ifelse(label == "within", "within peaks", "out of peaks"))
  p <- summary_tib3 %>% filter(significant == "significant") %>%
    ggplot(aes(x = TF, y = ratio, fill = label2)) +
    geom_boxplot() +
    scale_fill_manual(values = c( "#696969", "#DC143C")) +
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
    ylab("Ratio of SNPs with significant dMOCCS2score") 
  plot(p)
}
