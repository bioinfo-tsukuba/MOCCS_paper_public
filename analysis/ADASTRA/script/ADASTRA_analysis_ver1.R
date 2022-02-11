rm(list = ls())
target_tf <- "CTCF"
ID_list <- read_tsv(paste0("/home/s-tahara/singularity_ADASTRA/ID_list/", target_tf, ".txt"), col_names = F)
#df_all <- readRDS(paste0("/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/data/dMOCCS2score_output_all/ADASTRA_dMOCCS2score_", target_tf, "_all.rds"))
df_all <- readRDS(paste0("/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/data/dMOCCS2score_output_all/ADASTRA_dMOCCS2score_", target_tf, "_all_part4.rds")) #CTCF

# 1. number of IDs 
length(ID_list$X1)
length(unique(df_all$ID))

# 2. number of snps
df_all2 <- df_all %>% mutate(position = end - 6 + posi) 
snp <- df_all2 %>% unite("snp", c(chr, position)) %>% .$snp %>% unique() %>% as.character()
length(unique(snp))

# 3. ASB significance summary
df3 <- df_all2 %>% mutate(fdrp_bh_ref_log = -log10(fdrp_bh_ref),
                          fdrp_bh_alt_log = -log10(fdrp_bh_alt), 
                          ASB = ifelse(fdrp_bh_ref_log <= fdrp_bh_alt_log, fdrp_bh_alt_log, fdrp_bh_ref_log), 
                          x_axis = ifelse(fdrp_bh_ref_log <= fdrp_bh_alt_log, "positive", "negative"))
df4 <- df3 %>% filter(x_axis == "positive") %>% mutate(ASB_plot = ASB)
df5 <- df3 %>% filter(x_axis == "negative") %>% mutate(ASB_plot = -ASB)
df6 <- rbind(df4, df5)
summary(df6$ASB_plot)

# 4. dMOCCS2score summary
summary(df6$dMOCCS2score)

# 5. number of differential k-mer
total_kmer_pair <- df_all2  %>% select(ID, kmer_before_tou, kmer_after_tou, dMOCCS2score, q_value) %>% distinct()
dif_kmer_pair <- df_all2 %>% filter(q_value < 0.05) %>% select(ID, kmer_before_tou, kmer_after_tou, dMOCCS2score, q_value) %>% distinct()
nrow(total_kmer_pair)
nrow(dif_kmer_pair)


# 6. filter ASB significance and number of concordant and disconcordant snps (dMOCCS2score q<0.05, threshold q<0.05)
rm(df3, df4, df5, df_all, total_kmer_pair, dif_kmer_pair, snp)
threshold <- -log10(0.05)
df11 <- df6 %>% mutate(color1 = ifelse((ASB_plot < -threshold  & q_value < 0.05 & dMOCCS2score > 0) | (ASB_plot > threshold  & q_value < 0.05 & dMOCCS2score < 0), "concordant", "disconcordant")) 
df12 <- df11 %>% filter(color1 == "concordant") %>% mutate(color2 = color1) %>%
  select(-fdrp_bh_ref, -fdrp_bh_alt, -MOCCS2score_before, -MOCCS2score_after, -fdrp_bh_ref_log, -fdrp_bh_alt_log)
df13 <- df11 %>%  filter(color1 == "disconcordant") %>% mutate(color2 = ifelse((q_value < 0.05 & abs(ASB_plot) > threshold), "disconcordant", "nondifferent")) %>%
  select(-fdrp_bh_ref, -fdrp_bh_alt, -MOCCS2score_before, -MOCCS2score_after, -fdrp_bh_ref_log, -fdrp_bh_alt_log)
rm(df11)
df14 <- rbind(df12, df13)
rm(df12, df13)

concordant <- df14 %>% filter(color2 == "concordant")
concordant_num <- nrow(concordant)
concordant_num
rm(concordant_num)

disconcordant <- df14 %>% filter(color2 == "disconcordant") 
disconcordant_num <- nrow(disconcordant)
disconcordant_num
rm(disconcordant_num)

p1 <- df14 %>% ggplot(aes(x = ASB_plot, y = dMOCCS2score, color = color2)) +
  geom_point(size = 0.5, alpha = 0.7) +
  xlab("ASB significance") +
  ylab("dMOCCS2score") +
  scale_colour_manual(
    values = c(
      concordant = "red3",
      disconcordant = "blue3", 
      nondifferent = "gray"
    )
  )+
  ggtitle(paste0(target_tf, " ", "dMOCCS2score q<0.05, threshold q<0.05"  ))+
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
#p1
#ggsave(paste0("/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/plot/", target_tf, "_dMOCCS005_ASB005.pdf"), p1, width = 7, height = 7)
rm(df14)


# 7. filter ASB significance and number of concordant and disconcordant snps (dMOCCS2score q<0.05, threshold q<0.01)
rm(df3, df4, df5, df_all,df14, total_kmer_pair, dif_kmer_pair, snp)
threshold <- -log10(0.01)
df11 <- df6 %>% mutate(color1 = ifelse((ASB_plot < -threshold  & q_value < 0.05 & dMOCCS2score > 0) | (ASB_plot > threshold  & q_value < 0.05 & dMOCCS2score < 0), "concordant", "disconcordant")) 
df12 <- df11 %>% filter(color1 == "concordant") %>% mutate(color2 = color1)
df13 <- df11 %>%  filter(color1 == "disconcordant") %>% mutate(color2 = ifelse((q_value < 0.05 & abs(ASB_plot) > threshold), "disconcordant", "nondifferent")) 
rm(df11)
df14 <- rbind(df12, df13)
rm(df12, df13)

concordant <- df14 %>% filter(color2 == "concordant")
concordant_num <- nrow(concordant)
concordant_num
rm(concordant, concordant_num)

disconcordant <- df14 %>% filter(color2 == "disconcordant") 
disconcordant_num <- nrow(disconcordant)
disconcordant_num
rm(disconcordant, disconcordant_num)

p2 <- df14 %>% ggplot(aes(x = ASB_plot, y = dMOCCS2score, color = color2)) +
  geom_point(size = 0.5, alpha = 0.7) +
  xlab("ASB significance") +
  ylab("dMOCCS2score") +
  scale_colour_manual(
    values = c(
      concordant = "red3",
      disconcordant = "blue3", 
      nondifferent = "gray"
    )
  )+
  ggtitle(paste0(target_tf, " ", "all samples", " ", "dMOCCS2score q<0.05, threshold q<0.01"))+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_line(colour="gray"),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=10,face="bold"),
        axis.title=element_text(size=14)
  )
#p2
#ggsave(paste0("/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/plot/", target_tf, "_dMOCCS005_ASB001.pdf"), p2, width = 7, height = 7)
rm(df14)

# 8.filter ASB significance and number of concordant and disconcordant snps (dMOCCS2score q<0.01, threshold q<0.05)
threshold <- -log10(0.05)
df11 <- df6 %>% mutate(color1 = ifelse((ASB_plot < -threshold  & q_value < 0.01 & dMOCCS2score > 0) | (ASB_plot > threshold  & q_value < 0.01 & dMOCCS2score < 0), "concordant", "disconcordant")) 
df12 <- df11 %>% filter(color1 == "concordant") %>% mutate(color2 = color1)
df13 <- df11 %>%  filter(color1 == "disconcordant") %>% mutate(color2 = ifelse((q_value < 0.01 & abs(ASB_plot) > threshold), "disconcordant", "nondifferent")) 
rm(df11)
df14 <- rbind(df12, df13)
rm(df12,df13)

concordant <- df14 %>% filter(color2 == "concordant")
concordant_num <- nrow(concordant)
concordant_num

disconcordant <- df14 %>% filter(color2 == "disconcordant") 
disconcordant_num <- nrow(disconcordant)
disconcordant_num
rm(concordant, concordant_num, disconcordant, disconcordant_num)

p3 <- df14 %>% ggplot(aes(x = ASB_plot, y = dMOCCS2score, color = color2)) +
  geom_point(size = 0.5, alpha = 0.7) +
  xlab("ASB significance") +
  ylab("dMOCCS2score") +
  scale_colour_manual(
    values = c(
      concordant = "red3",
      disconcordant = "blue3", 
      nondifferent = "gray"
    )
  )+
  ggtitle(paste0(target_tf, " ","dMOCCS2score q<0.01, threshold q<0.05"  ))+
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

#p3
#ggsave(paste0("/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/plot/", target_tf, "_dMOCCS001_ASB005.pdf"), p2, width = 7, height = 7)
rm(df14)

# 9.filter ASB significance and number of concordant and disconcordant snps (dMOCCS2score q<0.01, threshold q<0.01)
threshold <- -log10(0.01)
df11 <- df6 %>% mutate(color1 = ifelse((ASB_plot < -threshold  & q_value < 0.01 & dMOCCS2score > 0) | (ASB_plot > threshold  & q_value < 0.01 & dMOCCS2score < 0), "concordant", "disconcordant")) 
df12 <- df11 %>% filter(color1 == "concordant") %>% mutate(color2 = color1)
df13 <- df11 %>%  filter(color1 == "disconcordant") %>% mutate(color2 = ifelse((q_value < 0.01 & abs(ASB_plot) > threshold), "disconcordant", "nondifferent")) 
rm(df11)
df14 <- rbind(df12, df13)
rm(df12, df13)

concordant <- df14 %>% filter(color2 == "concordant")
concordant_num <- nrow(concordant)
concordant_num

disconcordant <- df14 %>% filter(color2 == "disconcordant") 
disconcordant_num <- nrow(disconcordant)
disconcordant_num
rm(concordant, concordant_num, disconcordant, disconcordant_num)

p4 <- df14 %>% ggplot(aes(x = ASB_plot, y = dMOCCS2score, color = color2)) +
  geom_point(size = 0.5, alpha = 0.7) +
  xlab("ASB significance") +
  ylab("dMOCCS2score") +
  scale_colour_manual(
    values = c(
      concordant = "red3",
      disconcordant = "blue3", 
      nondifferent = "gray"
    )
  )+
  ggtitle(paste0(target_tf, " ", "dMOCCS2score q<0.01, threshold q<0.01"))+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_line(colour="gray"),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=10,face="bold"),
        axis.title=element_text(size=14)
  )
#p4
#ggsave(paste0("/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/plot/", target_tf, "_dMOCCS001_ASB001.pdf"), p4, width = 7, height = 7)
rm(df14)


# 10. negative control summary
df_all <- readRDS(paste0("/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/data/dMOCCS2score_output_all/ADASTRA_dMOCCS2score_", target_tf, "_all_part1.rds")) #CTCF
df_all2 <- df_all %>% mutate(position = end - 6 + posi) 
rm(df_all)
target_df <- df_all2  %>% drop_na(dMOCCS2score)
rm(df_all2)
dMOCCS2score_rand <- sample(target_df$dMOCCS2score, nrow(target_df))
target_df_rand <- target_df %>% select(-dMOCCS2score) %>% mutate(dMOCCS2score = dMOCCS2score_rand)
rm(target_df, dMOCCS2score_rand)

df3 <- target_df_rand %>% mutate(fdrp_bh_ref_log = -log10(fdrp_bh_ref),fdrp_bh_alt_log = -log10(fdrp_bh_alt), ASB = ifelse(fdrp_bh_ref_log <= fdrp_bh_alt_log, fdrp_bh_alt_log, fdrp_bh_ref_log), x_axis = ifelse(fdrp_bh_ref_log <= fdrp_bh_alt_log, "positive", "negative"))
df4 <- df3 %>% filter(x_axis == "positive") %>% mutate(ASB_plot = ASB)
df5 <- df3 %>% filter(x_axis == "negative") %>% mutate(ASB_plot = -ASB)
rm(df3, target_df_rand)
df6 <- rbind(df4, df5)
rm(df4, df5)

threshold <- -log10(0.01)
df11 <- df6 %>% mutate(color1 = ifelse((ASB_plot < -threshold  & q_value < 0.01 & dMOCCS2score > 0) | (ASB_plot > threshold  & q_value < 0.01 & dMOCCS2score < 0), "concordant", "disconcordant")) 
rm(df6)
df12 <- df11 %>% filter(color1 == "concordant") %>% mutate(color2 = color1)
df13 <- df11 %>%  filter(color1 == "disconcordant") %>% mutate(color2 = ifelse((q_value < 0.01 & abs(ASB_plot) > threshold), "disconcordant", "nondifferent")) 
rm(df11)
df14 <- rbind(df12, df13)
rm(df12, df13)

concordant <- df14 %>% filter(color2 == "concordant")
concordant_num <- nrow(concordant)
concordant_num

disconcordant <- df14 %>% filter(color2 == "disconcordant") 
disconcordant_num <- nrow(disconcordant)
disconcordant_num
rm(concordant, concordant_num, disconcordant, disconcordant_num)

p5 <- df14 %>% ggplot(aes(x = ASB_plot, y = dMOCCS2score, color = color2)) +
  geom_point(size = 0.5, alpha = 0.5) +
  xlab("ASB significance") +
  ylab("dMOCCS2score") +
  scale_colour_manual(
    values = c(
      concordant = "red3",
      disconcordant = "blue3", 
      nondifferent = "gray"
    )
  )+
  ggtitle(paste0(target_tf, " ", "all samples", " ", "dMOCCS2score q<0.01, threshold q<0.05"  ))+
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
#p5
#ggsave(paste0("/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/plot/", target_tf, "_negacon_dMOCCS001_ASB001.pdf"), p5, width = 7, height = 7)


p6 <- df14 %>% filter(color2 != "nondifferent") %>%
  ggplot(aes(x = ASB_plot, y = dMOCCS2score, color = color2)) +
  geom_point(size = 0.5, alpha = 0.7) +
  xlab("ASB significance") +
  ylab("dMOCCS2score") +
  scale_colour_manual(
    values = c(
      concordant = "red3",
      disconcordant = "blue3", 
      nondifferent = "gray"
    )
  )+
  ggtitle(paste0(target_tf, " ", "dMOCCS2score q<0.01, threshold q<0.05"  ))+
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
#p6


# 11. p value from negative control's experimental distribution
tf_list <- c("BRD4", "ESR1", "SPI1", "EP300", "JUND", "MYCN", "MAX", "EZH2", "SUMO2", "CEBPB", "GATA2", "GATA3", "HDAC2", "FOS")
threshold <- -log10(0.01)
ecdf_pval_list <- list()
for (i in 1:length(tf_list)) {
  
  target_tf <- tf_list[i]
  df <- readRDS(paste0("/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/data/dMOCCS2score_output_all/ADASTRA_dMOCCS2score_", target_tf, "_all.rds"))
  df3 <- df %>% mutate(fdrp_bh_ref_log = -log10(fdrp_bh_ref),fdrp_bh_alt_log = -log10(fdrp_bh_alt), ASB = ifelse(fdrp_bh_ref_log <= fdrp_bh_alt_log, fdrp_bh_alt_log, fdrp_bh_ref_log), x_axis = ifelse(fdrp_bh_ref_log <= fdrp_bh_alt_log, "positive", "negative"))
  df4 <- df3 %>% filter(x_axis == "positive") %>% mutate(ASB_plot = ASB)
  df5 <- df3 %>% filter(x_axis == "negative") %>% mutate(ASB_plot = -ASB)
  df6 <- rbind(df4, df5)
  rm(df3, df4, df5)
  
  df11 <- df6 %>% mutate(color1 = ifelse((ASB_plot < -threshold  & q_value < 0.01 & dMOCCS2score > 0) | (ASB_plot > threshold  & q_value < 0.01 & dMOCCS2score < 0), "concordant", "disconcordant")) 
  rm(df6)
  df12 <- df11 %>% filter(color1 == "concordant") %>% mutate(color2 = color1)
  df13 <- df11 %>%  filter(color1 == "disconcordant") %>% mutate(color2 = ifelse((q_value < 0.01 & abs(ASB_plot) > threshold), "disconcordant", "nondifferent")) 
  rm(df11)
  df14 <- rbind(df12, df13)
  rm(df12, df13)
  
  concordant <- df14 %>% filter(color2 == "concordant")
  concordant_num <- nrow(concordant)
  disconcordant <- df14 %>% filter(color2 == "disconcordant") 
  disconcordant_num <- nrow(disconcordant)
  concordant_ratio <- concordant_num/(concordant_num+disconcordant_num)
  
  
  # shuffle dMOCCS2score and build empirical distribution
  target_df <- df  %>% drop_na(dMOCCS2score)
  concordant_ratio_rand_list <- c()
  for (j in 1:100) {
    dMOCCS2score_rand <- sample(target_df$dMOCCS2score, nrow(target_df))
    target_df_rand <- target_df %>% select(-dMOCCS2score) %>% mutate(dMOCCS2score = dMOCCS2score_rand)
    df3_rand <- target_df_rand %>% mutate(fdrp_bh_ref_log = -log10(fdrp_bh_ref),fdrp_bh_alt_log = -log10(fdrp_bh_alt), ASB = ifelse(fdrp_bh_ref_log <= fdrp_bh_alt_log, fdrp_bh_alt_log, fdrp_bh_ref_log), x_axis = ifelse(fdrp_bh_ref_log <= fdrp_bh_alt_log, "positive", "negative"))
    df4_rand <- df3_rand %>% filter(x_axis == "positive") %>% mutate(ASB_plot = ASB)
    df5_rand <- df3_rand %>% filter(x_axis == "negative") %>% mutate(ASB_plot = -ASB)
    df6_rand <- rbind(df4_rand, df5_rand)
    rm(df3_rand, df4_rand, df5_rand)
    
    df11_rand <- df6_rand %>% mutate(color1 = ifelse((ASB_plot < -threshold  & q_value < 0.01 & dMOCCS2score > 0) | (ASB_plot > threshold  & q_value < 0.01 & dMOCCS2score < 0), "concordant", "disconcordant")) 
    rm(df6_rand)
    df12_rand <- df11_rand %>% filter(color1 == "concordant") %>% mutate(color2 = color1)
    df13_rand <- df11_rand %>%  filter(color1 == "disconcordant") %>% mutate(color2 = ifelse((q_value < 0.01 & abs(ASB_plot) > threshold), "disconcordant", "nondifferent")) 
    rm(df11_rand)
    df14_rand <- rbind(df12_rand, df13_rand)
    rm(df12_rand, df13_rand)
    
    concordant_rand <- df14_rand %>% filter(color2 == "concordant")
    concordant_num_rand <- nrow(concordant_rand)
    disconcordant_rand <- df14_rand %>% filter(color2 == "disconcordant") 
    disconcordant_num_rand <- nrow(disconcordant_rand)
    concordant_ratio_rand <- concordant_num_rand/(concordant_num_rand + disconcordant_num_rand)
    if(j == 1){
      concordant_ratio_rand_list <- concordant_ratio_rand
    }else{
      concordant_ratio_rand_list <- c(concordant_ratio_rand_list, concordant_ratio_rand)
    }
    
  }
  
  negative_samples <- concordant_ratio_rand_list
  func_pval <- ecdf(negative_samples)
  ecdf_pval_list[[target_tf]] <- 1 - func_pval(concordant_ratio)
  
  
  
  
}
saveRDS(ecdf_pval_list, "/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/data/ecdf_pval_list.rds")

# check p val from empircal distribution
for (t in 1:length(tf_list)) {
  target_tf <- tf_list[t]
  print(target_tf)
  print(ecdf_pval_list[[target_tf]])
}