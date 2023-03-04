library(tidyverse)
target_phenotype <- "CD"

if(target_phenotype == "CD"){
  #target_df <- readRDS(url("https://figshare.com/ndownloader/files/34065827", "rb"))
  target_df <- readRDS(url("https://figshare.com/ndownloader/files/39276440", "rb"))
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

head(target_df)
dim(target_df)

df_phenotype_binded_all_selected_annotated <- target_df %>% left_join(annotation, by = "ID")

tmp <- df_phenotype_binded_all_selected_annotated %>% filter(Antigen == target_tf)  %>% as_tibble()
head(tmp)


tmp2 <- tmp %>% select(-Trait, -snp) %>% filter(alt_bp != "?") %>% drop_na(dMOCCS2score)%>% distinct()
tmp2_within <- tmp2 %>% filter(peak == "within peaks")
tmp2_without <- tmp2 %>% filter(peak == "without peaks")


tmp3_within <- tmp2_within %>% mutate(significant = ifelse(q_val < 0.05, "significant", "non-significant")) %>% 
  unite("chr_posi", c(chr, chr_position), sep = "_") %>% group_by(ID, chr_posi) %>% summarise(max_abs_dMO = max(abs(dMOCCS2score)))
join_key_within <- tmp3_within %>% unite("key", c(ID, chr_posi), sep = "_") %>% .$key
tmp3_2_within <- tmp2_within %>% unite("chr_posi", c(chr, chr_position), sep = "_") %>% unite("key", c(ID, chr_posi), sep = "_") %>% select(key, q_val) %>%
  group_by(key) %>% summarise(min_qval = min(q_val))
tmp3_3_within <- tmp3_within %>% ungroup() %>% mutate(key = join_key_within) %>% left_join(tmp3_2_within, by = "key") %>% mutate(peak = "within peaks")

tmp3_without <- tmp2_without %>% mutate(significant = ifelse(q_val < 0.05, "significant", "non-significant")) %>% 
  unite("chr_posi", c(chr, chr_position), sep = "_") %>% group_by(ID, chr_posi) %>% summarise(max_abs_dMO = max(abs(dMOCCS2score)))
join_key_without <- tmp3_without %>% unite("key", c(ID, chr_posi), sep = "_") %>% .$key
tmp3_2_without <- tmp2_without %>% unite("chr_posi", c(chr, chr_position), sep = "_") %>% unite("key", c(ID, chr_posi), sep = "_") %>% select(key, q_val) %>%
  group_by(key) %>% summarise(min_qval = min(q_val))
tmp3_3_without <- tmp3_without %>% ungroup() %>% mutate(key = join_key_without) %>% left_join(tmp3_2, by = "key")%>% mutate(peak = "without peaks")

tmp3_3 <- rbind(tmp3_3_within, tmp3_3_without)

tmp3_4 <- tmp3_3 %>% mutate(significant = ifelse(min_qval < 0.05, "significant", "non-significant")) %>% 
  group_by(ID, peak) %>% summarise(snp_num = n()) 

# すべてのSNPの中の、within/withoutの割合
tmp3_4 %>%ggplot(aes(x = reorder(ID, -snp_num),  y = snp_num, fill = peak)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#DC143C", "#696969")) +
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
  ylab("Number of GWAS-SNPs") +
  labs(fill = "") 

# within peaks
tmp3_3$max_abs_dMO <- as.numeric(tmp3_3$max_abs_dMO)
tmp3_3$min_qval <- as.numeric(tmp3_3$min_qval)
tmp3_3 %>% filter(peak == "within peaks") %>% 
  mutate(significant = ifelse(min_qval < 0.05, "significant", "non-significant")) %>% 
  ggplot(aes(x = reorder(key, -max_abs_dMO), y = max_abs_dMO, color= significant)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#696969", "#DC143C")) +
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_line(colour="gray"),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=4,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        #legend.position = 'none',
        legend.title = element_blank(),
        aspect.ratio = 1
  )

# あるCHIP-seqサンプルに絞る
# FOXA1: SRX4815198
# GATA3: SRX1041802
# FOS: SRX150478
tmp3_3 %>% filter(peak == "within peaks" & ID == "SRX4815198") %>% 
  mutate(significant = ifelse(min_qval < 0.05, "significant", "non-significant")) %>% 
  ggplot(aes(x = reorder(key, -max_abs_dMO), y = max_abs_dMO, color= significant)) +
  geom_point(size = 5) +
  ggtitle("SRX4815198")+
  scale_color_manual(values = c("#696969", "#DC143C")) +
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_line(colour="gray"),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=4,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        #legend.position = 'none',
        legend.title = element_blank(),
        aspect.ratio = 1
  )

# without peaks

tmp3_3 %>% filter(peak == "without peaks") %>% 
  mutate(significant = ifelse(min_qval < 0.05, "significant", "non-significant")) %>% 
  ggplot(aes(x = reorder(key, -max_abs_dMO), y = max_abs_dMO, color= significant)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("#696969", "#DC143C")) +
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=0,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        #legend.position = 'none',
        legend.title = element_blank(),
        aspect.ratio = 1
  )

#tmp2 <- tmp %>% mutate(significant = ifelse(q_val < 0.05, "significant", "non-significant"))
#tmp3 <- tmp2 %>% unite("chr_posi", c(chr, chr_position), sep = "_") %>% group_by(ID, chr_posi) %>% summarise(max_abs_dMO = max(abs(dMOCCS2score)))
#join_key <- tmp3 %>% unite("key", c(ID, chr_posi), sep = "_") %>% .$key
#tmp3_2 <- tmp2 %>% unite("chr_posi", c(chr, chr_position), sep = "_") %>% unite("key", c(ID, chr_posi), sep = "_") %>% select(key, q_val) %>%
  #group_by(key) %>% summarise(min_qval = min(q_val))
#tmp3_3 <- tmp3 %>% ungroup() %>% mutate(key = join_key) %>% left_join(df3_2, by = "key")
#tmp3_4 <- tmp3_3 %>% mutate(significant = ifelse(min_qval < 0.05, "significant", "non-significant")) %>% 
  #group_by(ID, significant) %>% summarise(snp_num = n()) 
