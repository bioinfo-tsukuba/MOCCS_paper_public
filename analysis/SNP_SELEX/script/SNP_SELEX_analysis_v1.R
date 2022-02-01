rm(list=ls())
library(tidyverse)

# 1. summarise data 
target_tf <- "NR3C1"
print(target_tf)
target_df <- read_tsv(paste0("/home/s-tahara/allele_binding_SNP_hg38/SNP_SELEX/data/dMOCCS2score/dMOCCS2score_hg38_", target_tf, "_all.tsv" ))
totalization <- readRDS("/home/s-tahara/MOCCS-DB/WORK/MOCCS/MOCCS-ANALYSIS/output_hg38_v2/MOCCSout_hg38_all_qval_annotated.rds")
ID_hard <- readRDS("/home/s-tahara/DROMPA/code_hg38_2/hg38_hard_filter_ID.rds")
annotation <- totalization %>% select(ID, Antigen, Cell_type_class, Cell_type) %>% filter(Antigen == target_tf) %>% distinct()
target_df <- target_df %>% left_join(annotation, by = "ID")

## statistics
sample_num <- totalization %>% filter(ID %in% ID_hard) %>% filter(Antigen == target_tf) %>% .$ID %>% unique() %>% as.character()
length(sample_num)

snp_list <- target_df %>% mutate(snp_posi = end -6 + posi) %>% unite("snp", c(chr, snp_posi)) %>% .$snp %>% unique() %>% as.character()
length(snp_list)

pbs <- target_df %>% .$pbs %>% as.numeric()
summary(pbs)

dMOCCS2score <- target_df %>% .$dMOCCS2score %>% as.numeric()
summary(dMOCCS2score)

## cor
out <- c()
p_list <- c()
for(x in 1:6){
  target_df_posi <- target_df %>% filter(posi == x)
  a <- target_df_posi$dMOCCS2score %>% as.numeric()
  b <- target_df_posi$pbs %>% as.numeric()
  cor_out <- cor(a, b, use = "complete.obs", method = "spearman")
  out <- c(out, cor_out)
  out_p <- cor.test(a, b, method = "spearman")$p.value
  p_list <- c(p_list, out_p)
}
out
p_list


out2 <- c()
p_list2 <- c()
for(x in 1:6){
  target_df_posi <- target_df %>% filter(posi == x)
  a <- target_df_posi$dMOCCS2score %>% as.numeric()
  b <- target_df_posi$pbs %>% as.numeric()
  cor_out <- cor(a, b, use = "complete.obs", method = "pearson")
  out2 <- c(out2, cor_out)
  out_p2 <- cor.test(a, b, method = "pearson")$p.value
  p_list2 <- c(p_list2, out_p2)
}
out2
p_list2


x <- target_df$dMOCCS2score %>% as.numeric()
y <- target_df$pbs %>% as.numeric()
cor(x, y, use = "complete.obs", method = "spearman")
cor.test(x, y,  method = "spearman")$p.value
cor(x, y, use = "complete.obs", method = "pearson")
cor.test(x, y,  method = "pearson")$p.value

## plot
p1 <- target_df %>% ggplot(aes(x = dMOCCS2score, y = pbs))+
  geom_point(size = 0.5, alpha = 0.3)+
  geom_smooth(se = FALSE, method = lm) +
  ggtitle(target_tf)+
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
p1
ggsave(paste0("/home/s-tahara/allele_binding_SNP_hg38/SNP_SELEX/plot/", target_tf, "_dMOCCS2score_pbs_plot.pdf"))
