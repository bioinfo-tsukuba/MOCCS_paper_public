# GWASと1kGPでSNPのかぶりがないか確認 ----
kGP <- readRDS("/Users/saeko/MOCCS_paper_public/data/Fig6/1kGP/CD_1kGP_binded.rds")
GWAS <- readRDS("/Users/saeko/MOCCS_paper_public/data/Fig6/CD_peak_rand_binded_all_qval.rds")
dim(GWAS)

kGP_snp <- kGP %>% unite("key", c(chr, start, end), sep = "_") %>% distinct() %>% .$key %>% unique()
GWAS_snp <- GWAS %>% unite("key", c(chr, start, end), sep = "_") %>% filter(peak == "within peaks") %>% distinct() %>% .$key %>% unique()
intersect(GWAS_snp, kGP_snp)


# dMOCCS2scoreの大きさの違いをplot----
target_phenotype <- "SLE"
kGP <- readRDS(paste0("/Users/saeko/MOCCS_paper_public/data/Fig6/1kGP/",target_phenotype,"_1kGP_binded.rds"))
GWAS <- readRDS(paste0("/Users/saeko/MOCCS_paper_public/data/Fig6/", target_phenotype, "_peak_rand_binded_all_qval.rds"))
colnames(kGP)
colnames(GWAS)
kGP <- kGP %>% mutate(label = "100genome")
GWAS_selected <- GWAS %>% select(chr, start, end, alt_bp, snp_posi, kmer_before_tou, kmer_after_tou, dMOCCS2score, p_value, q_value, ID) %>%
  mutate(lable = "GWAS")
colnames(GWAS_selected) <- c(colnames(kGP), "label")
df_all <- rbind(kGP, GWAS_selected)

df_all %>% ggplot(aes(x = label, y =dMOCCS2score, fill = label)) +
  #geom_violin() +
  geom_boxplot() +
  scale_fill_manual(values = c("#696969", "#DC143C")) +
  theme(plot.title = element_text(face="bold",hjust = 0.5, size = 8), 
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_line(colour="gray"),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=10,face="bold"),
        axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=10,face="bold"),
        axis.title=element_text(size=10,face="bold"),
        #legend.position = 'none',
        legend.title = element_blank(),
        aspect.ratio = 1
  )+
  ggtitle(target_phenotype) +
  ylab("dMOCCS2score") 

df_all %>% ggplot(aes(x = label, y = abs(dMOCCS2score), fill = label)) +
  #geom_violin() +
  geom_boxplot() +
  scale_fill_manual(values = c("#696969", "#DC143C")) +
  theme(plot.title = element_text(face="bold",hjust = 0.5, size = 8), 
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_line(colour="gray"),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=10,face="bold"),
        axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=10,face="bold"),
        axis.title=element_text(size=10,face="bold"),
        #legend.position = 'none',
        legend.title = element_blank(),
        aspect.ratio = 1
  )+
  ggtitle(target_phenotype) +
  ylab("abs_dMOCCS2score") 

#GWASにあるkmerとalt kmerのペアは、1kGPから除去する----
target_phenotype <- "IBD"
kGP <- readRDS(paste0("/Users/saeko/MOCCS_paper_public/data/Fig6/1kGP/", target_phenotype,"_1kGP_binded.rds"))
GWAS <- readRDS(paste0("/Users/saeko/MOCCS_paper_public/data/Fig6/", target_phenotype,"_peak_rand_binded_all_qval.rds"))
GWAS2 <- GWAS %>% filter(peak == "within peaks")

GWAS_ref <- GWAS2$kmer_before_tou %>% unique()
GWAS_alt <- GWAS2$kmer_after_tou %>% unique()

kGP2 <- kGP %>% filter(!ref_kmer %in% GWAS_ref)

saveRDS(kGP2, paste0("/Users/saeko/MOCCS_paper_public/data/Fig6/1kGP/", target_phenotype,"_1kGP_binded_rmGWASrefkmer.rds"))


# Fig5B用にplot ----
# based on the results; @NIG, /home/s-tahara/allele_binding_SNP_hg38/GWAS/script_20230108/1000genome_make_summary_tib_Fig5B.R

library(tidyverse)
library(patchwork)
target_phenotype_list <- c("SLE", "MS", "IBD", "CD")
p_patch <- c()

for (target_phenotype in target_phenotype_list) {
  print(target_phenotype)
  summary_df_1kGP <- read_tsv(paste0("~/MOCCS_paper_public/data/Fig6/1kGP/", target_phenotype, "_summary_TF42_rmGWASrefkmer.tsv"))
  summary_df_1kGP2 <- summary_df_1kGP %>% mutate(label2 = "1kGP")
  shareTF <- summary_df_1kGP2$TF %>% unique()
  
  summary_df_GWAS <- read_tsv(paste0("~/MOCCS_paper_public/data/Fig6/sig_snp_ratio_", target_phenotype, "_TF42.tsv"))
  summary_df_GWAS2 <- summary_df_GWAS %>% mutate(label2 = "GWAS") %>% filter(label == "within" & TF %in% shareTF)
  
  summary_df2 <- rbind(summary_df_GWAS2,  summary_df_1kGP2)
  
  p <- summary_df2 %>% filter(significant == "significant") %>%
    ggplot(aes(x = TF, y = ratio, fill = label2)) +
    geom_boxplot(notchwidth = 0.3, width = 0.7) +
    scale_fill_manual(values = c("#696969", "#DC143C")) +
    theme(plot.title = element_text(face="bold",hjust = 0.5, size = 8), 
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_line(colour="gray"),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=4,face="bold"),
          axis.text.x =element_text(size=5,face="bold", angle = 45, hjust = 1),
          axis.text.y =element_text(size=3,face="bold"),
          axis.title=element_text(size=5,face="bold"),
          #legend.position = 'none',
          legend.title = element_blank(),
          aspect.ratio = 0.1
    )+
    ggtitle(target_phenotype) +
    ylab("Ratio of significant SNPs in peak region") 
  plot(p)
  ggsave(paste0("~/MOCCS_paper_public/plot/Fig6/snp_sig_ratio_1kGP_", target_phenotype, "_TF42_rmGWASrefkmer.pdf"), p)
  
  p2 <- summary_df2 %>% filter(significant == "significant") %>%
    ggplot(aes(x = label2, y = ratio, fill = label2))+
    geom_violin(notchwidth = 0.3, width = 0.7) +
    geom_point() +
    scale_fill_manual(values = c("#696969", "#DC143C")) +
    theme(plot.title = element_text(face="bold",hjust = 0.5, size = 8), 
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_line(colour="gray"),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=4,face="bold"),
          axis.text.x =element_text(size=5,face="bold", angle = 45, hjust = 1),
          axis.text.y =element_text(size=3,face="bold"),
          axis.title=element_text(size=5,face="bold"),
          #legend.position = 'none',
          legend.title = element_blank(),
          aspect.ratio = 1
    )+
    ggtitle(target_phenotype) +
    ylab("Ratio of significant SNPs in peak region") 
  plot(p2)
  ggsave(paste0("~/MOCCS_paper_public/plot/Fig6/snp_sig_ratio_1kGP_", target_phenotype, "_TF42_allsample_rmGWASrefkmer.pdf"), p2)
  
  #if(target_phenotype == target_phenotype_list[1]){
    #p_patch <- p
  #}else{
    #p_patch <- p_patch + p
  #}
}#for (target_phenotype in target_phenotype_list)
