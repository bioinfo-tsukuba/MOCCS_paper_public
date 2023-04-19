library(tidyverse)
library(patchwork)
target_phenotype_list <- c("SLE", "MS", "IBD", "CD")
p_patch <- c()
totalization <- readRDS("~/MOCCS_paper_public/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds")
annotation <- totalization %>% select(ID, Antigen, Cell_type_class, Cell_type) %>% distinct()

patch_list <- list()
for (target_phenotype in target_phenotype_list){
  print(target_phenotype)
  #summary_df <- read_tsv(paste0("~/MOCCS_paper_public/data/Fig6/sig_snp_ratio_", target_phenotype, ".tsv"))
  summary_df <- read_tsv(paste0("~/MOCCS_paper_public/data/Fig6/sig_snp_ratio_", target_phenotype, "_TF42.tsv"))
  summary_df2 <- summary_df %>% mutate(label2 = ifelse(label == "within", "within peaks", "out of peaks"))
  ID_list <- summary_df2$ID %>% unique()
  annotation2 <- annotation %>% filter(ID %in% ID_list)
  
  for (target_tf in unique(summary_df2$TF)) {
    print(target_tf)
    target_df <- summary_df2 %>% filter(TF == target_tf)
    target_annotation <- annotation2 %>% filter(ID %in% unique(target_df$ID))
    target_df2 <- target_df %>% left_join(target_annotation, by = "ID")
    
    p <- target_df2 %>% filter(significant == "significant")%>%
      ggplot(aes(x = Cell_type_class, y = ratio, fill = label2)) +
      geom_boxplot(notchwidth = 0.3, width = 0.7) +
      ggtitle(paste0(target_phenotype, "   ", target_tf)) +
      scale_fill_manual(values = c( "#696969", "#DC143C")) +
      theme(plot.title = element_text(face="bold",hjust = 0.5, size = 8), 
            panel.grid.major = element_line(colour = "gray"),
            panel.grid.minor = element_line(colour="gray"),
            panel.background = element_blank(), 
            axis.line = element_line(colour="black"),
            axis.text=element_text(size=4,face="bold"),
            axis.text.x =element_text(size=5,face="bold", angle = 45, hjust = 1),
            axis.text.y =element_text(size=3,face="bold"),
            axis.title=element_text(size=5,face="bold"),
            legend.position = 'none',
            legend.title = element_blank(),
            aspect.ratio = 0.3
      )+
      ylab("Ratio of significant SNPs") 
    if(target_tf == unique(summary_df2$TF)[[1]]){
      patch_list[[target_phenotype]] <- p
    }else{
      patch_list[[target_phenotype]] <- patch_list[[target_phenotype]] + p
    }
  } # for target_tf
}

summary_df2$TF %>% unique()
patch_list[[1]][[8]] 
patch_list[[2]][[29]]
patch_list[[3]][[40]]
patch_list[[4]][[16]]


# CD, FOS 検定
target_phenotype <- "CD"
target_tf <- "SMARCA4"

target_df <- summary_df2 %>% filter(TF == target_tf)
target_annotation <- annotation2 %>% filter(ID %in% unique(target_df$ID))
target_df2 <- target_df %>% left_join(target_annotation, by = "ID")

# wilcoxon test per TFs
q_list_all <- list()
p_list <- list()
for (tgt_ctc in unique(target_df2$Cell_type_class)) {
  print(tgt_ctc)
  print(target_tf)
  tgt_df <- target_df2 %>% filter(Cell_type_class == tgt_ctc)
  within_df <- tgt_df %>% filter(label == "within") %>% filter(significant == "significant") 
  out_of_df <- tgt_df %>% filter(label == "without") %>% filter(significant == "significant")
  share_ID <- intersect(unique(within_df$ID), unique(out_of_df$ID))
  if(length(share_ID) != 0){
    within_ratio <- within_df %>% filter(ID %in% share_ID) %>% .$ratio %>% as.numeric()
    out_of_ratio <- out_of_df %>% filter(ID %in% share_ID) %>% .$ratio %>% as.numeric()
    
    res <- wilcox.exact(within_ratio, out_of_ratio, paired = T, alternative = "greater")
    p <- res$p.value %>% as.numeric()
    print(p)
    p_list[[tgt_ctc]] <- p
  }
  q_list_all[[target_tf]] <- p.adjust(p_list)
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
