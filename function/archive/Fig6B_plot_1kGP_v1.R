Fig6B_plot_1kGP_v1 <- function(){
  
  # based on the results; @NIG, /home/s-tahara/allele_binding_SNP_hg38/GWAS/script_20230108/1000genome_make_summary_tib_Fig5B.R
  
  library(tidyverse)
  library(patchwork)
  target_phenotype_list <- c("SLE", "MS", "IBD", "CD")
  p_patch <- c()
  
  for (target_phenotype in target_phenotype_list) {
    print(target_phenotype)
    summary_df_1kGP <- read_tsv(paste0("~/MOCCS_paper_public/data/Fig6/1kGP/", target_phenotype, "_summary_TF42.tsv"))
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
    ggsave(paste0("~/MOCCS_paper_public/plot/Fig6/snp_sig_ratio_1kGP_", target_phenotype, "_TF42.pdf"), p)
    
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
    ggsave(paste0("~/MOCCS_paper_public/plot/Fig6/snp_sig_ratio_1kGP_", target_phenotype, "_TF42_allsample.pdf"), p2)
    
    if(target_phenotype == target_phenotype_list[1]){
      p_patch <- p
    }else{
      p_patch <- p_patch + p
    }
  }#for (target_phenotype in target_phenotype_list)
  return(p_patch)
}