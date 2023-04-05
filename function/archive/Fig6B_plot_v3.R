Fig6B_plot <- function( annotation_path){
  
  library(tidyverse)
  Fig6B_df <- read_csv(paste0(annotation_path, "Fig6B_df_20230110.csv"))
  Fig6B_df <- Fig6B_df %>% unite("color", c(peak, significant), sep = "_")
  
  # Fisher
  pheno_list <- unique(Fig6B_df$phenotype)
  p_list <- list()
  for (tgt_pheno in pheno_list) {
    tgt_df <- Fig6B_df %>% filter(phenotype == tgt_pheno)
    a <- tgt_df %>% filter(color == "within_peak_sig") %>% .$snp_num %>% as.numeric()
    b <- tgt_df %>% filter(color == "within_peak_non-sig") %>% .$snp_num %>% as.numeric()
    c <- tgt_df %>% filter(color == "out_peak_sig") %>% .$snp_num %>% as.numeric()
    d <- tgt_df %>% filter(color == "out_peak_non-sig") %>% .$snp_num %>% as.numeric()
    mx <- matrix(c(a, b, c, d), nrow=2, byrow=T)
    colnames(mx) <- c("dMOCCS2score_sig", "dMOCCS2score_nonsig")
    rownames(mx) <- c("within peaks", "out of peaks")
    print(list(tgt_pheno, mx))
    
    #フィッシャーの正確確率検定
    fisher_res <- fisher.test(mx)
    p_list[tgt_pheno] <- fisher_res$p.value
  }
  
  p <- c(rep(p_list[1], 4), rep(p_list[2], 4), rep(p_list[3], 4), rep(p_list[4], 4)) %>% as.numeric()
  Fig6B_df2 <- Fig6B_df %>% mutate(fisher_p = p)
  Fig6B_df2$color <- factor(Fig6B_df2$color, levels = c("out_peak_non-sig", "out_peak_sig", "within_peak_non-sig", "within_peak_sig"))
  Fig6B_plot <- Fig6B_df2 %>% ggplot(aes(x = phenotype, y = snp_num, fill = color)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    scale_fill_manual(values = c("#696969", "green4", "#4169E1", "#DC143C")) +
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_blank(),
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
    ggtitle("Number of SNPs") +
    labs(fill = "") +
    coord_flip()
  print(Fig6B_df2)
  
  return(Fig6B_plot)
  
  
}