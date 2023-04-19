FigS10_plot_v5 <- function(){
  
  # based on the results; ~/MOCCS_paper_public/function/Sub_function/make_summary_tib_Fig5B.R 
  
  library(tidyverse)
  target_phenotype_list <- c("SLE", "MS", "IBD", "CD")
  
  patch_list <- list()
  for (target_phenotype in target_phenotype_list) {
    print(target_phenotype)
    
    summary_df <- read_tsv(paste0("~/MOCCS_paper_public/data/Fig6/1kGP/", target_phenotype, "_summary_TF42.tsv"))
    summary_df2 <- summary_df %>% mutate(label2 = ifelse(label == "within", "within peaks", "out of peaks"))
    
    patch <- c()
    for (tgt_tf in as.character(unique(summary_df2$TF))) {
      print(tgt_tf)
      p <- summary_df2 %>% drop_na(snp_num) %>% filter(label == "within" & TF == tgt_tf) %>% 
        ggplot(aes(x = reorder(ID, -snp_num),  y = snp_num, fill = significant)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c("#696969", "#DC143C"), labels = c(ratio_nonsig = "non significant", ratio_sig ="significant")) +
        ggtitle(tgt_tf) +
        xlab("ChIP-seq sample IDs") +
        ylab("Number of peak-overlapping snps") +
        theme(plot.title = element_text(face="bold",hjust = 0.5), 
              panel.grid.major = element_line(colour = "gray"),
              panel.grid.minor = element_line(colour="gray"),
              panel.background = element_blank(), 
              axis.line = element_line(colour="black"),
              axis.text=element_text(size=5,face="bold"),
              axis.text.x =element_text(size=0,face="bold", angle = 45, hjust = 1),
              axis.title.x = element_text(size=3,face="bold"),
              axis.text.y =element_text(size=1.5,face="bold"),
              axis.title=element_text(size=5,face="bold"),
              legend.position = 'none',
              legend.title = element_blank(),
              aspect.ratio = 1
        )+
        labs(fill = "") 
      #plot(p)
      
      if(length(patch) == 0){
        patch <- p
      }else{
        patch <- patch + p
      }
    }
    plot(patch)
    #ggsave(paste0("~/MOCCS_paper_public/plot/FigS10/", target_phenotype, "_snp_count_sig_1kGP.pdf"), patch)
    patch_list[[target_phenotype]] <- patch
    
  } #for (target_phenotype in target_phenotype_list) 
  return(patch_list)
}