Fig1G_plot <- function(target_TF, annotation_path){
  
  library(tidyverse)
  library(RColorBrewer)
  library(Biostrings)
  
  totalization_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds"
  totalization <- readRDS(totalization_path)
  ID_hard <- readRDS(paste0(annotation_path, "hg38_hard_filter_ID.rds"))
  
  # filter target TF PWM table
  PWM_table_all <- readRDS("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/PWM_likelihood_HOMER.rds")
  target_PWM <- PWM_table_all[[target_TF]]
  target_PWM2 <- target_PWM %>% group_by(kmer) %>% summarise(max_PWMscore = max(PWMscore))
  
  # join MOCCS output and PWM table
  target_MOCCS <- totalization %>% filter(ID %in% ID_hard) %>% filter(Antigen == target_TF)
  
    # PWMのk-merをMOCCSの片方のk-merに絞る (相補鎖を考慮する)
    target_kmer <- unique(target_MOCCS$kmer)
    non_target_kmer <- setdiff(unique(target_PWM2$kmer), target_kmer)
    target_PWM3 <- target_PWM2 %>% mutate(kmer_class = ifelse(kmer %in% target_kmer, "MOCCS_kmer", "other_kmer"))
    for (r in 1:nrow(target_PWM3)) {
      target_row <- target_PWM3[r,]
      if(target_row$kmer_class == "other_kmer"){
        target_row$kmer <- reverseComplement(DNAStringSet(target_row$kmer)) %>% as.character()
      }
      if(r == 1){
        target_PWM4 <- target_row
      }else{
        target_PWM4 <- target_PWM4 %>% add_row(target_row)
      }
    }
  target_PWM5 <- target_PWM4 %>% group_by(kmer) %>% summarise(max_PWMscore = max(max_PWMscore))
  
  df_join1 <- target_MOCCS %>% left_join(target_PWM5, by = "kmer")
  df_join1[is.na(df_join1)] <- 0 #NAを0に置換
  
  
  
    
  # PWM scoreをsignificant k-merとそれ以外で分けてplot
  df_join2 <- df_join1 %>% mutate(sig_kmer = ifelse(q_value < 0.05, "significant k-mer", "non significant k-mer")) 
  df_join3 <- transform(df_join2, sig_kmer = factor(sig_kmer, levels = c("non significant k-mer", "significant k-mer")))
  #p <- df_join3 %>% ggplot(aes(x=sig_kmer, y = max_PWMscore, fill = sig_kmer)) +
    #geom_violin()  +
    #geom_beeswarm(size = 0.2) +
    #xlab("k-mer")+
    #theme(plot.title = element_text(face="bold",hjust = 0.5), 
          #panel.grid.major = element_line(colour = "gray"),
          #panel.grid.minor = element_line(colour="gray"),
          #panel.background = element_blank(), 
          #axis.line = element_line(colour="black"),
          #axis.text=element_text(size=12,face="bold"),
          #axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
          #axis.text.y =element_text(size=10,face="bold"),
          #axis.title=element_text(size=14,face="bold")
    #)
  
  x <- df_join3 %>% filter(sig_kmer == "significant k-mer") %>% .$max_PWMscore %>% as.numeric()
  y <- df_join3 %>% filter(sig_kmer == "non significant k-mer") %>% .$max_PWMscore %>% as.numeric()
  
  ks <- suppressMessages(ks.test(x,y,alternative="two.sided"))
  p_value <- ks$p.value
  
  
  p2 <- df_join3 %>% ggplot(aes(x=-log2(max_PWMscore), color = sig_kmer)) +
    stat_ecdf(size = 0.7) +
    xlab("-log2(PWMscore)")+
    ylab("probability") +
    scale_color_manual(values = c("#696969", "#DC143C")) +
    guides(color = guide_legend(reverse = TRUE)) +
    annotate("text",x=-Inf,y=Inf,label=paste0("p-value = ", p_value),hjust=-.2,vjust=2, size=5, fontface = "bold") +
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_line(colour="gray"),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
          axis.text.y =element_text(size=10,face="bold"),
          axis.title=element_text(size=14,face="bold"),
          legend.position = c(0.85,0.2),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.box.background = element_rect(colour = "black")
    )
  
  return(p2)
  
}