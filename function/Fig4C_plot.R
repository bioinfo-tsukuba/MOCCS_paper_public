Fig4C_plot <- function(target_ID1 = "", 
                       target_ID2 = "", 
                       target_TF = ""){
  
  totalization <- readRDS("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/MOCCSout_hg38_hard_filter_annotated.rds")
  
  target_MOCCS1 <- totalization %>% filter(ID == target_ID1)
  colnames(target_MOCCS1) <- c()
  target_MOCCS2 <- totalization %>% filter(ID == target_ID2)
  colnames(target_MOCCS2) <- c()
  
  df_join <- target_MOCCS1 %>% left_join(target_MOCCS2, by = "kmer")
  
  # calculate dMOCCS2score and p value
  
  
  # plot
  p1 <- df_join %>% ggplot(aes(x = MOCCS2score1, y = MOCCS2score2, color = q_value)) +
    geom_point() +
    ggtitle(target_TF) +
    xlab(target_ID1) +
    ylab(target_ID2) +
    scale_fill_manual(values = c("#ff0000", "#000080")) +
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
  
  
  
  return(p1)
}