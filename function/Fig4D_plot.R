Fig4D_plot <- function(target_ID1 = "", 
                       target_ID2 = "", 
                       target_TF = ""){
  
  totalization <- readRDS("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/MOCCSout_hg38_hard_filter_annotated.rds")
  
  target_MOCCS1 <- totalization %>% filter(ID == target_ID1)
  colnames(target_MOCCS1) <- c()
  target_MOCCS2 <- totalization %>% filter(ID == target_ID2)
  colnames(target_MOCCS2) <- c()
  
  df_join <- target_MOCCS1 %>% left_join(target_MOCCS2, by = "kmer")
  
  p1 <- df_join %>% ggplot(aes(x = MOCCS2score1, y = MOCCS2score2)) +
    geom_point() +
    ggtitle(target_CL) +
    xlab(target_ID1) +
    ylab(target_ID2)
  
  
  
  return(p1)
}