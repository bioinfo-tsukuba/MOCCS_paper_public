Fig5C_plot <- function(target_TF, annotation_path){
  
  library(tidyverse)
  target_df <- read_tsv(paste0(annotation_path, "SNP_SELEX/dMOCCS2score_hg38_", target_TF, "_all.tsv"))
  totalization <- readRDS("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds")
  annotation <- totalization %>% select(ID, Antigen, Cell_type_class, Cell_type) %>% filter(Antigen == target_TF) %>% distinct()
  target_df <- target_df %>% left_join(annotation, by = "ID")
  
  p1 <- target_df %>% ggplot(aes(x = dMOCCS2score, y = pbs))+
    geom_point(size = 0.5, alpha = 0.3)+
    geom_smooth(se = FALSE, method = lm) +
    ggtitle(target_TF)+
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