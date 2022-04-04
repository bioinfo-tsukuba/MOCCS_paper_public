Fig5C_plot <- function(target_TF, annotation_path){
  
  library(tidyverse)
  if(target_TF == "HOXB13"){
    target_df <- read_tsv(url("https://figshare.com/ndownloader/files/34336358", "rb"))
  }else if(target_TF == "GATA3"){
    target_df <- read_tsv(url("https://figshare.com/ndownloader/files/34336355", "rb"))
  }else if(target_TF == "FOXA1"){
    target_df <- read_tsv(url("https://figshare.com/ndownloader/files/34336352", "rb"))
  }
  totalization <- readRDS(url("https://figshare.com/ndownloader/files/34065686","rb")) #MOCCSout_hg38_all_qval_annotated.rds
  annotation <- totalization %>% select(ID, Antigen, Cell_type_class, Cell_type) %>% filter(Antigen == target_TF) %>% distinct()
  target_df <- target_df %>% left_join(annotation, by = "ID")
  
  p1 <- target_df  %>% ggplot(aes(x = dMOCCS2score, y = pbs))+
    geom_point(size = 1, alpha = 1)+
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