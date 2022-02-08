FigS2_plot <- function(path){
  
  library(tidyverse)
  library(GGally)
  df_tidy <- readRDS(paste0(path, "DROMPA_SUMMARY_hg38.rds"))
  experimentList_tab4 <- readRDS("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/experimentList_tab4.rds")
  experimentList_tab5 <- experimentList_tab4 %>% filter(Genome == "hg38" & Antigen_class == "TFs and others") %>% select(ID, peaks)
  df_tidy_joined <- df_tidy %>% left_join(experimentList_tab5, by = "ID")  
  df_tidy_joined_selected <- df_tidy_joined %>% filter(Library_complexity > 0.8 & total_reads > 10000000 & GC_content < 60 & NSC > 2.0 & Bu  > 0.8)
  
  g_pair <- df_tidy_joined %>% select(-ID) %>% filter(Bu != Inf) %>% ggpairs(lower=list(continuous=wrap("points",size=0.05)))+
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
  
  g_pair_filtered <- df_tidy_joined_selected %>% select(-ID) %>% filter(Bu != Inf) %>% ggpairs(lower=list(continuous=wrap("points",size=0.05)))+
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
  
  return(list(g_pair, g_pair_filtered))
}