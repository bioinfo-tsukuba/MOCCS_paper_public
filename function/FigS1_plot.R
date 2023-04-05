FigS1_plot <- function(path){
  
  library(tidyverse)
  library(GGally)
  options(timeout=1000)
  df_tidy <- readRDS(paste0(path, "DROMPA_SUMMARY_hg38.rds"))
  experimentList_tab4 <- readRDS(url("https://figshare.com/ndownloader/files/34065671", "rb"))
  experimentList_tab5 <- experimentList_tab4 %>% filter(Genome == "hg38" & Antigen_class == "TFs and others") %>% select(ID, peaks)
  df_tidy_joined <- df_tidy %>% left_join(experimentList_tab5, by = "ID")  
  df_tidy_joined_selected <- df_tidy_joined %>% filter(Library_complexity > 0.8 & total_reads > 10000000 & GC_content < 60 & NSC > 2.0 & Bu  > 0.8)
  hard_filter_ID <- df_tidy_joined_selected$ID %>% unique() %>% as.character()
  
  
  g_pair <- df_tidy_joined %>% 
    mutate(filter = ifelse(ID %in% hard_filter_ID, "hard", "others")) %>% 
    select(-ID) %>% filter(Bu != Inf) %>% 
    ggpairs(lower=list(continuous=wrap("points",size=0.05, alpha = 0.5)), 
            diag = list(continuous = wrap("densityDiag", alpha = 0.5)), 
            mapping = aes(color = filter))+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_line(colour="gray"),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
          axis.text.y =element_text(size=10,face="bold"),
          axis.title=element_text(size=14,face="bold")
    ) +
    theme_classic()

  
  return(g_pair)
}