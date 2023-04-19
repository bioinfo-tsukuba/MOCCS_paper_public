Fig3B_plot <- function(simu_name){
  
  library(tidyverse)
  df <- readRDS(paste0("~/MOCCS_paper_public/data/Fig4/", simu_name, "_new2.rds"))
  
  p <- df %>% mutate(color = ifelse(q_value < 0.05, "differential", "nondifferential")) %>%
    ggplot(aes(x = MOCCS2score, y = MOCCS2score2, color = color)) +
    geom_point(size = 1) +
    scale_colour_manual(
      values = c(
        differential = "red",
        nondifferential = "gray"
      )
    )+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_line(colour="gray"),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
          axis.text.y =element_text(size=10,face="bold"),
          axis.title=element_text(size=14,face="bold"),
          legend.title = element_blank(),
          aspect.ratio = 1
    ) +
    xlab("MOCCS2score1") +
    ggtitle(paste0(simu_name, "_q005"))
  
  return(p)
  
}