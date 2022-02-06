Fig4B_plot <- function(simulation_path){
  
  library(tidyverse)
  df <- readRDS(paste0(simulation_path,"summary_table_for_plot.rds"))
  df_differential <- df %>% filter(sig_or_dif == "differential")
  df_differential2 <- df_differential %>% pivot_longer(cols = c(-a, -N, -sigma, -q, -FDR, -sig_or_dif), names_to = "sensi_or_spe", values_to = "value")
  
  
  # alpha
  df_differential2$FDR <- format(df_differential2$FDR, digits = 2)
  p1 <- df_differential2 %>% filter(N == 12000 & sigma == "W5" & q == 0.05) %>%
    ggplot(aes(x = a, y = value, fill = sensi_or_spe))+
    geom_bar(width = 0.5, stat = "identity", position = "dodge") +
    geom_text(aes(x=a, y=value+0.05, label=FDR), position=position_dodge2(0.7)) +
    ylim(0,1.2)+
    ylab("Sensitivity or Specificity") +
    xlab("α") +
    labs(fill="") +
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
  
  # N
  df_differential2$N <- as.character(df_differential2$N)
  p2 <- df_differential2 %>% filter(a == "[0.1,0.2]" & sigma == "W5" & q == 0.05) %>%
    ggplot(aes(x = N, y = value, fill = sensi_or_spe))+
    geom_bar(width = 0.5, stat = "identity", position = "dodge") +
    geom_text(aes(x=N, y=value+0.05, label=FDR), position=position_dodge2(0.6)) +
    ylim(0,1.2)+
    ylab("Sensitivity or Specificity") +
    labs(fill="") +
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
  
  # sigma
  p3 <- df_differential2 %>% filter(a == "[0.1,0.2]" & N == 12000 & q == 0.05) %>%
    ggplot(aes(x = sigma, y = value, fill = sensi_or_spe)) +
    geom_bar(width = 0.5, stat = "identity", position = "dodge") +
    geom_text(aes(x=sigma, y=value+0.05, label=FDR), position=position_dodge2(0.6)) +
    ylim(0,1.2)+
    ylab("Sensitivity or Specificity") +
    xlab("σ") +
    labs(fill="") +
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
  
  return(list(p1, p2, p3))
  
}