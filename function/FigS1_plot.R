FigS1_plot <- function(simulation_path){
  
  library(tidyverse)
  library(patchwork)
  df <- readRDS(paste0(simulation_path,"summary_table_for_plot.rds"))
  df_significant <- df %>% filter(sig_or_dif == "significant")
  df_significant2 <- df_significant %>% pivot_longer(cols = c(-a, -N, -sigma, -q, -FDR, -sig_or_dif), names_to = "sensi_or_spe", values_to = "value")
  
  
  # sensitivity
  # alpha
  df_significant2$FDR <- format(df_significant2$FDR, digits = 2)
  p1 <- df_significant2 %>% filter(N == 12000 & sigma == "W5" & q == 0.05) %>%
    ggplot(aes(x = a, y = value, fill = sensi_or_spe))+
    geom_bar(width = 0.5, stat = "identity", position = "dodge") +
    ylim(0,1.2)+
    labs(fill="") +
    scale_fill_manual(values = c("#ff0000", "#000080")) +
    ylab("Sensitivity or Specificity") +
    xlab("α") +
    geom_hline(yintercept=1, linetype="dashed") +
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
          axis.text.y =element_text(size=10,face="bold"),
          axis.title=element_text(size=14,face="bold")
    )
  
  
  # N
  df_significant2$N <- as.character(df_significant2$N)
  dfsignificant2 <- transform(df_significant2, N= factor(N, levels = c("6000", "12000")))
  p2 <- df_significant2   %>% filter(a == "[0.1,0.2]" & sigma == "W5" & q == 0.05) %>%
    ggplot(aes(x = N, y = value, fill = sensi_or_spe))+
    geom_bar(width = 0.5, stat = "identity", position = "dodge") +
    ylim(0,1.2)+
    labs(fill="") +
    scale_fill_manual(values = c("#ff0000", "#000080")) +
    ylab("Sensitivity or Specificity") +
    xlab("N") +
    geom_hline(yintercept=1, linetype="dashed") +
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
          axis.text.y =element_text(size=10,face="bold"),
          axis.title=element_text(size=14,face="bold")
    )
  
  # sigma
  df_significant2 <- transform(df_significant2, sigma= factor(sigma, levels = c("W5", "W2")))
  p3 <- df_significant2 %>% filter(a == "[0.1,0.2]" & N == 12000 & q == 0.05) %>%
    ggplot(aes(x = sigma, y = value, fill = sensi_or_spe)) +
    geom_bar(width = 0.5, stat = "identity", position = "dodge") +
    ylim(0,1.2)+
    labs(fill="") +
    scale_fill_manual(values = c("#ff0000", "#000080")) +
    ylab("Sensitivity or Specificity") +
    xlab("σ") +
    geom_hline(yintercept=1, linetype="dashed") +
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
          axis.text.y =element_text(size=10,face="bold"),
          axis.title=element_text(size=14,face="bold")
    )
  
  # FDR
  # alpha
  df_significant2$FDR<- as.numeric(df_significant2$FDR)
  p4 <- df_significant2 %>% filter(N == 12000 & sigma == "W5" & q == 0.05) %>%
    ggplot(aes(x = a, y = FDR, fill = sensi_or_spe))+
    geom_bar(width = 0.5, stat = "identity", position = "dodge") +
    ylim(c(0, 0.001)) +
    #geom_hline(yintercept=0.05, linetype="dashed") +
    ylab("FDR") +
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
  df_significant2$N <- as.character(df_significant2$N)
  df_significant2 <- transform(df_significant2, N= factor(N, levels = c("6000", "12000")))
  p5 <- df_significant2 %>% filter(a == "[0.1,0.2]" & sigma == "W5" & q == 0.05) %>%
    ggplot(aes(x = N, y = FDR, fill = sensi_or_spe))+
    geom_bar(width = 0.5, stat = "identity", position = "dodge") +
    ylim(c(0, 0.001)) +
    #geom_hline(yintercept=0.05, linetype="dashed") +
    ylab("FDR") +
    xlab("N") +
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
  df_significant2 <- transform(df_significant2, sigma= factor(sigma, levels = c("W5", "W2")))
  p6 <- df_significant2 %>% filter(a == "[0.1,0.2]" & N == 12000 & q == 0.05) %>%
    ggplot(aes(x = sigma, y = FDR, fill = sensi_or_spe)) +
    geom_bar(width = 0.5, stat = "identity", position = "dodge") +
    ylim(c(0, 0.001)) +
    #geom_hline(yintercept=0.05, linetype="dashed") +
    ylab("FDR") +
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
  
  p_all_1 <- (p1 | p2 | p3) 
  p_all_2 <- (p4 | p5 | p6)
  return(list(p_all_1, p_all_2))
  
}