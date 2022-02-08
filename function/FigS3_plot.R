FigS3_plot <- function(annotation_path, path){
  
  library(tidyverse)
  totalization_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds"
  totalization <- readRDS(totalization_path)
  ID_hard <- readRDS(paste0(annotation_path, "hg38_hard_filter_ID.rds"))
  
  experimentList_tab4 <- readRDS("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/experimentList_tab4.rds")
  experimentList_tab5 <- experimentList_tab4 %>% filter(Genome == "hg38" & Antigen_class == "TFs and others" & ID %in% ID_hard) %>% select(ID, peaks)
  totalization_joined <- totalization %>% filter(ID %in% ID_hard) %>% left_join(experimentList_tab5, by = "ID") %>% distinct()
  
  # 1. peak and sample count
  peak_count_df <- totalization_joined %>% select(ID, peaks) %>% distinct()
  length(unique(peak_count_df$ID))
  p1 <- peak_count_df %>% ggplot(aes(x = log(peaks))) +
    geom_histogram(bins = 100)+
    ylab("sample count") +
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
  
  
  
  # 2. peak and number of significant k-mer
  p2_df <- totalization_joined %>% mutate(signigicance = ifelse(q_value < 0.05, "significant", "non significant")) %>%
    filter(signigicance == "significant") %>%
    group_by(ID, peaks) %>%
    summarise(sig_num = n())
  
  p2 <- p2_df %>% ggplot(aes(x = log(peaks), y = sig_num)) +
    geom_bin2d(bins = 50) +
    scale_fill_continuous(type = "viridis") +
    ylab("significant k-mer count") +
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
  
  p3_df <- totalization_joined %>% mutate(signigicance = ifelse(q_value < 0.01, "significant", "non significant")) %>%
    filter(signigicance == "significant") %>%
    group_by(ID, peaks) %>%
    summarise(sig_num = n())
  
  p3 <- p3_df %>% ggplot(aes(x = log(peaks), y = sig_num)) +
    geom_bin2d(bins = 50) +
    scale_fill_continuous(type = "viridis") +
    ylab("significant k-mer count (q<0.01)") +
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