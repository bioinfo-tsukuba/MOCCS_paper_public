Fig1G_plot <- function(annotation_path){
  
  library(tidyverse)
  library(RColorBrewer)

  totalization_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds"
  totalization <- readRDS(totalization_path)
  ID_hard <- readRDS(paste0(annotation_path, "hg38_hard_filter_ID.rds"))
  #ID_soft <- readRDS(paste0(annotation_path, "ID_soft_filter_hg38.rds"))
  Antigen_list <- read_tsv(paste0(annotation_path, "Antigen_list_hg38.txt"), col_names = FALSE)
  Antigen_list <- Antigen_list$X1 %>% as.character() %>% unique()
  
  kmer_num_all <- totalization %>% filter(q_value < 0.05) %>% group_by(ID) %>% summarise(sig_kmer_num = n()) 
  kmer_num_all2 <- kmer_num_all %>% mutate(filter =  rep("All", nrow(kmer_num_all)))
  #kmer_num_soft <- totalization %>% filter(ID %in% ID_soft) %>% filter(q_value < 0.05) %>% group_by(ID) %>% summarise(sig_kmer_num = n())
  #kmer_num_soft2 <- kmer_num_soft %>% mutate(filter =  rep("Soft", nrow(kmer_num_soft)))
  kmer_num_hard <- totalization %>% filter(ID %in% ID_hard) %>% filter(q_value < 0.05) %>% group_by(ID) %>% summarise(sig_kmer_num = n())
  kmer_num_hard2 <- kmer_num_hard %>% mutate(filter =  rep("Hard", nrow(kmer_num_hard)))
  
  #df_plot <- rbind(kmer_num_all2, kmer_num_soft2)
  #df_plot <- rbind(df_plot, kmer_num_hard2)
  df_plot <- rbind(kmer_num_all2, kmer_num_hard2)
  
  #df_plot2 <- transform(df_plot, filter= factor(filter, levels = c("All", "Soft", "Hard")))
  df_plot2 <- transform(df_plot, filter= factor(filter, levels = c("All", "Hard")))
  p <- df_plot2 %>% ggplot(aes(x = filter, y = sig_kmer_num, fill = filter)) +
    geom_boxplot() +
    ylab("number of significant k-mer") +
    scale_fill_manual(values = c("#F8766D","#619CFF"))+ #added
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        aspect.ratio = 1,
        legend.position = "none"
    )
  
  
  return(p)
  
  
}