Fig6B_plot <- function( annotation_path){
  
  library(tidyverse)
  Fig6B_df <- read_csv(paste0(annotation_path, "Fig6B_df.csv"))
  Fig6B_df <- Fig6B_df %>% unite("color", c(peak, significant), sep = "_")
  
  #Fig6B_df$color[is.na(Fig6B_df$color)] 
  Fig6B_df$color <- factor(Fig6B_df$color, levels = c("without_non-sig", "within_non-sig", "within_sig"))
  Fig6B_plot <- Fig6B_df %>% ggplot(aes(x = phenotype, y = snp_num, fill = color)) +
    geom_bar(stat = "identity", position = "fill") +
    #geom_text(aes(label = snp_num), size = 2,position = "stack") +
    scale_fill_manual(values = c("#696969", "#DC143C", "#FF7F50")) +
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_line(colour="gray"),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=6,face="bold", angle = 45, hjust = 1),
          axis.text.y =element_text(size=10,face="bold"),
          axis.title=element_text(size=14,face="bold"),
          #legend.position = 'none',
          legend.title = element_blank(),
          aspect.ratio = 1
    )+
    ggtitle("Number of SNPs") +
    labs(fill = "") +
    coord_flip()
  
  return(Fig6B_plot)
  
  
}