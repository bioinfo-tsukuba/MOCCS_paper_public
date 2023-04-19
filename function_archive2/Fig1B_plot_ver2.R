Fig1B_plot_ver2 <- function(){
  
  library(tidyverse)
  library(RColorBrewer)
  library(colorspace)
  
  annotation_path <- "~/MOCCS_paper_public/data/Fig1/"
  annotation_all <- readRDS(url("https://figshare.com/ndownloader/files/34065671","rb")) #experimentList_tab4.rds
  ID_hard <- readRDS(paste0(annotation_path, "hg38_hard_filter_ID.rds"))
  Antigen_list <- read_tsv(paste0(annotation_path, "Antigen_list_hg38.txt"), col_names = FALSE)
  Antigen_list <- Antigen_list$X1 %>% as.character() %>% unique()
  
  annotation_hg38 <- annotation_all %>% filter(Genome == "hg38" & Antigen_class == "TFs and others" & ID %in% Antigen_list) %>% select(ID, Antigen_class, Antigen, Cell_type_class, Cell_type) %>% distinct()
  df_all <- annotation_hg38 %>% mutate(filter = "All")
  df_hard <- annotation_hg38 %>% filter(ID %in% ID_hard) %>% mutate(filter = "Hard")
  df_plot <- rbind(df_all, df_hard)
  
  tmp <- df_hard %>% group_by(Antigen) %>% summarise(n = n()) %>% arrange(desc(n))
  target_tf_list <- tmp$Antigen[1:20]%>% as.character()
  #target_tf_list <- tmp[1:20, 1] %>% .$Antigen %>% as.character()
  df_hard2 <- df_hard %>% filter(Antigen %in% target_tf_list)  %>% left_join(tmp , by = "Antigen")
  color_list <- c(brewer.pal(10,"Spectral"),brewer.pal(10,"BrBG"))
  
  p1 <- df_hard2 %>% ggplot(aes(x = reorder(Antigen, -n), fill = Cell_type_class))+
    geom_bar(width = 0.8)+
    xlab("TF") +
    scale_fill_manual(values = color_list) +
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
          axis.text.y =element_text(size=10,face="bold"),
          axis.title=element_text(size=14,face="bold"),
          #legend.position = "none",
          aspect.ratio = 0.8
    )
  
  
  tmp2 <- df_hard %>% group_by(Cell_type_class) %>% summarise(n = n()) %>% arrange(desc(n))
  df_hard3 <- df_hard %>% left_join(tmp2, by = "Cell_type_class")
  color_list <- c(brewer.pal(10,"Spectral"),brewer.pal(10,"BrBG"))
  
  p2 <- df_hard3 %>% ggplot(aes(x = reorder(Cell_type_class, -n), fill = Cell_type_class))+
    geom_bar(width = 0.8)+
    xlab("Cell type class") +
    scale_fill_manual(values = color_list) +
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
          axis.text.y =element_text(size=10,face="bold"),
          axis.title=element_text(size=14,face="bold"),
          legend.position = "none",
          aspect.ratio = 0.7
    )
  
  return(list(p1, p2))
}