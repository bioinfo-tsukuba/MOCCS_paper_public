Fig1E_plot <- function(annotation_path){
  
  library(tidyverse)
  library(RColorBrewer)
  library(colorspace)
  #library(reshape2)

  #annotation_all <- readRDS(paste0(annotation_path, "experimentList_tab4.rds"))
  annotation_all <- readRDS(url("https://figshare.com/ndownloader/files/34065671","rb")) #experimentList_tab4.rds
  ID_hard <- readRDS(paste0(annotation_path, "hg38_hard_filter_ID.rds"))
  Antigen_list <- read_tsv(paste0(annotation_path, "Antigen_list_hg38.txt"), col_names = FALSE)
  Antigen_list <- Antigen_list$X1 %>% as.character() %>% unique()
  
  annotation_hg38 <- annotation_all %>% filter(Genome == "hg38" & Antigen_class == "TFs and others" & ID %in% Antigen_list) %>% select(ID, Antigen_class, Antigen, Cell_type_class, Cell_type) %>% distinct()
  df_all <- annotation_hg38 %>% mutate(filter = "All")
  df_hard <- annotation_hg38 %>% filter(ID %in% ID_hard) %>% mutate(filter = "Hard")
  
  df_plot <- rbind(df_all, df_hard)
  
  df_plot2 <- transform(df_plot, filter= factor(filter, levels = c("All", "Hard")))
  color_list <- c(brewer.pal(10,"Spectral"),brewer.pal(10,"BrBG"), "gray")
  p1 <- df_plot2 %>% ggplot(aes(x= filter, fill = Cell_type_class)) +
    geom_bar(width = 0.6, position = "fill")+
    scale_fill_manual(values = color_list) +
    xlab("ChIP-seq sample filtering")+
    ylab("Ratio")+
    labs(fill="Cell type class") +
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=15,face="bold"),
          axis.text.y =element_text(size=15,face="bold"),
          axis.title=element_text(size=15,face="bold"),
          aspect.ratio = 1
    )
  
  
  # 少ないサンプル数のTFはothersとしてまとめる
  df_all2  <- df_all %>% filter(Antigen != "Epitope tags" & Antigen != "GFP") %>% group_by(Antigen) %>% summarise(n_sample = n()) %>% arrange(desc(n_sample)) 
  TF_selected_list1 <- df_all2[1:20,]$Antigen %>% as.character()
  df_hard2  <- df_all %>% group_by(Antigen) %>% summarise(n_sample = n()) %>% arrange(desc(n_sample)) 
  TF_selected_list3 <- df_all2[1:20,]$Antigen %>% as.character()
  TF_selected_list <- unique(c(TF_selected_list1, TF_selected_list3))
  df_plot3 <- c()
  for (i in 1:nrow(df_plot)) {
    target_row <- df_plot[i,]
    if(!target_row$Antigen %in% TF_selected_list){
      target_row$Antigen <- "others"
    }
    if(i == 1){
      df_plot3 <- target_row
    }else{
      df_plot3 <- df_plot3 %>% add_row(target_row)
    }
  }
  
  df_plot3$Antigen <- factor(df_plot3$Antigen, levels=rev(c( "ESR1","CTCF","AR","BRD4","FOXA1","EP300","RELA","TP53","MYC","NR3C1","SPI1","SMARCA4","RAD21","EZH2","MED1", "PGR","GATA3","REST","RUNX1","HSF1","others")))
  df_plot3 <- transform(df_plot3, filter= factor(filter, levels = c("All", "Soft", "Hard")))
  color_list2 <- rev(color_list)
  p2 <- df_plot3 %>% ggplot(aes(x=filter, fill = Antigen)) +
    geom_bar(width = 0.6, position = "fill")+
    xlab("ChIP-seq sample filtering")+
    ylab("Ratio")+
    labs(fill="TF") +
    guides(fill = guide_legend(reverse = TRUE))+
    scale_fill_manual(values = color_list2) +
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=15,face="bold"),
          axis.text.y =element_text(size=15,face="bold"),
          axis.title=element_text(size=15,face="bold"),
          aspect.ratio = 1
    )
  
  return(list(p1, p2))
  
}