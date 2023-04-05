###### plot -----
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
#ant_df <- read_tsv("~/MOCCS_paper_public/data/TableS1/Supplementary_table_S1_parts_7.tsv")
#target_ant <- ant_df %>% filter(ctc_depedent == "Y" | ctc_depedent == "N") %>% .$TF %>% unique()
exec_Large_Heatmap <- function(){
  library(tidyverse)
  library(ComplexHeatmap)
  library(colorspace)
  large_mat_wider_tib3 <- readRDS("~/MOCCS_paper_public/plot/tmp/large_mat_wider_tib3.rds")
  rowname <- large_mat_wider_tib3$ID1
  large_mat <- large_mat_wider_tib3 %>% select(-ID1) %>% as.matrix()
  rownames(large_mat) <- rowname
  
  df_heatmap_2 <- large_mat
  df_heatmap_2[is.na(df_heatmap_2)] <- 0
  dim(df_heatmap_2)
  View(df_heatmap_2[1:100,1:100])
  
  # Complexheatmap ------
  # Antigenごとに色を分けてplot
  ant_list <- readRDS("~/MOCCS_paper_public/plot/tmp/ant_list.rds")
  annot <- data.frame(Antigen = ant_list)
  
  col_20_2 <- sample(rainbow(length(unique(as.character(annot$Antigen)))), length(unique(as.character(annot$Antigen))))
  if (nrow(df_heatmap_2) == 0){
    return(NULL)
  } else {
    #mtxt <- paste0(target_ant, "")
    #col_20_2 <- qualitative_hcl(length(unique(annot$Antigen)), c = 50, l = 70)
    
    png(paste0("~/MOCCS_paper_public/plot/Fig1/heatmap_all_jaccard.png"),
    #png(paste0("~/MOCCS_paper_public/plot/tmp/heatmap_all_jaccard.png"),
        width = 1800, height = 1800)
    
    #pdf(paste0("~/MOCCS_paper_public/plot/tmp/heatmap_all_jaccard.pdf"))
    NMF::aheatmap(df_heatmap_2, Rowv = NA, Colv = NA,
                  labRow = NA, labCol = NA,
                  annColors = list(col_20_2),
                  annCol = list(Antigen = as.character(annot$Antigen)), 
                  #annRow = list(as.character(annot$Antigen)),
                  #main = mtxt,
                  color = "Reds",
                  #annLegend = FALSE,
                  cexRow = 0,
                  cexCol = 0,
                  fontsize = 15
    )
    dev.off()
    
  } 
  
  
  if (nrow(df_heatmap_2) == 0){
    return(NULL)
  } else {
    png(paste0("~/MOCCS_paper_public/plot/Fig1/heatmap_all_jaccard_rowlable.png"),
    #png(paste0("~/MOCCS_paper_public/plot/tmp/heatmap_all_jaccard_rowlable.png"),
        width = 1800, height = 1800)
    
    #pdf(paste0("~/MOCCS_paper_public/plot/tmp/heatmap_all_jaccard_rowlabel.pdf"))
    NMF::aheatmap(df_heatmap_2, Rowv = NA, Colv = NA,
                  labRow = NA, labCol = NA,
                  annColors = list(col_20_2),
                  #annCol = list(Antigen = as.character(annot$Antigen)), 
                  annRow = list(Antigen = as.character(annot$Antigen)),
                  #main = mtxt,
                  color = "Reds",
                  #annLegend = FALSE,
                  cexRow = 0,
                  cexCol = 0,
                  fontsize = 15
    )
    dev.off()
    
  } 
  
  
  
  # heatmap.2 ------
  library(gplots)
  col_69 <- sample(rainbow(length(unique(as.character(annot$Antigen)))),   length(unique(as.character(annot$Antigen))))
  #col_69 <- rainbow(length(unique(as.character(annot$Antigen))))
  #col_69 <- qualitative_hcl(length(unique(annot$Antigen)), c = 50, l = 70)
  annot_tf_list <- unique(annot$Antigen) %>% as.character()
  annot_color_list <- c()
  for(i in 1:length(unique(annot$Antigen))){
    tgt_ant <- annot_tf_list[i]
    tmp <- annot %>% as_tibble() %>% filter(Antigen == tgt_ant)
    tgt_ant_num <- nrow(tmp)
    annot_color_list <- c(annot_color_list, rep(col_69[i], tgt_ant_num))
  }
  annot2 <- data.frame(Antigen = ant_list, color = annot_color_list)
  
  
  library(RColorBrewer)
  Colors= brewer.pal(5,"Spectral") %>% rev()
  
  if (nrow(df_heatmap_2) == 0){
    return(NULL)
  } else {
    col_20_2 <- sample(rainbow(length(unique(as.character(annot$Antigen)))),   length(unique(as.character(annot$Antigen))))
    
    png(paste0("~/MOCCS_paper_public/plot/Fig1/heatmap_all_jaccard_v2.png"),
    #png(paste0("~/MOCCS_paper_public/plot/tmp/heatmap_all_jaccard_v2.png"),
        width = 3300, height = 3300)
    
    #pdf(paste0("~/MOCCS_paper_public/plot/tmp/heatmap_all_jaccard.pdf"))
    heatmap.2(df_heatmap_2, 
              Rowv = FALSE, Colv = FALSE,
              #col=Colors, 
              #ColSideColors= annot2$color,
              #RowSideColors = annot2$color,
              labRow = annot$Antigen,
              labCol = annot$Antigen,
              #colRow = unique(annot2$Antigen),
              #colCol= annot2$Antigen,
              keysize = 0.4,
              trace = "none"
    )
    dev.off()
    
  } 
  
}



# TF family ごとに色を分けてplot