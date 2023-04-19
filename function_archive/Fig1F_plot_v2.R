Fig1F_plot <- function(target_TF, annotation_path){
  
  library(tidyverse)
  library(pROC)
  library(RColorBrewer)
  
  ID_hard <- readRDS(paste0(annotation_path, "hg38_hard_filter_ID.rds"))
  #ID_soft <- readRDS(paste0(annotation_path, "ID_soft_filter_hg38.rds"))
  Antigen_list <- read_tsv(paste0(annotation_path, "Antigen_list_hg38.txt"), col_names = FALSE)
  Antigen_list <- Antigen_list$X1 %>% as.character() %>% unique()
  totalization <- readRDS(url("https://figshare.com/ndownloader/files/34065686","rb")) #MOCCSout_hg38_all_qval_annotated.rds

  # filter q value < 0.05 k-merに
  hg38_selected <- totalization %>%
    filter(Cell_type_class != "Unclassified") %>% 
    filter(q_value < 0.05)%>%
    select(ID, Antigen, Cell_type_class, Cell_type,kmer, MOCCS2score)
  
  # filter target TF table
  target_MOCCS <- hg38_selected %>%
    select(ID, Antigen, Cell_type_class, kmer, MOCCS2score) %>%
    filter(Antigen == target_TF) 
  
  # filter target TF PWM table
  PWM_table_all <- readRDS(url("https://figshare.com/ndownloader/files/34065698","rb")) #PWM_likelihood_HOMER.rds
  target_PWM <- PWM_table_all[[target_TF]]
  
  ## calculate per sample and plot
  sample_list <- unique(target_MOCCS$ID)
  color_list <-  brewer.pal(12,"Set3")
  color_list <- rep(color_list, 100000)
  AUC_table <- c()
  
  for (z in seq_along(sample_list)) {
    
    # get PWM top PWMscore k-mer (using the number of k-mer in target_MOCCS k-mer)
    target_sample <- sample_list[z]
    target_MOCCS_kmer <- target_MOCCS %>%
      arrange(desc(MOCCS2score)) %>%
      filter(ID == target_sample) %>%
      select(kmer, MOCCS2score)
    kmer_N <- nrow(target_MOCCS_kmer)
    
    target_PWM_kmer <- target_PWM %>%
      arrange(desc(PWMscore)) 
    
    # PWMsoreは、k-merに被りがあるので、各k-merの最大値を採用する
    target_PWM_kmer <- target_PWM_kmer %>%
      group_by(kmer) %>%
      summarise(max_PWMscore = max(PWMscore)) 
    target_PWM_kmer <- target_PWM_kmer %>%
      arrange(desc(max_PWMscore))
    
    s <- round(kmer_N*10/100)
    target_PWM_kmer <- target_PWM_kmer[1:s,] #PWM top 10%
    
    df <- left_join(target_MOCCS_kmer, target_PWM_kmer, by = "kmer")
    df[is.na(df)] <- 0 #NAを0に置換
    
    #AUROC用の列を準備
    df <- df %>%
      mutate(ROC = ifelse(max_PWMscore == 0, FALSE, TRUE))
    
    if(length(unique(df$ROC)) == 1){
      print("no common k-mer")
    }else{
      
      #roc()でROC曲線のオブジェクトを作成する
      ROC <- roc(ROC ~ MOCCS2score, data = df, ci = TRUE) #Xが連続値のMOCCS2score, YがMOCCSのkmerがPWMに含まれているかどうか  
      target_row <- tibble(ID = target_sample, AUC = as.numeric(ROC$auc))
      
      if(z == 1){
        AUC_table <- target_row
        #plot(ROC,col=colors()[z])
      }else{
        AUC_table <- AUC_table %>% add_row(target_row)
        #plot(ROC, add = TRUE, col=colors()[z])
      }#zのifの終わり
    }#ifelseのifの終わり
  }
  
  AUC_hard <- AUC_table %>% mutate(filter = ifelse(ID %in% ID_hard, "hard", "others")) %>% filter(filter == "hard")
  #AUC_soft <- AUC_table %>% mutate(filter = ifelse(ID %in% ID_soft, "soft", "others")) %>% filter(filter == "soft")
  AUC_table <- AUC_table %>% mutate(filter = rep("all", nrow(AUC_table)))
  AUC_plot <- rbind(AUC_table, AUC_hard)
  #AUC_plot <- rbind(AUC_plot, AUC_soft)
  
  #AUC_plot2 <- transform(AUC_plot, filter= factor(filter, levels = c("all", "soft", "hard")))
  AUC_plot2 <- transform(AUC_plot, filter= factor(filter, levels = c("all", "hard")))
  p <- AUC_plot2 %>% ggplot(aes(x = filter, y = AUC, fill = filter)) +
    geom_violin()+
    geom_jitter(size = 0.5) +
    ggtitle(target_TF) +
    scale_fill_manual(values = c("#F8766D","#619CFF"))+ #added
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_line(colour="gray"),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
          axis.text.y =element_text(size=10,face="bold"),
          axis.title=element_text(size=14,face="bold"),
          aspect.ratio = 1
    )
  return(p)
}