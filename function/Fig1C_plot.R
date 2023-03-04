Fig1C_plot <- function(target_TF, load, filter){
  
  library(tidyverse)
  library(pROC)
  library(colorspace)
  
  # localから読み込む場合
  if(load == "local"){
    totalization <- readRDS("~/MOCCS_paper_public/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds")
  }else{
    # figshareから読み込む場合
    totalization <- readRDS(url("https://figshare.com/ndownloader/files/34065686","rb")) #MOCCSout_hg38_all_qval_annotated.rds
  }
  if(filter == "hard"){
    ID_hard <- readRDS("~/MOCCS_paper_public/data/Fig1/hg38_hard_filter_ID.rds")
    totalization2 <- totalization %>% filter(ID %in% ID_hard)
  }else if(filter == "soft"){
    ID_soft <- readRDS("~/MOCCS_paper_public/data/Fig1/ID_soft_filter_hg38.rds")
    totalization2 <- totalization %>% filter(ID %in% ID_soft)
  }else{
    totalization2 <- totalization
  }
  
  # filter q value < 0.05 k-merに
  hg38_selected <- totalization2 %>%
    filter(Cell_type_class != "Unclassified") %>% 
    filter(q_value < 0.05)%>%
    select(ID, Antigen, Cell_type_class, Cell_type,kmer, MOCCS2score)
  
  # filter target TF table
  target_MOCCS <- hg38_selected %>%
    select(ID, Antigen, Cell_type_class, kmer, MOCCS2score) %>%
    filter(Antigen == target_TF) 
  
  # filter target TF PWM table
  if(load == "local"){
    #from local repository
    PWM_table_all <- readRDS("~/MOCCS_paper_public/data/Fig1/PWM_likelihood_HOMER.rds")
  }else{
    #from figshare
    PWM_table_all <- readRDS(url("https://figshare.com/ndownloader/files/34065698","rb"))
  }
  target_PWM <- PWM_table_all[[target_TF]]
  
  ## calculate per sample and plot
  sample_list <- unique(target_MOCCS$ID)
  AUC_list <- list()
  color_list <- rep(c(brewer.pal(10,"Spectral"),brewer.pal(10,"BrBG")), 15)
  
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
      AUC_list[[target_TF]][[z]] <- "no common kmer"
    }else{
      
      #roc()でROC曲線のオブジェクトを作成する
      ROC <- roc(ROC ~ MOCCS2score, data = df, ci = TRUE) #Xが連続値のMOCCS2score, YがMOCCSのkmerがPWMに含まれているかどうか  
      AUC_list[[target_TF]][[z]] <- ROC$auc
      
      #png(paste0("~/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig1/Fig1C_", target_TF, ".png" ))
      if(z == 1){
        #plot(ROC,col=colors()[z])
        plot(ROC,col=color_list[z])
      }else{
        #plot(ROC, add = TRUE, col=colors()[z])
        plot(ROC, add = TRUE, col=color_list[z])
      }#zのifの終わり
      #dev.off()
    }#ifelseのifの終わり
  }
  
  print(paste0("AUC = ", AUC_list[[target_TF]]))
  return(AUC_list[[target_TF]])
}