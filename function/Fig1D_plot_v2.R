Fig1D_plot <- function(TF_list, load, filter){
  
  #############################################
  ############# New Fig. 1D ##################
  #############################################
  
  library(tidyverse)
  library(pROC)
  library(colorspace)
  library(RColorBrewer)
  
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
  
  AUC_list <- list()
  
  for(i in 1:length(TF_list)){
    target_TF <- TF_list[i]
    print(target_TF)
    
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
        
      }#ifelseのifの終わり
    }
  }
  
  #print(paste0("AUC = ", AUC_list[[target_TF]]))
  
  tib <- tibble()
  for (i in 1:length(TF_list)) {
    target_TF <- TF_list[i]
    df <- tibble(tf = rep(target_TF, length(AUC_list[[target_TF]])), 
                 auc = AUC_list[[target_TF]]
                 )
    if(i == 1){
      tib <- df
    }else{
      tib <- rbind(tib, df)
    }
  }
  
  tib$auc <- as.numeric(tib$auc)
  p <- tib %>% ggplot(aes(x = reorder(tf, -auc), y = auc, fill = tf)) +
    geom_violin() +
    xlab("TF")+
    ylab("AUC")+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=12,face="bold", angle = 45, hjust = 1),
          axis.text.y =element_text(size=15,face="bold"),
          axis.title=element_text(size=15,face="bold"),
          aspect.ratio = 1
    )
  
  
  # Add shuffle data ----------
  AUC_list_shuf <- list()
  for(i in 1:length(TF_list)){
    target_TF <- TF_list[i]
    print(target_TF)
    
    # filter target TF table
    target_MOCCS <- hg38_selected %>%
      select(ID, Antigen, Cell_type_class, kmer, MOCCS2score) %>%
      filter(Antigen == target_TF) 
    MOCCS2score <- sample(target_MOCCS$MOCCS2score, nrow(target_MOCCS))
    target_MOCCS_shuf <- target_MOCCS %>% select(-MOCCS2score) %>% mutate(MOCCS2score = MOCCS2score)
    
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
    color_list <- rep(c(brewer.pal(10,"Spectral"),brewer.pal(10,"BrBG")), 15)
    
    for (z in seq_along(sample_list)) {
      
      # get PWM top PWMscore k-mer (using the number of k-mer in target_MOCCS k-mer)
      target_sample <- sample_list[z]
      target_MOCCS_kmer <- target_MOCCS_shuf %>%
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
        AUC_list_shuf[[target_TF]][[z]] <- ROC$auc
        
      }#ifelseのifの終わり
    }
  }
  
  tib2 <- tibble()
  for (i in 1:length(TF_list)) {
    target_TF <- TF_list[i]
    print(target_TF)
    df <- tibble(tf = rep(target_TF, length(unlist(AUC_list_shuf[[target_TF]]))), 
                 auc = unlist(AUC_list_shuf[[target_TF]])
    )
    if(i == 1){
      tib2 <- df
    }else{
      tib2 <- rbind(tib2, df)
    }
  }
  
  tib <- tib %>% mutate(label = "original")
  tib2 <- tib2 %>% mutate(label = "shuffle")
  tib3 <- rbind(tib, tib2)
  
  #tib3$auc <- unlist(tib3$auc)
  p2 <- tib3 %>% ggplot(aes(x = reorder(tf, -auc), y = auc, fill = label)) +
    geom_violin() +
    scale_fill_manual(values = c("red", "gray")) +
    xlab("TF")+
    ylab("AUC")+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          #legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=12,face="bold"),
          axis.text.y =element_text(size=15,face="bold"),
          axis.title=element_text(size=15,face="bold"),
          aspect.ratio = 1
    )
  
  return(p2)
}
