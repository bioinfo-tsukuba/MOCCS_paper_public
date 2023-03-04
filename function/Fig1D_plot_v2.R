Fig1D_plot <- function(TF_list, annotation_path){
  
  library(tidyverse)
  library(RColorBrewer)
  library(Biostrings)
  
  totalization <- readRDS(url("https://figshare.com/ndownloader/files/34065686","rb")) #MOCCSout_hg38_all_qval_annotated.rds
  ID_hard <- readRDS(paste0(annotation_path, "hg38_hard_filter_ID.rds"))
  
  # filter target TF PWM table
  PWM_table_all <- readRDS(url("https://figshare.com/ndownloader/files/34065698","rb")) #PWM_likelihood_HOMER.rds
  
  AUC_list <- list()
  for (i in 1:length(TF_list)) {
    target_TF <- TF_list[i]
    print(target_TF)
    target_PWM <- PWM_table_all[[target_TF]]
    target_PWM2 <- target_PWM %>% group_by(kmer) %>% summarise(max_PWMscore = max(PWMscore))
    
    # join MOCCS output and PWM table
    target_MOCCS <- totalization %>% filter(ID %in% ID_hard) %>% filter(Antigen == target_TF)
    
    # PWMのk-merをMOCCSの片方のk-merに絞る (相補鎖を考慮する)
    target_kmer <- unique(target_MOCCS$kmer)
    non_target_kmer <- setdiff(unique(target_PWM2$kmer), target_kmer)
    target_PWM3 <- target_PWM2 %>% mutate(kmer_class = ifelse(kmer %in% target_kmer, "MOCCS_kmer", "other_kmer"))
    
    for (r in 1:nrow(target_PWM3)) {
      target_row <- target_PWM3[r,]
      if(target_row$kmer_class == "other_kmer"){
        target_row$kmer <- reverseComplement(DNAStringSet(target_row$kmer)) %>% as.character()
      }
      if(r == 1){
        target_PWM4 <- target_row
      }else{
        target_PWM4 <- target_PWM4 %>% add_row(target_row)
      }
    }
    target_PWM5 <- target_PWM4 %>% group_by(kmer) %>% summarise(max_PWMscore = max(max_PWMscore))
    
    df_join1 <- target_MOCCS %>% left_join(target_PWM5, by = "kmer")
    df_join1[is.na(df_join1)] <- 0 #NAを0に置換
    
    ## calculate per sample and plot (PWMscoreが連続値になる)
    sample_list <- unique(df_join1$ID)
    color_list <- rep(c(brewer.pal(10,"Spectral"),brewer.pal(10,"BrBG")), 15)
    
    for (z in seq_along(sample_list)) {
      
      # get PWM top PWMscore k-mer (using the number of k-mer in target_MOCCS k-mer)
      target_sample <- sample_list[z]
      target_PWM_kmer <- df_join1 %>%
        arrange(desc(max_PWMscore)) %>%
        filter(ID == target_sample) %>%
        select(kmer, max_PWMscore)
      
      target_MOCCS_kmer <- target_MOCCS %>%
        arrange(desc(MOCCS2score)) %>%
        filter(ID == target_sample) %>% 
        select(kmer,  MOCCS2score, q_value)
      kmer_N <- nrow(target_MOCCS_kmer)
      
      target_MOCCS_kmer <- target_MOCCS_kmer %>% filter(q_value < 0.05) #MOCCS2score significant
      
      df <- left_join(target_PWM_kmer, target_MOCCS_kmer, by = "kmer")
      df[is.na(df)] <- 0 #NAを0に置換
      
      #AUROC用の列を準備
      df <- df %>%
        mutate(ROC = ifelse(MOCCS2score == 0, FALSE, TRUE))
      
      if(length(unique(df$ROC)) == 1){
        print("no common k-mer")
        AUC_list[[target_TF]][[z]] <- "no common kmer"
      }else{
        
        #roc()でROC曲線のオブジェクトを作成する
        ROC <- roc(ROC ~ max_PWMscore, data = df, ci = TRUE) #Xが連続値のPWMscore, YがMOCCSのkmerがsignificantかどうか  
        AUC_list[[target_TF]][[z]] <- ROC$auc
        
      }#ifelseのifの終わり
    } #for (z in seq_along(sample_list)) 
  } #for (i in 1:length(TF_list)) 
  
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
  for (i in 1:length(TF_list)) {
    target_TF <- TF_list[i]
    print(target_TF)
    target_PWM <- PWM_table_all[[target_TF]]
    target_PWM2 <- target_PWM %>% group_by(kmer) %>% summarise(max_PWMscore = max(PWMscore))
    
    # join MOCCS output and PWM table
    target_MOCCS <- totalization %>% filter(ID %in% ID_hard) %>% filter(Antigen == target_TF)
    
    # PWMのk-merをMOCCSの片方のk-merに絞る (相補鎖を考慮する)
    target_kmer <- unique(target_MOCCS$kmer)
    non_target_kmer <- setdiff(unique(target_PWM2$kmer), target_kmer)
    target_PWM3 <- target_PWM2 %>% mutate(kmer_class = ifelse(kmer %in% target_kmer, "MOCCS_kmer", "other_kmer"))
    
    for (r in 1:nrow(target_PWM3)) {
      target_row <- target_PWM3[r,]
      if(target_row$kmer_class == "other_kmer"){
        target_row$kmer <- reverseComplement(DNAStringSet(target_row$kmer)) %>% as.character()
      }
      if(r == 1){
        target_PWM4 <- target_row
      }else{
        target_PWM4 <- target_PWM4 %>% add_row(target_row)
      }
    }
    target_PWM5 <- target_PWM4 %>% group_by(kmer) %>% summarise(max_PWMscore = max(max_PWMscore))
    
    df_join1 <- target_MOCCS %>% left_join(target_PWM5, by = "kmer")
    df_join1[is.na(df_join1)] <- 0 #NAを0に置換
    
    ## calculate per sample and plot (PWMscoreが連続値になる)
    sample_list <- unique(df_join1$ID)
    color_list <- rep(c(brewer.pal(10,"Spectral"),brewer.pal(10,"BrBG")), 15)
    
    for (z in seq_along(sample_list)) {
      
      # get PWM top PWMscore k-mer (using the number of k-mer in target_MOCCS k-mer)
      target_sample <- sample_list[z]
      target_PWM_kmer <- df_join1 %>%
        arrange(desc(max_PWMscore)) %>%
        filter(ID == target_sample) %>%
        select(kmer, max_PWMscore)
      PWM_shuf <-sample(target_PWM_kmer$max_PWMscore, nrow(target_PWM_kmer))
      target_PWM_kmer_shuf <- target_PWM_kmer %>% select(-max_PWMscore) %>% mutate(max_PWMscore = PWM_shuf)
      
      target_MOCCS_kmer <- target_MOCCS %>%
        arrange(desc(MOCCS2score)) %>%
        filter(ID == target_sample) %>% 
        select(kmer,  MOCCS2score, q_value)
      kmer_N <- nrow(target_MOCCS_kmer)
      
      target_MOCCS_kmer <- target_MOCCS_kmer %>% filter(q_value < 0.05) #MOCCS2score significant
      
      df <- left_join(target_PWM_kmer_shuf, target_MOCCS_kmer, by = "kmer")
      df[is.na(df)] <- 0 #NAを0に置換
      
      #AUROC用の列を準備
      df <- df %>%
        mutate(ROC = ifelse(MOCCS2score == 0, FALSE, TRUE))
      
      if(length(unique(df$ROC)) == 1){
        print("no common k-mer")
        AUC_list_shuf[[target_TF]][[z]] <- "no common kmer"
      }else{
        
        #roc()でROC曲線のオブジェクトを作成する
        ROC <- roc(ROC ~ max_PWMscore, data = df, ci = TRUE) #Xが連続値のPWMscore, YがMOCCSのkmerがsignificantかどうか  
        AUC_list_shuf[[target_TF]][[z]] <- ROC$auc
        
      }#ifelseのifの終わり
    } #for (z in seq_along(sample_list))
  }
  
  tib2 <- tibble()
  for (i in 1:length(TF_list)) {
    target_TF <- TF_list[i]
    print(target_TF)
    df <- tibble(tf = rep(target_TF, length(AUC_list_shuf[[target_TF]])), 
                 auc = AUC_list_shuf[[target_TF]]
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
  
  tib3$auc <- as.numeric(tib3$auc)
  p2 <- tib3 %>% ggplot(aes(x = reorder(tf, -auc), y = auc, fill = label)) +
    geom_violin() +
    scale_x_discrete(limit=c("CTCF", "SPI1", "FOXA1")) +
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