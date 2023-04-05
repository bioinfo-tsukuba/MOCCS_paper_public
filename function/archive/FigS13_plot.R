FigS13_plot <- function(TF_list){
  
  #############################################
  ############# New Fig. S13##################
  #############################################
  
  library(tidyverse)
  library(pROC)
  library(colorspace)
  library(RColorBrewer)
  
  
  
  ID_hard <- readRDS("~/MOCCS_paper_public/data/Fig1/hg38_hard_filter_ID.rds")
  tib <- tibble()
  
  for (x in c(10, 5, 2.5, 1.25)) {
    print(x)
    for (k in 6:8) {
      
      print(k)
      
      # localから読み込む場合
      if(k == 6){
        totalization <- readRDS("~/MOCCS_paper_public/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds")
        totalization2 <- totalization %>% filter(ID %in% ID_hard & Antigen %in% TF_list)
      }else{
        totalization <- readRDS(paste0("~/MOCCS_paper_public/data/Fig1/totalization_CTCF_FOXA1_SPI1_", k, "mer.rds"))
        totalization2 <- totalization %>% filter(ID %in% ID_hard & Antigen %in% TF_list)
        colnames(totalization2) <- c("kmer", "auc", "count", "MOCCS2score", "p_value", "q_value", "ID", 
                                     "Antigen_class", "Antigen", "Cell_type_class", "Cell_type")
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
        if(k == 6){
          #from local repository
          PWM_table_all <- readRDS("~/MOCCS_paper_public/data/Fig1/PWM_likelihood_HOMER.rds")
        }else{
          #from figshare
          PWM_table_all <- readRDS(paste0("~/MOCCS_paper_public/data/Fig1/PWM_kmer_list_", k, "mer.rds"))
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
          target_PWM_kmer <- target_PWM_kmer %>% drop_na(PWMscore) %>%
            group_by(kmer) %>%
            summarise(max_PWMscore = max(PWMscore)) 
          target_PWM_kmer <- target_PWM_kmer %>%
            arrange(desc(max_PWMscore))
          
          #kmer_N <- 4^k
          s <- round(kmer_N*x/100)
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
        }#for (z in seq_along(sample_list))
      }#for(i in 1:length(TF_list))
      
      for (j in 1:length(TF_list)) {
        target_TF <- TF_list[j]
        df <- tibble(tf = rep(target_TF, length(AUC_list[[target_TF]])), 
                     auc = AUC_list[[target_TF]],
                     k = rep(k, length(AUC_list[[target_TF]])),
                     top = rep(x, length(AUC_list[[target_TF]]))
        )
        if(j == 1 & k == 6 & x == 10){
          tib <- df
        }else{
          tib <- rbind(tib, df)
        }
      }#for (j in 1:length(TF_list))
      
    }#for (k in 6:8)
  }#for (x in c(10, 5, 2.5, 1.25))
  
  
  
  
  tib$auc <- as.numeric(tib$auc)
  tib$top <- as.numeric(tib$top)
  
  p_list <- list()
  for (target_TF in TF_list) {
    p <- tib %>% filter(tf == target_TF) %>%
      ggplot(aes(x = factor(top), y = auc, fill = factor(k))) + 
      geom_violin() +
      xlab("Top ratio")+
      ylab("AUROC")+
      ggtitle(target_TF)+
      theme(plot.title = element_text(face="bold",hjust = 0.5), 
            #legend.position = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour="black"),
            axis.text=element_text(size=12,face="bold"),
            axis.text.x =element_text(size=12,face="bold", angle = 45, hjust = 1),
            axis.text.y =element_text(size=15,face="bold"),
            axis.title=element_text(size=15,face="bold"),
            aspect.ratio = 0.2
      )
    plot(p)
    p_list[[target_TF]] <- p
  }
  return(p_list)
}