# calculate AUC in CTCF samples (Fig1C) ---------------
target_TF <- "CTCF"
load <- "local" #local or figshare
filter <- "all" #soft or hard or all

library(tidyverse)
library(pROC)
library(colorspace)

# localから読み込む場合
if(load == "local"){
  totalization <- readRDS("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds")
}else{
  # figshareから読み込む場合
  totalization <- readRDS(url("https://figshare.com/ndownloader/files/34065686","rb")) #MOCCSout_hg38_all_qval_annotated.rds
}
if(filter == "hard"){
  ID_hard <- readRDS("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/hg38_hard_filter_ID.rds")
  totalization2 <- totalization %>% filter(ID %in% ID_hard)
}else if(filter == "soft"){
  ID_soft <- readRDS("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/ID_soft_filter_hg38.rds")
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
  PWM_table_all <- readRDS("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/PWM_likelihood_HOMER.rds")
}else{
  #from figshare
  PWM_table_all <- readRDS(url("https://figshare.com/ndownloader/files/34065698","rb"))
}
target_PWM <- PWM_table_all[[target_TF]]

## calculate per sample and plot
sample_list <- unique(target_MOCCS$ID)
AUC_list <- list()
color_list <- qualitative_hcl(length(sample_list), "Dark2")

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
print(paste0("AUC = ", AUC_list[[target_TF]]))

sample_N_all <- length(AUC_list[[target_TF]])
AUC_tb <- tibble(AUC = as.numeric(AUC_list[[target_TF]])) %>% drop_na(AUC)
hist(AUC_tb$AUC, breaks = 30)
AUC_tb_08 <- AUC_tb %>% filter(AUC > 0.8)
AUC_tb_09 <- AUC_tb %>% filter(AUC > 0.9)
ratio_08 <- nrow(AUC_tb_08) / sample_N_all
ratio_08
ratio_09 <- nrow(AUC_tb_09) / sample_N_all
ratio_09


# calculate number of TFs and Cell types after filtering ------------------
annotation_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/"
library(RColorBrewer)

annotation_all <- readRDS(paste0(annotation_path, "experimentList_tab4.rds"))
ID_hard <- readRDS(paste0(annotation_path, "hg38_hard_filter_ID.rds"))
ID_soft <- readRDS(paste0(annotation_path, "ID_soft_filter_hg38.rds"))
Antigen_list <- read_tsv(paste0(annotation_path, "Antigen_list_hg38.txt"), col_names = FALSE)
Antigen_list <- Antigen_list$X1 %>% as.character() %>% unique()

annotation_hg38 <- annotation_all %>% filter(Genome == "hg38" & Antigen_class == "TFs and others" & ID %in% Antigen_list) %>% select(ID, Antigen_class, Antigen, Cell_type_class, Cell_type) %>% distinct()
df_all <- annotation_hg38 %>% mutate(filter = "All")
df_hard <- annotation_hg38 %>% filter(ID %in% ID_hard) %>% mutate(filter = "Hard")
df_soft <- annotation_hg38 %>% filter(ID %in% ID_soft) %>% mutate(filter = "Soft")

# all
length(unique(df_all$Antigen))
length(unique(df_all$Cell_type_class))

# hard
length(unique(df_hard$Antigen))
length(unique(df_hard$Cell_type_class))
length(ID_hard)/10534

# soft
length(unique(df_soft$Antigen))
length(unique(df_soft$Cell_type_class))
length(ID_soft)/10534


# calculate ratio of significant k-mers(Fig1F)-------------
annotation_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/"
totalization_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds"
totalization <- readRDS(totalization_path)
ID_hard <- readRDS(paste0(annotation_path, "hg38_hard_filter_ID.rds"))
ID_soft <- readRDS(paste0(annotation_path, "ID_soft_filter_hg38.rds"))
Antigen_list <- read_tsv(paste0(annotation_path, "Antigen_list_hg38.txt"), col_names = FALSE)
Antigen_list <- Antigen_list$X1 %>% as.character() %>% unique()

kmer_num_all <- totalization  %>% group_by(ID) %>% summarise(kmer_num = n()) 
mean(kmer_num_all$kmer_num)
kmer_num_all_sig <- totalization %>% filter(q_value < 0.05) %>% group_by(ID) %>% summarise(sig_kmer_num = n()) 
mean(kmer_num_all_sig$sig_kmer_num)

kmer_num_soft_sig <- totalization %>% filter(ID %in% ID_soft) %>% filter(q_value < 0.05) %>% group_by(ID) %>% summarise(sig_kmer_num = n())
mean(kmer_num_soft_sig$sig_kmer_num)
mean(kmer_num_soft_sig$sig_kmer_num)/mean(kmer_num_all$kmer_num)

kmer_num_hard_sig <- totalization %>% filter(ID %in% ID_hard) %>% filter(q_value < 0.05) %>% group_by(ID) %>% summarise(sig_kmer_num = n())
mean(kmer_num_hard_sig$sig_kmer_num)
mean(kmer_num_hard_sig$sig_kmer_num)/mean(kmer_num_all$kmer_num)


