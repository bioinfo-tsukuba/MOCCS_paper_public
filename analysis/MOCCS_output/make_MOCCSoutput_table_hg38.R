# rbind MOCCS outputs (only hard filter)
rm (list = ls())
library(tidyverse)
path <-  "/home/s-tahara/MOCCS-DB/WORK/MOCCS/MOCCS-OUTPUT-HG38-6mer-ANALYSIS_v2"

ID_hard <- readRDS("/home/s-tahara/DROMPA/code_hg38_2/hg38_hard_filter_ID.rds")
df_binded_all <- c()
for (i in 1:length(ID_hard)) {
  
  target_ID <- ID_hard[i]
  target_file_path <- paste0(path, "/", target_ID, "_6mer_v2.auc_count.txt") 
  
  if(file.exists(target_file_path) == TRUE){
    df1 <- suppressMessages(read_tsv(target_file_path, col_names=TRUE))
    df2 <- df1 %>% mutate(ID = rep(target_ID, nrow(df1)))
  
    if(i == 1){
      df_binded_all <- df2
    }else{
      df_binded_all <- rbind(df_binded_all, df2)
    }
  }else{
    print("no target file")
  }
}
dim(df_binded_all)
length(unique(df_binded_all$ID))
saveRDS(df_binded_all, "/home/s-tahara/MOCCS-DB/WORK/MOCCS/MOCCS-ANALYSIS/output_hg38/MOCCSout_hg38_hard_filter_rbinded.rds")
MOCCSout_ID <- unique(df_binded_all$ID) %>% as.character()
retry_MOCCS_ID <- setdiff(ID_hard, MOCCSout_ID)

# add annotation
experimentList_tab4 <- readRDS("/home/s-tahara/DROMPA/code_hg38_2/experimentList_tab4.rds")
annotation_hg38_hard <- experimentList_tab4 %>% filter(ID %in% ID_hard) %>% distinct() %>% select(ID, Antigen_class, Antigen, Cell_type_class, Cell_type)
df_binded_all_annotated <- df_binded_all %>% left_join(annotation_hg38_hard)
saveRDS(df_binded_all_annotated, "/home/s-tahara/MOCCS-DB/WORK/MOCCS/MOCCS-ANALYSIS/output_hg38/MOCCSout_hg38_hard_filter_annotated.rds")

# calculate q value
library(pbapply)
list_files <- list.files(path,recursive=T,pattern = "_6mer_v2.auc_count.txt")
file_N <- length(list_files)
W <- 350

df_q_value <- bind_rows(
  pblapply(1:file_N, function(file_num){
    
    # ファイルを読み込む
    target_file <- paste(path,list_files[file_num], sep="/") #このファイルリストの番号と行番号が一致
    df1 <- suppressMessages(read_tsv(target_file, col_names=TRUE))
    
    if(nrow(df1) != 0){
      # kmerごと(行ごと)に p valueを計算する
      p_list <- lapply(1:nrow(df1), function(y){
        target_AUC <- as.numeric(df1[y,2])
        target_kmer_count <- as.numeric(df1[y,3])
        target_p <- 1-pnorm(target_AUC, mean = 0, sd = sqrt(W^2/12/target_kmer_count))
        return(target_p)
      })
      
      # sampleごとに多重検定補正
      p_list <- unlist(p_list)
      q_list <- p.adjust(p_list)
      
      # IDとpvalueとqvalueを足したtableにして返す 
      pre_ID <- gsub("_6mer_v2.auc_count.txt", "", target_file)
      ID <- gsub("/home/s-tahara/MOCCS-DB/WORK/MOCCS/MOCCS-OUTPUT-HG38-6mer-ANALYSIS_v2/", "", pre_ID)
      df2 <- df1 %>% mutate(p_value = p_list, q_value = q_list, ID = rep(ID, nrow(df1)))
    } #if
    return(df2)
  }) #lapply
) #bind_rows
dim(df_q_value)
head(df_q_value)
length(unique(df_q_value$ID))
saveRDS(df_q_value, "/home/s-tahara/MOCCS-DB/WORK/MOCCS/MOCCS-ANALYSIS/output_hg38/MOCCSout_hg38_all_qval.rds")
