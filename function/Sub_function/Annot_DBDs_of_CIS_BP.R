Annot_DBDs_of_CIS_BP <- function(df_raw,
                                 TF_info_dir = "~/MOCCS_paper_public/annotation/CIS-BP/Homo_sapiens_2021_04_21_1_47_am/TF_Information.txt"){
  
  library("tidyverse")
  
  df_raw[, c("ID", "Antigen", "Cell_type", "Cell_type_class")] %>%
    dplyr::distinct(ID, .keep_all = TRUE) -> df5
  
  input_TF_symbol_list <- as.character(df5$Antigen)
  
  ## Import db data
  TF_info <- read.csv(TF_info_dir, sep = "\t")
  TF_info_symbol <- TF_info$TF_Name
  
  TF_info_symbol_uniq <- TF_info_symbol %>% unique()
  uniq_row_list <- c()
  
  for (uniq_index in 1:length(TF_info_symbol_uniq)){
  
    uniq_row_num <- match(TF_info_symbol_uniq[uniq_index],
                         TF_info_symbol)[1]
    
    uniq_row_list <- append(uniq_row_list, uniq_row_num)
    
  }
  
  TF_info_uniq <- TF_info[uniq_row_list,]
  
  ## gene ID を mutate する
  library(org.Hs.eg.db)
  
  hs <- org.Hs.eg.db
  
  # DB 側
  TF_info_symbol_uniq_ID <- select(hs,
                                   keys = TF_info_uniq$TF_Name,
                                   column = c("ENTREZID", "SYMBOL"),
                                   keytype = "SYMBOL")
  
  TF_info_uniq <- TF_info_uniq %>%
                        mutate("ENTREZID" = TF_info_symbol_uniq_ID$ENTREZID)
  
  # input 側
  input_TF_symbol_list_ID <- select(hs,
                                 keys = input_TF_symbol_list,
                                 column = c("ENTREZID", "SYMBOL"),
                                 keytype = "SYMBOL") %>% as_tibble()
  
  NA_list_in_db <- filter(input_TF_symbol_list_ID, is.na(ENTREZID))
  NA_symbol_list_in_db <- NA_list_in_db[,"SYMBOL"] %>% unique()
  print(paste0("No conveted TFs: ", NA_symbol_list_in_db))
  
  TF_info_uniq_DBDs <- TF_info_uniq$DBDs
  TF_info_uniq_FamilyName <- TF_info_uniq$Family_Name
  
  
  ## Symbol ごとに検索する
  input_TF_DBDs_list_with_NA <- c()
  input_TF_FamilyName_list_with_NA <- c()
  i <- 0
  for (symbol_num_index in 1:length(input_TF_symbol_list)){
  
      input_TF_ID <- input_TF_symbol_list_ID[symbol_num_index, "ENTREZID"]
      info_row_num <- match(input_TF_ID, TF_info_symbol_uniq_ID[,"ENTREZID"])
      
      if (is.na(info_row_num)){
        i <- i + 1
      }
      
      input_TF_DBDs <- TF_info_uniq_DBDs[info_row_num]
      input_TF_FamilyName <- TF_info_uniq_FamilyName[info_row_num]
      
      input_TF_DBDs_list_with_NA <- append(input_TF_DBDs_list_with_NA,
                                           input_TF_DBDs)
      input_TF_FamilyName_list_with_NA <- append(input_TF_FamilyName_list_with_NA,
                                                  input_TF_FamilyName)
      
  }
  
  input_TF_FamilyName_list <- replace_na(input_TF_FamilyName_list_with_NA, "No_annotation")
  input_TF_DBDs_list <- replace_na(input_TF_DBDs_list_with_NA, "No_annotation")
  
  df_fam <- cbind(df5, input_TF_FamilyName_list)
  colnames(df_fam)[colnames(df_fam) == "input_TF_FamilyName_list"] <- "Family"
  
  return(df_fam)
  
}