df_p_3_gp <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_3_gp.rds")


TF_list <- unique(union(df_p_3_gp$ID1_Antigen, df_p_3_gp$ID2_Antigen))

target_tf_list <- c()
for (i in 1:length(TF_list)) {
  print(paste0(i, "/", length(TF_list)))
  tgt_tf <- TF_list[i]
  target_df <- df_p_3_gp %>% filter(ID1_Antigen == tgt_tf & ID2_Antigen == tgt_tf)
  tmp1 <- target_df %>% select(ID1, ID1_Cell_type_class)
  tmp2 <- target_df %>% select(ID2, ID2_Cell_type_class)
  colnames(tmp1) <- c("ID", "Cell_type_class")
  colnames(tmp2) <- c("ID", "Cell_type_class")
  tmp3 <- rbind(tmp1, tmp2) %>% distinct()
  rm_CTC <- tmp3 %>% group_by(Cell_type_class) %>% summarise(n = n()) %>% filter(n == 1) %>% .$Cell_type_class
  target_df_old <- target_df
  target_df <- target_df_old %>% filter(!ID1_Cell_type_class %in% rm_CTC & !ID2_Cell_type_class %in% rm_CTC) 
  
  if(nrow(target_df) != 0){
    target_tf_list <- c(target_tf_list, tgt_tf)
  }
}

length(target_tf_list)
saveRDS()