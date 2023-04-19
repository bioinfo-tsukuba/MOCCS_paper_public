Annot_pairs <- function(df_fam){
  
  df_fam[, c("ID", "Antigen", "Cell_type", "Cell_type_class", "Family")] %>%
    dplyr::distinct(ID, .keep_all = TRUE) -> df_fam_2

  id_list <- unique(df_fam_2$ID)
  pair_num <- length(id_list) * (length(id_list) - 1) / 2
  
  df_p_1_colnames <- c("ID1", "ID1_Antigen", "ID1_Cell_type_class", "ID1_Family", "ID1_Cell_type",
                    "ID2", "ID2_Antigen", "ID2_Cell_type_class", "ID2_Family", "ID2_Cell_type")
  df_p_1 <- matrix(0, nrow = pair_num, ncol = length(df_p_1_colnames))
  colnames(df_p_1) <- df_p_1_colnames

  pair_id <- 0
  for (id1_ind in 1:(length(id_list) - 1)){
    id1 <- df_fam_2$ID[id1_ind]
    id1_ant <- as.character(df_fam_2$Antigen[id1_ind])
    id1_ctc <- as.character(df_fam_2$Cell_type_class[id1_ind])
    id1_fam <- as.character(df_fam_2$Family[id1_ind])
    id1_ct <- as.character(df_fam_2$Cell_type[id1_ind])
    
    for (id2_ind in (id1_ind + 1):length(id_list)){
      id2 <- df_fam_2$ID[id2_ind]
      id2_ant <- as.character(df_fam_2$Antigen[id2_ind])
      id2_ctc <- as.character(df_fam_2$Cell_type_class[id2_ind])
      id2_fam <- as.character(df_fam_2$Family[id2_ind])
      id2_ct <- as.character(df_fam_2$Cell_type[id2_ind])
      
      pair_id <- pair_id + 1
      
      if (pair_id %% 100000 == 0){
        print(pair_id)
      }      
      
      # Substitute ID1, ID1_Antigen and ID1_Cell_type_class
      df_p_1[pair_id, "ID1"] <- id1
      df_p_1[pair_id, "ID1_Antigen"] <- id1_ant
      df_p_1[pair_id, "ID1_Cell_type_class"] <- id1_ctc
      df_p_1[pair_id, "ID1_Family"] <- id1_fam
      df_p_1[pair_id, "ID1_Cell_type"] <- id1_ct
      
      # Substitute ID2, ID2_Antigen and ID2_Cell_type_class
      df_p_1[pair_id, "ID2"] <- id2
      df_p_1[pair_id, "ID2_Antigen"] <- id2_ant
      df_p_1[pair_id, "ID2_Cell_type_class"] <- id2_ctc
      df_p_1[pair_id, "ID2_Family"] <- id2_fam
      df_p_1[pair_id, "ID2_Cell_type"] <- id2_ct
      
    }  
      
  }
  
  df_p_1 <- as_tibble(df_p_1)
  
  return(df_p_1)
  
}