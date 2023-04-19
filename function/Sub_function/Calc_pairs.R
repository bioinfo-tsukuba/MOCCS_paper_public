Calc_pairs <- function(df_raw){
  
  library(dplyr)
  
  df_raw_sg <- df_raw
  df_raw_sg$MOCCS2score[df_raw_sg$q_value >= 0.05] <- 0
  
  df_raw_sg[, c("kmer", "MOCCS2score", "ID")] %>%
    dplyr::group_by(ID) %>%
    dplyr::distinct(kmer, .keep_all = TRUE) %>%
    tidyr::pivot_wider(names_from = "kmer", values_from = "MOCCS2score" ) -> df_raw_2
  df_raw_2[is.na(df_raw_2)] <- 0
  
  id_list <- unique(df_raw_2$ID)
  pair_num <- length(id_list) * (length(id_list) - 1) / 2
  
  df_p_2_colnames <- c("ID1", "ID2", "k_sim_1", "k_sim_2")
  df_p_2 <- matrix(0, nrow = pair_num, ncol = length(df_p_2_colnames))
  colnames(df_p_2) <- df_p_2_colnames

  pair_id <- 0
  for (id1_ind in 1:(length(id_list) - 1)){
  # for (id1_ind in 1:100){
    id1 <- id_list[id1_ind]
    id1_score <- unlist(df_raw_2[id1_ind, colnames(df_raw_2) != "ID"])

    for (id2_ind in (id1_ind + 1):length(id_list)){
    # for (id2_ind in (id1_ind + 1):100){
      id2 <- id_list[id2_ind]
      id2_score <- unlist(df_raw_2[id2_ind, colnames(df_raw_2) != "ID"])
      
      pair_id <- pair_id + 1

      if (pair_id %% 10000 == 0){
        print(pair_id)
      }      
      
      # Substitute ID1, ID1_Antigen and ID1_Cell_type_class
      df_p_2[pair_id, "ID1"] <- id1

      # Substitute ID2, ID2_Antigen and ID2_Cell_type_class
      df_p_2[pair_id, "ID2"] <- id2

      # Calculate k-sim 1/2
      id1_sig_k_mer <- names(id1_score[abs(id1_score) > 0])
      id2_sig_k_mer <- names(id2_score[abs(id2_score) > 0])
      
      k_sim_1 <- length(intersect(id1_sig_k_mer, id2_sig_k_mer)) / length(union(id1_sig_k_mer, id2_sig_k_mer))
      k_sim_2 <- cor(id1_score, id2_score)
      
      # Substitute k-sim 1/2
      df_p_2[pair_id, "k_sim_1"] <- k_sim_1
      df_p_2[pair_id, "k_sim_2"] <- k_sim_2
      
    }
    
  }
  
  # Pre-processing
  df_p_2[,"k_sim_2"][is.na(df_p_2[,"k_sim_2"])] <- 0
  df_p_2 <- as_tibble(df_p_2)
  df_p_2[,"k_sim_1"] <- as.double(df_p_2$k_sim_1)
  df_p_2[,"k_sim_2"] <- as.double(df_p_2$k_sim_2)
  
  return(df_p_2)
  
}