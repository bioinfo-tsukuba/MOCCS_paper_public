Get_ans <- function(df_raw){
  
  df_raw_sg <- df_raw
  df_raw_sg$MOCCS2score[df_raw_sg$q_value >= 0.05] <- 0
  
  df_raw_sg[, c("kmer", "MOCCS2score", "ID", "Antigen", "Cell_type_class")] %>%
    dplyr::group_by(ID) %>%
    dplyr::distinct(kmer, .keep_all = TRUE) %>%
    tidyr::pivot_wider(names_from = "kmer", values_from = "MOCCS2score" ) -> df_raw_2
  df_raw_2[is.na(df_raw_2)] <- 0
  
  df_raw_3 <- df_raw_2[, 4:ncol(df_raw_2)]
  
  all_ns <- df_raw_2$ID[rowSums(df_raw_3) == 0]
  
  return(all_ns)

}