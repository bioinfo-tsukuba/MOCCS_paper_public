Group_df <- function(df_p_3){
  
  s_ant <- as.numeric(df_p_3$ID1_Antigen == df_p_3$ID2_Antigen)

  s_ctc <- as.numeric(df_p_3$ID1_Cell_type_class == df_p_3$ID2_Cell_type_class)
  unk_flag <- df_p_3$ID1_Cell_type_class == "Unclassified" | df_p_3$ID1_Cell_type_class == "Others" | df_p_3$ID2_Cell_type_class == "Unclassified" | df_p_3$ID2_Cell_type_class == "Others"
  s_ctc[unk_flag] <- 0.5

  s_ct <- as.numeric(df_p_3$ID1_Cell_type == df_p_3$ID2_Cell_type)
  unk_flag <- df_p_3$ID1_Cell_type == "Unclassified" | df_p_3$ID2_Cell_type == "Unclassified"
  s_ct[unk_flag] <- 0.5
  
  s_ant_f <- as.numeric(df_p_3$ID1_Family == df_p_3$ID2_Family)
  mult_f_vec <- union(grep(",", df_p_3$ID1_Family), grep(",", df_p_3$ID2_Family))
  for (i in 1:length(mult_f_vec)){
    pair_id <- mult_f_vec[i]
    id1_f <- df_p_3$ID1_Family[pair_id]
    id2_f <- df_p_3$ID2_Family[pair_id]
    id1_f_div <- strsplit(id1_f, ",")[[1]]
    id2_f_div <- strsplit(id2_f, ",")[[1]]
    jd <- c()
    for (j in 1:length(id1_f_div)){
      for (k in 1:length(id2_f_div)){
        if(id1_f_div[j] == id2_f_div[[k]]){
          # print(pair_id)
          s_ant_f[pair_id] <- 1
        }
      }
    }
  }
  unk_flag <- df_p_3$ID1_Family == "No_annotation" | df_p_3$ID2_Family == "No_annotation" | df_p_3$ID1_Family == "Unknown" | df_p_3$ID2_Family == "Unknown"
  s_ant_f[unk_flag] <- 0.5
  
  df_p_3 %>% mutate(s_ctc = s_ctc) %>%
             mutate(s_ant = s_ant) %>%
             mutate(s_ant_f = s_ant_f) %>%
             mutate(s_ct = s_ct) -> df_p_3_gp
  
  return(df_p_3_gp)
  
}