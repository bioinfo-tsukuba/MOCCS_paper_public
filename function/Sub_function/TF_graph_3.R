TF_graph_3 <- function(df_p_3){ # k-sim 1
  
  library(ggraph)
  library(tidyverse)
  library(igraph)
  
  df_p_3[is.na(df_p_3)] <- 0
  
  # receiver_tf_list <- readRDS(paste0("~/MOCCS_paper_public/results/sig_symbol_list.rds"))
  receiver_tf_list <- c("JUN", "GATA2")
  
  sender_tf_list <- union(unique(df_p_3$ID1_Antigen), unique(df_p_3$ID2_Antigen))
  
  for (i in 1:length(receiver_tf_list)){
    
    at_1 <- receiver_tf_list[i]
    target_flag_1 <- df_p_3$ID1_Antigen == at_1 | df_p_3$ID2_Antigen == at_1
    df_at_1 <- df_p_3[target_flag_1, ]
    ctc_list <- union(df_at_1$ID1_Cell_type_class, df_at_1$ID2_Cell_type_class)
    
    k_sim_1_mat <- matrix(0, nrow = length(ctc_list), ncol = length(sender_tf_list))
    rownames(k_sim_1_mat) <- ctc_list
    colnames(k_sim_1_mat) <- sender_tf_list
    
    for (ctc_ind in 1:length(ctc_list)){

      ctc_1 <- ctc_list[ctc_ind]

      for (j in 1:length(sender_tf_list)){
        
        at_2 <- sender_tf_list[j]
        target_flag_2 <- df_p_3$ID1_Antigen == at_1 & df_p_3$ID2_Antigen == at_2 & df_p_3$ID1_Cell_type_class == ctc_1 & df_p_3$ID2_Cell_type_class == ctc_1
        target_flag_3 <- df_p_3$ID1_Antigen == at_2 & df_p_3$ID2_Antigen == at_1 & df_p_3$ID1_Cell_type_class == ctc_1 & df_p_3$ID2_Cell_type_class == ctc_1
        target_flag_4 <- target_flag_2 | target_flag_3
        df_pair <- df_p_3[target_flag_4, ]
        
        print(paste0(at_1, " - ", ctc_1, " - ",  at_2))
        
        if (nrow(df_pair) > 0){
          
          k_sim_1_mat[ctc_ind, j] <- mean(df_pair$k_sim_1)
          # print(mean(df_pair$k_sim_1))
          
        } else {
          
          k_sim_1_mat[ctc_ind, j] <- NA

        }
        
      }
      
    }

    saveRDS(k_sim_1_mat, paste0("~/MOCCS_paper_public/data/Fig3/obj/k_sim_jaccard_mat/k_sim_jaccard_mat_", at_1, ".rds"))
    
  }
  
  
  return()
  
}