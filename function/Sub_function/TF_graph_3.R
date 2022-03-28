TF_graph_3 <- function(df_p_3){ # k-sim 1
  
  library(ggraph)
  library(tidyverse)
  library(igraph)
  
  df_p_3[is.na(df_p_3)] <- 0
  
  receiver_tf_list <- readRDS(paste0("~/MOCCS-DB_paper/results/sig_symbol_list.rds"))[20:42]
  # receiver_tf_list <- union(unique(df_p_3$ID1_Antigen), unique(df_p_3$ID2_Antigen))[1]
  # receiver_tf_list <- c("GATA2")
  # receiver_tf_list <- c("JUNB", "GATA2", "EGR1", "JUND", "MYC", "FOS", "JUN")
  # receiver_tf_list <- union(unique(df_p_3$ID1_Antigen), unique(df_p_3$ID2_Antigen))
  # receiver_tf_list <- readRDS(paste0("~/MOCCS-DB_paper/results/sig_symbol_list.rds"))
  
  
  sender_tf_list <- union(unique(df_p_3$ID1_Antigen), unique(df_p_3$ID2_Antigen))
  
  system(paste0("mkdir ~/MOCCS-DB_paper/plot/Fig3/Fig3X/heatmap_k_sim_1_whole"))
  
  for (i in 1:length(receiver_tf_list)){
    
    at_1 <- receiver_tf_list[i]
    # system(paste0("mkdir ~/MOCCS-DB_paper/plot/Fig3/Fig3X/Fig3X_", at_1, "/"))
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

    saveRDS(k_sim_1_mat, paste0("~/MOCCS-DB_paper/results/k_sim_1_mat/k_sim_1_mat_", at_1, ".rds"))
    # k_sim_1_mat <- readRDS(paste0("~/MOCCS-DB_paper/results/k_sim_1_mat/k_sim_1_mat_", at_1, ".rds"))
    
    if (at_1 == "FOS"){
      
      k_sim_1_mat_2 <- k_sim_1_mat[c(2, 14, 16),]
      k_sim_1_mat_3 <- k_sim_1_mat_2[, !is.na(colSums(k_sim_1_mat_2))]
      
      png(paste0("~/MOCCS-DB_paper/plot/Fig3/Fig3X/heatmap_k_sim_1_", at_1, "_2.png"),
          width = 1600, height = 800)
      NMF::aheatmap(k_sim_1_mat_3,
                    breaks = 0.5,
                    main = paste0("k-sim 2 in ", at_1, " with the selected TFs"),
                    color = "Reds")
      dev.off()
      
    }
    mtxt <- paste0("k-sim 2 in ", at_1, " with the all TFs")
    
    k_sim_1_mat[is.na(k_sim_1_mat)] <- -1
    
    png(paste0("~/MOCCS-DB_paper/plot/Fig3/Fig3X/heatmap_k_sim_1_whole/heatmap_k_sim_1_", at_1, ".png"),
        width = 1600, height = 800)
    NMF::aheatmap(k_sim_1_mat,
                  labCol = NA,
                  breaks = 0.5,
                  main = mtxt)
    dev.off()
  
    
  }
  
  
  return()
  
}