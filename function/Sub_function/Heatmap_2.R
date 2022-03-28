Heatmap_2 <- function(df_p_3, target_ant_1, target_ant_2, plot_ant){
  
  # Heatmap_2(df_p_3, target_ant_1 = "FOS", target_ant_2 = "JUN", plot_ant = "FOS")
  
  target_ant_1 <- "FOS"
  target_ant_2 <- "JUN"
  
  flag_1 <- df_p_3$ID1_Antigen %in% target_ant_1
  flag_2 <- df_p_3$ID1_Antigen %in% target_ant_2
  flag_3 <- df_p_3$ID2_Antigen %in% target_ant_1
  flag_4 <- df_p_3$ID2_Antigen %in% target_ant_2
  
  target_flag_1 <- flag_1 & flag_3 # ant_1 - ant_1
  target_flag_2 <- flag_2 & flag_4 # ant_2 - ant_2
  target_flag_3 <- flag_1 & flag_4 # ant_1 - ant_2
  target_flag_4 <- flag_2 & flag_3 # ant_2 - ant_1
  
  target_flag_5 <- target_flag_1 | target_flag_2 | target_flag_3 | target_flag_4
  
  target_df <- df_p_3[target_flag_5, ]
  
  if (nrow(target_df) == 0){
    
    return(NULL)
    
  } else {
    
    id_list <- union(unique(target_df$ID1), unique(target_df$ID2))
    ctc_list <- c()
    ant_list <- c()
    
    for (id_ind in 1:length(id_list)){
      target_id <- id_list[id_ind]
      target_ind <- target_df$ID1 == target_id
      target_ctc <- unique(target_df$ID1_Cell_type_class[target_ind])
      target_ant <- unique(target_df$ID1_Antigen[target_ind])
      if (length(target_ctc) == 0){
        target_ind <- target_df$ID2 == target_id
        target_ctc <- unique(target_df$ID2_Cell_type_class[target_ind])
        target_ant <- unique(target_df$ID2_Antigen[target_ind])
      }
      ctc_list <- append(ctc_list, target_ctc)
      ant_list <- append(ant_list, target_ant)
    }
    
    order_ind_1 <- order(ant_list)
    ctc_list_sorted_1 <- ctc_list[order_ind_1]
    ant_list_sorted_1 <- ant_list[order_ind_1]
    id_list_sorted_1 <- id_list[order_ind_1]
    
    ant_1_ind <- ant_list_sorted_1 == target_ant_1
    ant_2_ind <- ant_list_sorted_1 == target_ant_2
    
    ant_1_ctc_list <- ctc_list_sorted_1[ant_1_ind]
    ant_2_ctc_list <- ctc_list_sorted_1[ant_2_ind]
    ant_1_id_list <- id_list_sorted_1[ant_1_ind]
    ant_2_id_list <- id_list_sorted_1[ant_2_ind]
    
    order_ind_2 <- order(ant_1_ctc_list)
    order_ind_3 <- order(ant_2_ctc_list)
    
    ant_1_ctc_list_sorted <- ant_1_ctc_list[order_ind_2]
    ant_2_ctc_list_sorted <- ant_2_ctc_list[order_ind_3]
    ant_1_id_list_sorted <- ant_1_id_list[order_ind_2]
    ant_2_id_list_sorted <- ant_2_id_list[order_ind_3]
    
    id_list_sorted <- append(ant_1_id_list_sorted, ant_2_id_list_sorted)
    ctc_list_sorted <- append(ant_1_ctc_list_sorted, ant_2_ctc_list_sorted)
    
    df_heatmap <- matrix(0, nrow = length(id_list_sorted), ncol = length(id_list_sorted))
    
    for (row_ind in 1:length(id_list_sorted)){
      
      row_id <- id_list_sorted[row_ind]
      
      for (col_ind in 1:length(id_list_sorted)){
        
        col_id <- id_list_sorted[col_ind]
        
        if (row_id == col_id){
          df_heatmap[row_ind, col_ind] <- 1 # same id
        } else {
          target_ind <- (target_df$ID1 == row_id & target_df$ID2 == col_id) | (target_df$ID1 == col_id & target_df$ID2 == row_id)
          df_heatmap[row_ind, col_ind] <- target_df[target_ind, ]$k_sim_1
        }
        
      }
    }
    
    df_heatmap_2 <- df_heatmap
    df_heatmap_2[is.na(df_heatmap_2)] <- 0
    
    non_zero_flag <- !(colSums(df_heatmap_2) == 1)
    df_heatmap_3 <- df_heatmap_2[non_zero_flag, non_zero_flag]
    
    annot <- data.frame(Cell_type = ctc_list_sorted[non_zero_flag])
    
    if (sum(non_zero_flag) == 0){
      
      return(NULL)
      
    } else {
      
      annot_list <- ctc_list_sorted[non_zero_flag]
      ant_list_sorted_2 <- ant_list_sorted_1[non_zero_flag]
      uniq_ctc_list <- unique(annot_list)
      uniq_ctc_num <- length(uniq_ctc_list)
      adjm <- matrix(0, nrow = uniq_ctc_num * 2, ncol = uniq_ctc_num * 2)
      rownames(adjm) <- c(paste(target_ant_1, " in ", uniq_ctc_list), paste(target_ant_2, " in ", uniq_ctc_list))
      colnames(adjm) <- c(paste(target_ant_1, " in ", uniq_ctc_list), paste(target_ant_2, " in ", uniq_ctc_list))
      
      ant_1_flag <- ant_list_sorted_2 == target_ant_1
      ant_2_flag <- ant_list_sorted_2 == target_ant_2
      
      # ant_1-ant_1 FOS-FOS
      for (uniq_ctc_ind_1 in 1:uniq_ctc_num){
        target_ctc_1 <- uniq_ctc_list[uniq_ctc_ind_1]
        target_row <- seq(1:length(annot_list))[target_ctc_1 == annot_list & ant_1_flag]
        for (uniq_ctc_ind_2 in 1:uniq_ctc_num){
          target_ctc_2 <- uniq_ctc_list[uniq_ctc_ind_2]
          target_col <- seq(1:length(annot_list))[target_ctc_2 == annot_list & ant_1_flag]
          adjm[uniq_ctc_ind_1, uniq_ctc_ind_2] <- mean(df_heatmap_3[target_row, target_col])
        }
      }
      
      # ant_2-ant_2 JUN-JUN
      for (uniq_ctc_ind_1 in 1:uniq_ctc_num){
        target_ctc_1 <- uniq_ctc_list[uniq_ctc_ind_1]
        target_row <- seq(1:length(annot_list))[target_ctc_1 == annot_list & ant_2_flag]
        for (uniq_ctc_ind_2 in 1:uniq_ctc_num){
          target_ctc_2 <- uniq_ctc_list[uniq_ctc_ind_2]
          target_col <- seq(1:length(annot_list))[target_ctc_2 == annot_list & ant_2_flag]
          adjm[uniq_ctc_ind_1 + uniq_ctc_num, uniq_ctc_ind_2 + uniq_ctc_num] <- mean(df_heatmap_3[target_row, target_col])
        }
      }
      
      # ant_1-ant_2 FOS-JUN
      for (uniq_ctc_ind_1 in 1:uniq_ctc_num){
        target_ctc_1 <- uniq_ctc_list[uniq_ctc_ind_1]
        target_row <- seq(1:length(annot_list))[target_ctc_1 == annot_list & ant_1_flag]
        for (uniq_ctc_ind_2 in 1:uniq_ctc_num){
          target_ctc_2 <- uniq_ctc_list[uniq_ctc_ind_2]
          target_col <- seq(1:length(annot_list))[target_ctc_2 == annot_list & ant_2_flag]
          adjm[uniq_ctc_ind_1, uniq_ctc_ind_2 + uniq_ctc_num] <- mean(df_heatmap_3[target_row, target_col])
        }
      }
      
      # ant_2-ant_1 JUN-FOS
      for (uniq_ctc_ind_1 in 1:uniq_ctc_num){
        target_ctc_1 <- uniq_ctc_list[uniq_ctc_ind_1]
        target_row <- seq(1:length(annot_list))[target_ctc_1 == annot_list & ant_2_flag]
        for (uniq_ctc_ind_2 in 1:uniq_ctc_num){
          target_ctc_2 <- uniq_ctc_list[uniq_ctc_ind_2]
          target_col <- seq(1:length(annot_list))[target_ctc_2 == annot_list & ant_1_flag]
          adjm[uniq_ctc_ind_1 + uniq_ctc_num, uniq_ctc_ind_2] <- mean(df_heatmap_3[target_row, target_col])
        }
      }
      
      adjm[is.na(adjm)] <- 0
      
      g <- igraph::graph_from_adjacency_matrix(adjmatrix = adjm, weighted = TRUE)
      igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, "red", "blue")
      igraph::E(g)$width <- abs(igraph::E(g)$weight) * 3
      igraph::V(g)$color <- "white"
      igraph::V(g)$x <- c(rep(1, uniq_ctc_num), rep(2, uniq_ctc_num))
      igraph::V(g)$y <- c(seq(1:uniq_ctc_num), seq(1:uniq_ctc_num))
      svg("~/MOCCS-DB_paper/plot/Fig3/Fig3X/FOS_JUN.pdf", width = 12)
      plot(g, edge.arrow.size = 0, vertex.label.color = "black")
      dev.off()
      
      if (nrow(df_heatmap_3) == 0){
        
        return(NULL)
        
      } else {
        
        mtxt <- paste0(target_ant_1, "-", target_ant_2, " k-sim 1 ")
        
        all_same_flag <- length(df_heatmap_3) == sum(df_heatmap_3 == df_heatmap_3[1,1])
        
        if (!all_same_flag){
          
          png(paste0("~/MOCCS-DB_paper/plot/Fig3/Fig3X/Fig3X_ALL/heatmap_2/png/heatmap_k_sim_1_",
                     target_ant_1, "_", target_ant_2, ".png"),
              width = 800, height = 720)
          NMF::aheatmap(df_heatmap_3, Rowv = NA, Colv = NA,
                        labRow = NA, labCol = NA,
                        annCol = annot, annRow = annot,
                        annColors = "Pastel1",
                        main = mtxt,
                        color = "Reds")
          dev.off()
          
          svg(paste0("~/MOCCS-DB_paper/plot/Fig3/Fig3X/Fig3X_ALL/heatmap_2/svg/heatmap_k_sim_1_",
                     target_ant_1, "_", target_ant_2, ".svg"),
              width = 800, height = 720)
          NMF::aheatmap(df_heatmap_3, Rowv = NA, Colv = NA,
                        labRow = NA, labCol = NA,
                        annCol = annot, annRow = annot,
                        annColors = "Pastel1",
                        main = mtxt,
                        color = "Reds")
          dev.off()
          
          
        } else {
          
          png(paste0("~/MOCCS-DB_paper/plot/Fig3/Fig3X/Fig3D_ALL/heatmap_2/png/heatmap_k_sim_1_",
                     target_ant_1, "_", target_ant_2, ".png"),
              width = 800, height = 720)
          NMF::aheatmap(df_heatmap_3, Rowv = NA, Colv = NA,
                        labRow = NA, labCol = NA,
                        annCol = annot, annRow = annot,
                        annColors = "Pastel1",
                        breaks = 0.5,
                        main = mtxt,
                        color = "Reds")
          dev.off()
          
          svg(paste0("~/MOCCS-DB_paper/plot/Fig3/Fig3X/Fig3D_ALL/heatmap_2/svg/heatmap_k_sim_1_",
                     target_ant_1, "_", target_ant_2, ".svg"),
              width = 800, height = 720)
          NMF::aheatmap(df_heatmap_3, Rowv = NA, Colv = NA,
                        labRow = NA, labCol = NA,
                        annCol = annot, annRow = annot,
                        annColors = "Pastel1",
                        breaks = 0.5,
                        main = mtxt,
                        color = "Reds")
          dev.off()
          
        }
        
      }  
      
    }
    
  }
  
}