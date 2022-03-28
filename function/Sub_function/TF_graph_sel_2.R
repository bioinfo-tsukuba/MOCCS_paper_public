TF_graph_sel_2 <- function(rt_list, t_ctc_list, t_num = 10){
  
  # rt_list <- c("FOS", "JUN")
  # t_ctc_list <- c("Blood", "Breast")
  # t_num <- 3
  
  library(ggraph)
  library(igraph)
  
  st_cand_list_1 <- vector("list", length(rt_list))
  t_row_non_na_2_list <- vector("list", length(t_ctc_list))
  
  for (t_ctc_ind in 1:length(t_ctc_list)){
    
    t_ctc <- t_ctc_list[t_ctc_ind]
    
    for (i in 1:length(rt_list)){
      
      at_1 <- rt_list[i]
      print(at_1)
      k_sim_1_mat <- readRDS(paste0("~/MOCCS-DB_paper/results/k_sim_1_mat/k_sim_1_mat_", at_1, ".rds"))
      
      t_row <- k_sim_1_mat[rownames(k_sim_1_mat) == t_ctc]
      names(t_row) <- colnames(k_sim_1_mat)
      
      t_row_non_na <- t_row[!is.na(t_row)]
      t_row_non_na_2 <- t_row_non_na[t_row_non_na != 0] 
      
      t_row_non_na_2_list[[t_ctc_ind]][[i]] <- t_row_non_na_2
      st_cand_list_1[[t_ctc_ind]][[i]] <- names(t_row_non_na_2[order(t_row_non_na_2, decreasing = TRUE)][1:t_num])
      print(paste0(st_cand_list_1[[t_ctc_ind]][[i]], " - ", t_ctc))
      
    }
    
  }

  st_cand_list_uni <- c()  
  for (t_ctc_ind in 1:length(t_ctc_list)){
    for (i in 1:length(rt_list)){
      print(st_cand_list_1[[t_ctc_ind]][[i]])
      st_cand_list_uni <- union(st_cand_list_uni, st_cand_list_1[[t_ctc_ind]][[i]])
    }
  }
  
  print(st_cand_list_uni)
  
  for (t_ctc_ind in 1:length(t_ctc_list)){
    for (i in 1:length(rt_list)){
      for (i_st in 1:length(st_cand_list_uni)){
        print(st_cand_list_uni[i_st] %in% names(t_row_non_na_2_list[[t_ctc_ind]][[i]]))
      }
    }
  }
  
  sender_tf_list <- st_cand_list_uni
  
  for (t_ctc_ind in 1:length(t_ctc_list)){
    
    t_ctc <- t_ctc_list[t_ctc_ind]
    
    for (i in 1:length(rt_list)){
      
      at_1 <- rt_list[i]
      system(paste0("mkdir ~/MOCCS-DB_paper/plot/Fig3/Fig3F/Fig3F_", at_1, "/"))
      
      k_sim_2_mat <- readRDS(paste0("~/MOCCS-DB_paper/results/k_sim_2_mat/k_sim_2_mat_", at_1, ".rds"))
      
      t_row <- k_sim_2_mat[rownames(k_sim_2_mat) == t_ctc]
      names(t_row) <- colnames(k_sim_2_mat)
      
      t_row_non_na <- t_row[!is.na(t_row)]
      
      t_row_non_na_top <- rep(0, length(sender_tf_list))
      names(t_row_non_na_top) <- sender_tf_list
      for (j in 1:length(sender_tf_list)){
        sbs_j <- t_row_non_na[names(t_row_non_na) == sender_tf_list[j]]
        if (length(sbs_j) > 0){
          t_row_non_na_top[j] <- sbs_j
        }
      }
      
      row_col_num <- length(sender_tf_list) + 1
      adjm <- matrix(0, nrow = row_col_num, ncol = row_col_num)
      colnames(adjm) <- c(names(t_row_non_na_top), at_1)
      rownames(adjm) <- c(names(t_row_non_na_top), at_1)
      
      for (k in 1:row_col_num - 1){
        adjm[k, row_col_num] <- t_row_non_na_top[k]
      }
      
      g <- graph_from_adjacency_matrix(adjmatrix = adjm, weighted = TRUE)
      V(g)$name <- colnames(adjm)
      g_dm <- graph_from_adjacency_matrix(adjmatrix = adjm[-1,-1])
      xy <- igraph::layout_in_circle(g_dm)
      xy <- rbind(xy, c(0, 0))
      p_gg <- ggraph(g, layout = "manual", x = xy[,1], y = xy[,2]) +
        geom_node_point(size = 12, colour = "grey") +
        geom_node_text(aes(label = name)) +
        geom_edge_link(aes(color = weight),
                       start_cap = circle(7, 'mm'),
                       end_cap = circle(7, 'mm'), width = 1) +
        scale_edge_colour_gradient(low = "white", high = "red") +
        theme_graph() +
        theme(legend.position = "none") +
        coord_equal()
      ggsave(paste0("~/MOCCS-DB_paper/plot/Fig3/Fig3F/Fig3F_", at_1, "/Fig3F_", at_1, "_", t_ctc, ".png"), plot = p_gg, width = 5, height = 5)
      ggsave(paste0("~/MOCCS-DB_paper/plot/Fig3/Fig3F/Fig3F_", at_1, "/Fig3F_", at_1, "_", t_ctc, ".pdf"), plot = p_gg, width = 5, height = 5)  
      
    }
    
  }
  
}