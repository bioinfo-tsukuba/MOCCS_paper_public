TF_graph_sel_db <- function(ctc_spe, df_p_3, receiver_tf_list, t_num = 10){
  
  library(ggraph)
  library(igraph)
  
  for (i in 1:length(receiver_tf_list)){
  
    at_1 <- receiver_tf_list[i]

    k_sim_2_mat <- readRDS(paste0("~/MOCCS_paper_public/data/Fig2/obj/k_sim_pearson_mat/k_sim_pearson_mat_", at_1, ".rds"))

    if (is.null(ctc_spe)){
      tf_non_na_num_list <- c()
      for (row_i in 1:nrow(k_sim_2_mat)){
        tf_non_na <- k_sim_2_mat[row_i, ][!is.na(k_sim_2_mat[row_i, ])]
        tf_non_na_num <- length(tf_non_na)
        tf_non_na_num_list <- append(tf_non_na_num_list, tf_non_na_num)
      }
      t_ctc <- rownames(k_sim_2_mat)[which.max(tf_non_na_num_list)]
    }
    
    t_row <- k_sim_2_mat[rownames(k_sim_2_mat) == t_ctc]
    names(t_row) <- colnames(k_sim_2_mat)
    
    t_row_non_na <- t_row[!is.na(t_row)]
    t_row_non_na_top <- t_row_non_na[order(t_row_non_na, decreasing = TRUE)][1:t_num]
    
    adjm <- matrix(0, nrow = t_num + 1, ncol = t_num + 1)
    colnames(adjm) <- c(names(t_row_non_na_top), at_1)
    rownames(adjm) <- c(names(t_row_non_na_top), at_1)
    
    for (j in 1:t_num){
      adjm[j, t_num + 1] <- t_row_non_na_top[j]
    }

    adjm_self <- t_row_non_na[names(t_row_non_na) == at_1]
    if (length(adjm_self) != 0){
      adjm[nrow(adjm), ncol(adjm)] <- t_row_non_na[names(t_row_non_na) == at_1]
    }
    
    min_val <- min(adjm[, t_num + 1])
    min_val_dc <- round(min_val * 100)
    
    max_val <- max(adjm[, t_num + 1])
    max_val_dc <- round(max_val * 100)
    
    col_gg <- colorRampPalette(RColorBrewer::brewer.pal(100, "Reds"))
    
    print(paste0("minval: ", min_val))
    print(paste0("maxval: ", max_val))
    
    if (min_val != 0){
      col_min <- col_gg(100)[min_val_dc]
    } else {
      col_min <- "white"
    }
    
    if (max_val != 1){
      col_max <- col_gg(100)[max_val_dc]
    } else {
      col_max <- "red"
    }
    pdf (paste0("~/MOCCS_paper_public/plot/Fig2/Fig2E/Fig2E_", at_1, "_color_legend.pdf"))
    RColorBrewer::display.brewer.pal(9, "Reds")
    dev.off()
    
    print(paste0("=== ", at_1, " ",
                 unique(df_p_3[df_p_3$ID1_Antigen == at_1,]$ID1_Family)," ",
                 t_ctc, " ==="))
    t_l <- colnames(adjm)
    for (t_i in 1:length(t_l)){
      t_s <- t_l[t_i]
      t_s_f <- unique(df_p_3[df_p_3$ID1_Antigen == t_s,]$ID1_Family)
      print(paste0(t_s, " - ", t_s_f))
    }
    print(paste0("==="))
    
    g <- graph_from_adjacency_matrix(adjmatrix = adjm, weighted = TRUE)
    V(g)$name <- colnames(adjm)
    xy <- graphlayouts::layout_with_stress(g)
    xy <- graphlayouts::layout_rotate(xy, 45)
    p_gg <- ggraph(g, "manual", x = xy[,1], y = xy[,2]) +
      geom_node_point(size = 12, colour = "grey") +
      geom_node_text(aes(label = name)) +
      geom_edge_link(aes(color = weight),
                     start_cap = circle(7, 'mm'),
                     end_cap = circle(7, 'mm'), width = 1) +
      scale_edge_colour_gradient(low = col_min, high = col_max) +
      theme_graph() +
      theme(legend.position = "none") +
      coord_equal()
    # ggsave(paste0("~/MOCCS_paper_public/plot/Fig2/Fig2E/Fig2E_top_10_", at_1, "_", t_ctc, ".png"), plot = p_gg, width = 5, height = 5)
    ggsave(paste0("~/MOCCS_paper_public/plot/Fig2/Fig2E/Fig2E_top_10_", at_1, "_", t_ctc, ".pdf"), plot = p_gg, width = 5, height = 5)  
    
  }
  
}