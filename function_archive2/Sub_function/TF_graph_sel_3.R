TF_graph_sel_3 <- function(a1, t_ctc_list, t_num){
  
  # a1 <- "FOS"
  # t_ctc_list <- c("Blood", "Breast")
  # t_num <- 15
  
  library(tidyverse)
  library(ggraph)
  library(igraph)
  
  a1_mat <- readRDS(paste0("~/MOCCS_paper_public/data/Fig3/obj/k_sim_jaccard_mat/k_sim_jaccard_mat_", a1, ".rds"))
  
  if (is.null(t_ctc_list)){
    s_num_list <- c()
    for (row_i in 1:nrow(a1_mat)){
      s_num <- sum(!is.na(a1_mat[row_i, ]))
      s_num_list <- append(s_num_list, s_num)
    }
    t_ctc_list <- rownames(a1_mat)[order(s_num_list, decreasing = TRUE)][1:2]
  }
  
  t_a1_1 <- a1_mat[t_ctc_list[1],][order(a1_mat[t_ctc_list[1],], decreasing = TRUE)]
  t_a1_2 <- a1_mat[t_ctc_list[2],][order(a1_mat[t_ctc_list[2],], decreasing = TRUE)]
  
  f_uni <- intersect(names(t_a1_1[!is.na(t_a1_1)]), names(t_a1_2[!is.na(t_a1_2)]))

  f_1 <- t_a1_1[names(t_a1_1) %in% f_uni]
  f_2 <- t_a1_2[names(t_a1_2) %in% f_uni]
  
  f_12 <- union(names(f_1), names(f_2))

  f_a1_1 <- t_a1_1[names(t_a1_1) %in% f_12]
  f_a1_2 <- t_a1_2[names(t_a1_2) %in% f_12]
  
  diff_f_a1 <- f_a1_1 - f_a1_2
  diff_f_2_a1 <- diff_f_a1[order(abs(diff_f_a1), decreasing = TRUE)][1:t_num]
  f_2_12 <- names(diff_f_2_a1)
  
  f_2_a1_1 <- t_a1_1[names(t_a1_1) %in% f_2_12]
  f_2_a1_2 <- t_a1_2[names(t_a1_2) %in% f_2_12]
  
  f_a1_list <- vector("list", length(t_ctc_list))
  f_a1_list[[1]] <- f_2_a1_1[order(f_2_a1_1)]
  f_a1_list[[2]] <- f_2_a1_2[names(f_a1_list[[1]])]
  f_3_12 <- names(f_a1_list[[1]])
  
  min_val <- min(c(f_a1_list[[1]], f_a1_list[[2]]))
  min_val_dc <- round(min_val * 100)

  max_val <- max(c(f_a1_list[[1]], f_a1_list[[2]]))
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
  #pdf (paste0("~/MOCCS_paper_public/plot/Fig3/Fig3F/Fig3F_", a1, "_color_legend.pdf"))
  pdf (paste0("~/MOCCS_paper_public/plot/Fig2/Fig2F/Fig2F_", a1, "_color_legend.pdf"))
  RColorBrewer::display.brewer.pal(9, "Reds")
  dev.off()

  for (t_ctc_ind in 1:length(t_ctc_list)){
    
    t_ctc <- t_ctc_list[t_ctc_ind]
    
    rc_num <- length(f_3_12) + 1
    adjm <- matrix(0, nrow = rc_num, ncol = rc_num)
    colnames(adjm) <- c(f_3_12, a1)
    rownames(adjm) <- c(f_3_12, a1)
      
    for (k in 1:rc_num - 1){
      adjm[k, rc_num] <- f_a1_list[[t_ctc_ind]][k]
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
      scale_edge_colour_gradient(low = col_min, high = col_max) +
      theme_graph() +
      theme(legend.position = "none") +
      coord_equal()
      # labs(title = paste0(a1 , " in ", t_ctc))
    #ggsave(paste0("~/MOCCS_paper_public/plot/Fig3/Fig3F/Fig3F_", a1, "_", t_ctc, ".pdf"), plot = p_gg, width = 5, height = 5)  
    ggsave(paste0("~/MOCCS_paper_public/plot/Fig2/Fig2F/Fig2F_", a1, "_", t_ctc, ".pdf"), plot = p_gg, width = 5, height = 5)  

  }
    
}