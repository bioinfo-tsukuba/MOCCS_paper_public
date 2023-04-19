TF_graph <- function(df_p_3){
  
  library(ggraph)
  library(tidyverse)
  library(igraph)
  
  df_p_3[is.na(df_p_3)] <- 0
  
  receiver_tf_list <- c("FOS", "JUN",
                        "ATF3", "E2F1", "EGR1",
                        "FOXA1", "FOXA2", "GABPA",
                        "JUND", "MAX", "NANOG", "REST", "TAF1")
  
  sender_tf_list <- c("FOS", "JUN",
                      "ATF3", "E2F1", "EGR1",
                      "FOXA1", "FOXA2", "GABPA",
                      "JUND", "MAX", "NANOG", "REST", "TAF1")
  
  for (i in 1:length(receiver_tf_list)){
    
    at_1 <- receiver_tf_list[i]
    system(paste0("mkdir ~/MOCCS-DB_paper/plot/Fig2/Fig2E/Fig2E_", at_1, "/"))
    target_flag_1 <- df_p_3$ID1_Antigen == at_1 | df_p_3$ID2_Antigen == at_1
    df_at_1 <- df_p_3[target_flag_1, ]
    ctc_list <- union(df_at_1$ID1_Cell_type_class, df_at_1$ID2_Cell_type_class)

    for (ctc_ind in 1:length(ctc_list)){
      
      ctc_1 <- ctc_list[ctc_ind]
      adjm <- matrix(0, nrow = length(sender_tf_list) + 1, ncol = length(sender_tf_list) + 1)
      del_list <- c()
      
      for (j in 1:length(sender_tf_list)){
        
        at_2 <- sender_tf_list[j]
        target_flag_2 <- df_p_3$ID1_Antigen == at_1 & df_p_3$ID2_Antigen == at_2 & df_p_3$ID1_Cell_type_class == ctc_1 & df_p_3$ID2_Cell_type_class == ctc_1
        target_flag_3 <- df_p_3$ID1_Antigen == at_2 & df_p_3$ID2_Antigen == at_1 & df_p_3$ID1_Cell_type_class == ctc_1 & df_p_3$ID2_Cell_type_class == ctc_1
        target_flag_4 <- target_flag_2 | target_flag_3
        df_pair <- df_p_3[target_flag_4, ]
        
        if (nrow(df_pair) > 0){
          
          adjm[j, length(sender_tf_list) + 1] <- mean(df_pair$k_sim_2)
          # print(mean(df_pair$k_sim_2))
          
        } else {
          
          adjm[j, length(sender_tf_list) + 1] <- 0
          del_list <- append(del_list, j)
          
        }
      }
      
      if (!is.null(del_list)){
        adjm_2 <- adjm[- del_list, - del_list]
        sender_tf_list_2 <- sender_tf_list[- del_list]
      } else {
        adjm_2 <- adjm
        sender_tf_list_2 <- sender_tf_list
      }
      
      if (sum(adjm_2) != 0){
        
        g <- graph_from_adjacency_matrix(adjmatrix = adjm_2, weighted = TRUE)
        V(g)$name <- c(sender_tf_list_2, receiver_tf_list[i])
        xy <- graphlayouts::layout_with_stress(g)
        xy <- graphlayouts::layout_rotate(xy, 45)
        p_gg <- ggraph(g, "manual", x = xy[,1], y = xy[,2]) +
          geom_node_point(size = 12, colour = "grey") +
          geom_node_text(aes(label = name)) +
          geom_edge_link(aes(color = weight),
                         start_cap = circle(7, 'mm'),
                         end_cap = circle(7, 'mm'), width = 1) +
          scale_edge_colour_gradient(low = "white", high = "red") +
          theme_graph() +
          theme(legend.position = "none") +
          coord_equal()
          # ggtitle(paste0(at_1, " in ", ctc_1))
        ggsave(paste0("~/MOCCS-DB_paper/plot/Fig2/Fig2E/Fig2E_", at_1, "/Fig2E_graph_", at_1, "_", ctc_1, ".png"), plot = p_gg, width = 5, height = 5)
        ggsave(paste0("~/MOCCS-DB_paper/plot/Fig2/Fig2E/Fig2E_", at_1, "/Fig2E_graph_", at_1, "_", ctc_1, ".pdf"), plot = p_gg, width = 5, height = 5)  
        
      }
      
    }
    
  }
  
  
  return()
  
}

TF_graph_all_ctc <- function(df_p_3){ # Not divide by cell type class
  
  library(ggraph)
  library(tidyverse)
  library(igraph)
  
  df_p_3[is.na(df_p_3)] <- 0
  
  receiver_tf_list <- c("ATF3", "CTCF", "E2F1", "EGR1",
                      "FOXA1", "FOXA2", "GABPA",
                      "JUND", "MAX", "NANOG", "REST", "TAF1",
                      "FOS", "JUN")
  
  sender_tf_list <- c("ATF3", "CTCF", "E2F1", "EGR1",
                      "FOXA1", "FOXA2", "GABPA",
                      "JUND", "MAX", "NANOG", "REST", "TAF1",
                      "FOS", "JUN")
  
  for (i in 1:length(receiver_tf_list)){
    
    at_1 <- receiver_tf_list[i]
    system(paste0("mkdir ~/MOCCS-DB_paper/plot/Fig2/Fig2E_", at_1, "/"))
    target_flag_1 <- df_p_3$ID1_Antigen == at_1 | df_p_3$ID2_Antigen == at_1
    df_at_1 <- df_p_3[target_flag_1, ]
    ctc_list <- union(df_at_1$ID1_Cell_type_class, df_at_1$ID2_Cell_type_class)
    adjm_list <- vector("list", length(ctc_list))
    
    for (ctc_ind in 1:length(ctc_list)){
      adjm <- matrix(0, nrow = length(sender_tf_list) + 1, ncol = length(sender_tf_list) + 1)
      adjm_list[[ctc_ind]] <- adjm
    }
    
    del_list <- c()
    
    for (j in 1:length(sender_tf_list)){
      
      at_2 <- sender_tf_list[j]
      target_flag_2 <- df_p_3$ID1_Antigen == at_1 & df_p_3$ID2_Antigen == at_2
      target_flag_3 <- df_p_3$ID1_Antigen == at_2 & df_p_3$ID2_Antigen == at_1
      target_flag_4 <- target_flag_2 | target_flag_3
      df_pair <- df_p_3[target_flag_4, ]
      
      if (nrow(df_pair) > 0){
        
        adjm[j, length(sender_tf_list) + 1] <- mean(df_pair$k_sim_2)
        print(mean(df_pair$k_sim_2))
        
      } else {
        
        adjm[j, length(sender_tf_list) + 1] <- 0
        del_list <- append(del_list, j)
        
      }
    }
    
    if (!is.null(del_list)){
      adjm_2 <- adjm[- del_list, - del_list]
    } else {
      adjm_2 <- adjm
    }
    
    g <- graph_from_adjacency_matrix(adjmatrix = adjm_2, weighted = TRUE)
    V(g)$name <- c(sender_tf_list, receiver_tf_list[i])
    xy <- graphlayouts::layout_with_stress(g)
    xy <- graphlayouts::layout_rotate(xy, 45)
    p_gg <- ggraph(g, "manual", x = xy[,1], y = xy[,2]) +
      geom_node_point(size = 12, colour = "grey") +
      geom_node_text(aes(label = name)) +
      geom_edge_link(aes(color = weight), arrow = arrow(length = unit(3, 'mm')),
                     start_cap = circle(5, 'mm'),
                     end_cap = circle(7, 'mm'), width = 1) +
      scale_edge_colour_gradient(low = "white", high = "red") +
      theme_graph() +
      theme(legend.position = "none") +
      # ggtitle(paste0(at_1, " in all the cell types")) +
      + coord_equal()
    ggsave(paste0("~/MOCCS-DB_paper/plot/Fig2/Fig2E/Fig2E_", at_1, "/Fig2E_graph_", at_1, "_all.png"), plot = p_gg, width = 5, height = 5)
    ggsave(paste0("~/MOCCS-DB_paper/plot/Fig2/Fig2E/Fig2E_", at_1, "/Fig2E_graph_", at_1, "_all.pdf"), plot = p_gg, width = 5, height = 5)  
    
  }
  
  
  return()
  
}