Heatmap_db <- function(df_p_3_gp, target_ant){
  
  sig_symbol_list <- readRDS("~/MOCCS_paper_public/data/Fig3/obj/sig_symbol_list.rds")
  if (target_ant %in% sig_symbol_list){
    sig_flag <- "*"
  } else{
    sig_flag <- ""
  }
  
  target_flag <- df_p_3_gp$ID1_Antigen == target_ant & df_p_3_gp$ID2_Antigen == target_ant & df_p_3_gp$s_ctc != 0.5
  target_df <- df_p_3_gp[target_flag, ]
  
  if (nrow(target_df) == 0){
    
    print("No data exist...")
    return(NULL)
    
  } else {
    
    id_list <- union(unique(target_df$ID1), unique(target_df$ID2))
    ctc_list <- c()
    
    for (id_ind in 1:length(id_list)){
      target_id <- id_list[id_ind]
      target_ind <- target_df$ID1 == target_id
      target_ctc <- unique(target_df$ID1_Cell_type_class[target_ind])
      if (length(target_ctc) == 0){
        target_ind <- target_df$ID2 == target_id
        target_ctc <- unique(target_df$ID2_Cell_type_class[target_ind])
      }
      ctc_list <- append(ctc_list, target_ctc)
    }
    
    order_ind <- order(ctc_list)
    ctc_list_sorted <- ctc_list[order_ind]
    id_list_sorted <- id_list[order_ind]
    
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
    
    annot <- data.frame(Cell_type_class = ctc_list_sorted)

    if (nrow(df_heatmap_2) == 0){
      
      return(NULL)
      
    } else {
      
      if (sig_flag == "*"){
        mtxt <- paste0(target_ant, " ", sig_flag)
      } else {
        mtxt <- paste0(target_ant, "")
      }
      
      col_8 <- RColorBrewer::brewer.pal(8, "Accent")
      col_20 <- colorRampPalette(col_8)(20)
      col_id <- c(seq(from = 1, to = 16, by = 5),
                  seq(from = 2, to = 17, by = 5),
                  seq(from = 3, to = 18, by = 5),
                  seq(from = 4, to = 19, by = 5),
                  seq(from = 5, to = 20, by = 5))
      col_20_2 <- col_20[col_id]
      
      png(paste0("~/MOCCS_paper_public/plot/FigS6/heatmap_k_sim_jaccard_", target_ant, ".png"),
          width = 900, height = 720)
      NMF::aheatmap(df_heatmap_2, Rowv = NA, Colv = NA,
                    labRow = NA, labCol = NA,
                    annCol = annot, annRow = annot,
                    annColors = list(Cell_type_class = col_20_2),
                    main = mtxt,
                    color = "Reds",
                    fontsize = 20)
      dev.off()
        
    } 
    
  }
  
}