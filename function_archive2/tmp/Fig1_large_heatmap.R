Large_Heatmap <- function(df_p_3_gp, target_ant, cal_opt){
  
  
  library(colorspace)
  if(cal_opt == TRUE){
    count <- 0
    for (i in 1:length(target_ant)) {
      print(paste0(i, "/", length(target_ant)))
      tgt_TF <- target_ant[i]
      target_flag <- df_p_3_gp$ID1_Antigen == tgt_TF & df_p_3_gp$ID2_Antigen == tgt_TF #& df_p_3_gp$s_ctc != 0.5
      target_df <- df_p_3_gp[target_flag, ]
      count <- count + nrow(target_df)
    }
    
    large_mat <- data.frame(matrix(nrow = count, ncol = ncol(df_p_3_gp)))
    num <- 0
    for (i in 1:length(target_ant)) {
      print(paste0(i, "/", length(target_ant)))
      tgt_TF <- target_ant[i]
      target_flag <- df_p_3_gp$ID1_Antigen == tgt_TF & df_p_3_gp$ID2_Antigen == tgt_TF #& df_p_3_gp$s_ctc != 0.5
      target_df <- df_p_3_gp[target_flag, ]
      
      num_start <- num+1
      num_end <- num+nrow(target_df)
      large_mat[num_start:num_end, ] <- target_df
      num <- num + nrow(target_df)
      
    }
    colnames(large_mat) <- colnames(target_df)
    large_mat <- as_tibble(large_mat)
    large_mat <- large_mat %>% arrange(ID1_Family, ID2_Family, ID1_Antigen, ID2_Antigen)
    saveRDS(large_mat, "/Users/saeko/MOCCS_paper_public/function/tmp/large_mat.rds")
  }
  
  
  # heatmap matrix----
  #large_mat <- df_p_3_gp %>% arrange(ID1_Family, ID2_Family, ID1_Antigen, ID2_Antigen)
  large_mat <- readRDS("/Users/saeko/MOCCS_paper_public/function/tmp/large_mat.rds")
  id_list <- union(unique(large_mat$ID1), unique(large_mat$ID2))
  df_heatmap <- matrix(0, nrow = length(id_list), ncol = length(id_list))
  
  for (row_ind in 1:length(id_list)){
    print(paste0(row_ind, "/", length(id_list)))
    row_id <- id_list[row_ind]
    for (col_ind in 1:length(id_list)){
      #print(col_ind)
      col_id <- id_list[col_ind]
      if (row_id == col_id){
        df_heatmap[row_ind, col_ind] <- 1 # same id
      } else {
        tmp <- large_mat %>% filter((ID1 == row_id & ID2 == col_id)|(ID1 == col_id & ID2 == row_id))
        if(nrow(tmp)==0){
          df_heatmap[row_ind, col_ind] <- NA
        }else{
          #large_mat %>% filter((ID1 == row_id & ID2 == col_id)|(ID1 == col_id & ID2 == row_id))
          target_ind <- (large_mat$ID1 == row_id & large_mat$ID2 == col_id) | (large_mat$ID1 == col_id & large_mat$ID2 == row_id)
          df_heatmap[row_ind, col_ind] <- large_mat[target_ind, ]$k_sim_1
        }
      }
    }
  }
  
  saveRDS(df_heatmap, "~/MOCCS_paper_public/plot/tmp/df_heatmap_all.rds")
  df_heatmap_2 <- df_heatmap
  df_heatmap_2[is.na(df_heatmap_2)] <- 0
  
  # id_listに対応するAntigen listを作成
  ant_list <- c()
  for (i in 1:length(id_list)) {
    tgt_id <- id_list[i]
    tgt_row <- large_mat %>% filter(ID1 == tgt_ant | ID2 == tgt_ant) 
    tgt_ant <- unique(tgt_row$ctc1, tgt_row$ctc2)
    ant_list <- c(ant_list, tgt_ant)
  }
  annot <- data.frame(Antigen = ant_list)
  
  
  if (nrow(df_heatmap_2) == 0){
    return(NULL)
  } else {
    mtxt <- paste0(target_ant, "")
    col_20_2 <- qualitative_hcl(length(unique(annot$Antigen)), c = 50, l = 70)
    
    png(paste0("~/MOCCS_paper_public/plot/tmp/heatmap_all_jaccard.png"),
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