Compare_cell_type <- function(df_p_3_gp, plot_ant){
  
  library(patchwork)
  library(ComplexHeatmap)
  library(colorspace)
  library(dplyr)
  library(ggplot2)
  
  ant_list <- union(unique(df_p_3_gp$ID1_Antigen), unique(df_p_3_gp$ID2_Antigen))
  
  df_p_3_gp[is.na(df_p_3_gp)] <- 0
  
  # ant_list <- c("FOS", "JUN", "JUND", "GATA2")
  
  sig_flag_list <- c()
  p_val_list <- c()
  for (ant_ind in 1:length(ant_list)){
    res.viol <- Violinplot(df_p_3_gp, target_ant = ant_list[ant_ind], plot_ant)
    Heatmap(df_p_3_gp, target_ant = ant_list[ant_ind], res.viol$sig_flag, plot_ant)
    sig_flag_list <- append(sig_flag_list, res.viol$sig_flag)
    p_val_list <- append(p_val_list, res.viol$p_val)
  }
  
  saveRDS(sig_flag_list, "~/MOCCS_paper_public/data/Fig3/obj/sig_flag_list.rds")
  saveRDS(p_val_list, "~/MOCCS_paper_public/data/Fig3/obj/p_val_list.rds")
  
  all_Violin_plot(df_p_3_gp, ant_list, sig_flag_list)
  
  
  # Make table of cell-type dependent TFs
  Make_TF_table(ant_list)
  Make_TF_table2()
  
  # Fig. 2E
  CircleChart(sig_flag_list)
  
  
  return()
  
}

Violinplot <- function(df_p_3_gp, target_ant, plot_ant){
  
  target_flag <- df_p_3_gp$ID1_Antigen == target_ant & df_p_3_gp$ID2_Antigen == target_ant
  target_df <- df_p_3_gp[target_flag, ]
  
  
  ####### 20230112 Added ######
  tmp1 <- target_df %>% select(ID1, ID1_Cell_type_class)
  tmp2 <- target_df %>% select(ID2, ID2_Cell_type_class)
  colnames(tmp1) <- c("ID", "Cell_type_class")
  colnames(tmp2) <- c("ID", "Cell_type_class")
  tmp3 <- rbind(tmp1, tmp2) %>% distinct()
  rm_CTC <- tmp3 %>% group_by(Cell_type_class) %>% summarise(n = n()) %>% filter(n == 1) %>% .$Cell_type_class
  target_df_old <- target_df
  target_df <- target_df_old %>% filter(!ID1_Cell_type_class %in% rm_CTC & !ID2_Cell_type_class %in% rm_CTC) 
  #############################
  
  target_df[target_df$s_ctc == 1, ] %>% #s_ctc=1でCell type classが同じことを表す
    mutate(Group = "Same") -> df_same
  
  target_df[target_df$s_ctc == 0, ] %>% #s_ctc=0でCell type classが異なることを表す
    mutate(Group = "Different") -> df_diff
  
  if (nrow(df_same) == 0 || nrow(df_diff) == 0){
    
    return(list(sig_flag = "NULL",
                p_val = "NULL"))
    
  } else {
    
    df_viol <- rbind(df_same, df_diff)
    
    p_violin <- ggplot(df_viol, aes(x = Group, y = k_sim_1)) +
      geom_violin(scale="width",
                  adjust = 1,
                  width = 0.5,
                  fill = "gray80") +
      theme_bw() +
      geom_boxplot(width = .1, fill = "white", outlier.colour = NA) +
      stat_summary(fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
      ylab("k-sim Jaccard") +
      xlab("") +
      ylim(0, 1.05) +
      theme(panel.grid.minor = element_blank()) +
      ggtitle(paste0(target_ant))
    
    res.wil <- exactRankTests::wilcox.exact(x = df_same$k_sim_1, y = df_diff$k_sim_1, paired = FALSE)
    p_val <- res.wil$p.value
    if (p_val < 0.05){
      sig_flag <- "*"
      p_violin <- p_violin + 
        # geom_text(x = 1.5, y = 1.06, label = "*") +
        geom_segment(x = 1, xend = 2, y = 1.05, yend = 1.05) +
        geom_segment(x = 1, xend = 1, y = 1.02, yend = 1.05) +
        geom_segment(x = 2, xend = 2, y = 1.02, yend = 1.05)
    } else {
      sig_flag <- ""
    }
    
    if (target_ant %in% plot_ant){
      # ggsave(paste0("~/MOCCS_paper_public/plot/Fig3/Fig3D/Fig3D_violin_plot_k_sim_jaccard_",
      #              target_ant, ".png"),
      #       plot = p_violin, width = 3, height = 6)
      ggsave(paste0("~/MOCCS_paper_public/plot/Fig3/Fig3D/Fig3D_violin_plot_k_sim_jaccard_",
                    target_ant, ".svg"),
             plot = p_violin, width = 3, height = 6)
    }
    
    return(list(sig_flag = sig_flag,
                p_val = p_val,
                p_violin = p_violin))
    
  }
  
}

Violinplot_2 <- function(df_p_3_gp, target_ant){
  
  target_flag <- df_p_3_gp$ID1_Antigen == target_ant & df_p_3_gp$ID2_Antigen == target_ant
  target_df <- df_p_3_gp[target_flag, ]
  
  ####### 20230112 Added ######
  tmp1 <- target_df %>% select(ID1, ID1_Cell_type_class)
  tmp2 <- target_df %>% select(ID2, ID2_Cell_type_class)
  colnames(tmp1) <- c("ID", "Cell_type_class")
  colnames(tmp2) <- c("ID", "Cell_type_class")
  tmp3 <- rbind(tmp1, tmp2) %>% distinct()
  rm_CTC <- tmp3 %>% group_by(Cell_type_class) %>% summarise(n = n()) %>% filter(n == 1) %>% .$Cell_type_class
  target_df_old <- target_df
  target_df <- target_df_old %>% filter(!ID1_Cell_type_class %in% rm_CTC & !ID2_Cell_type_class %in% rm_CTC) 
  #############################
  
  target_df[target_df$s_ctc == 1, ] %>%
    mutate(Group = "Same") -> df_same
  
  target_df[target_df$s_ctc == 0, ] %>%
    mutate(Group = "Different") -> df_diff
  
  if (nrow(df_same) == 0 || nrow(df_diff) == 0){
    
    return(list(sig_flag = "NULL",
                p_val = "NULL"))
    
  } else {
    
    df_viol <- rbind(df_same, df_diff)
    
    p_violin <- ggplot(df_viol, aes(x = Group, y = k_sim_1)) +
      geom_violin(scale="width",
                  adjust = 1,
                  width = 0.5,
                  fill = "gray80") +
      theme_bw() +
      geom_boxplot(width = .1, fill = "white", outlier.colour = NA) +
      stat_summary(fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
      ylab("k-sim Jaccard") +
      xlab("") +
      ylim(0, 1) +
      theme(panel.grid.minor = element_blank()) +
      ggtitle(paste0(target_ant, " *"))
    
    
    return(list(p_violin = p_violin))
    
  }
  
}

all_Violin_plot <- function(df_p_3_gp, ant_list, sig_flag_list){
  
  
  ########## 2023/01/12 Added ############
  p_violin_list <- list()
  for (tgt_ant in ant_list){
    res.viol <- Violinplot_2(df_p_3_gp, target_ant = tgt_ant)
    if(is.null(res.viol$p_violin)==TRUE){
      p_violin_list[[tgt_ant]] <- "NULL"
    }else{
      p_violin_list[[tgt_ant]] <- res.viol$p_violin 
    }
  }
  
  sig_ant <- ant_list[sig_flag_list == "*"]
  sig_ant_2 <- sig_ant[order(sig_ant)]
  sig_vec <- c()
  for (sig_i in 1:length(sig_ant_2)){
    sig_vec <- append(sig_vec, seq(1:length(ant_list))[sig_ant_2[sig_i] == ant_list])
  }
  
  p_patch_1 <- c()
  p_patch_1 <- p_violin_list[[sig_vec[1]]]
  for (i in 2:20){
    i_2 <- sig_vec[i]
    p_patch_1 <- p_patch_1 + p_violin_list[[i_2]]
  }
  
  p_patch_2 <- c()
  p_patch_2 <- p_violin_list[[sig_vec[21]]]
  for (i in 22:length(sig_vec)){
    i_2 <- sig_vec[i]
    p_patch_2 <- p_patch_2 + p_violin_list[[i_2]]
  }
  ########################################
  
  pdf("~/MOCCS_paper_public/plot/FigS7/FigS7_viol_1.pdf", paper = "a4")
  plot(p_patch_1 + plot_layout(nrow = 5, ncol = 4), height = 50)
  dev.off()
  
  pdf("~/MOCCS_paper_public/plot/FigS7/FigS7_viol_2.pdf", paper = "a4")
  plot(p_patch_2 + plot_layout(nrow = 5, ncol = 4), height = 50)
  dev.off()
  
}

Heatmap <- function(df_p_3_gp, target_ant, sig_flag, plot_ant){
  
  target_flag <- df_p_3_gp$ID1_Antigen == target_ant & df_p_3_gp$ID2_Antigen == target_ant & df_p_3_gp$s_ctc != 0.5
  target_df <- df_p_3_gp[target_flag, ]
  
  ####### 20230112 Added ######
  tmp1 <- target_df %>% select(ID1, ID1_Cell_type_class)
  tmp2 <- target_df %>% select(ID2, ID2_Cell_type_class)
  colnames(tmp1) <- c("ID", "Cell_type_class")
  colnames(tmp2) <- c("ID", "Cell_type_class")
  tmp3 <- rbind(tmp1, tmp2) %>% distinct()
  rm_CTC <- tmp3 %>% group_by(Cell_type_class) %>% summarise(n = n()) %>% filter(n == 1) %>% .$Cell_type_class
  target_df_old <- target_df
  target_df <- target_df_old %>% filter(!ID1_Cell_type_class %in% rm_CTC & !ID2_Cell_type_class %in% rm_CTC) 
  #############################
  
  if (nrow(target_df) == 0){
    
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
      
      col_20_2 <- qualitative_hcl(length(unique(annot$Cell_type_class)), c = 50, l = 70)
      
      if (sig_flag == "*"){
        
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
      
      
      if (target_ant %in% plot_ant){
        
        svg(paste0("~/MOCCS_paper_public/plot/Fig3/Fig3D/Fig3D_heatmap_k_sim_jaccard_", target_ant, ".svg"),
            width = 8, height = 7.2)
        NMF::aheatmap(df_heatmap_2, Rowv = NA, Colv = NA,
                      labRow = NA, labCol = NA,
                      annCol = annot, annRow = annot,
                      annColors = list(col_20_2),
                      main = mtxt,
                      color = "Reds",
                      width = 0.3, height = 0.3)
        dev.off()
        
        pdf(paste0("~/MOCCS_paper_public/plot/Fig3/Fig3D/Fig3D_heatmap_k_sim_jaccard_", target_ant, ".pdf"),
            width = 8, height = 7.2)
        NMF::aheatmap(df_heatmap_2, Rowv = NA, Colv = NA,
                      labRow = NA, labCol = NA,
                      annCol = annot, annRow = annot,
                      annColors = list(col_20_2),
                      main = mtxt,
                      color = "Reds",
                      width = 0.3, height = 0.3)
        dev.off()
      }
      
    } 
    
  }
  
}

CircleChart <- function(sig_flag_list){
  
  ############################################################################
  # 2023/02/04 added
  sig_flag_list_2 <- sig_flag_list
  sig_num <- sum(sig_flag_list_2 == "*")
  sig_ratio <- sig_num / length(sig_flag_list_2)
  nonsig_num <- sum(sig_flag_list_2 == "")
  non_sig_ratio <- nonsig_num / length(sig_flag_list_2)
  null_num <- sum(sig_flag_list_2 == "NULL")
  null_ratio <- null_num / length(sig_flag_list_2)
  
  df_circ_2 = data.frame(category  =  c("   Cell-type dependent \n    TFs",
                                        "     Cell-type \n    non-dependent TFs",
                                        "Null"),
                         rate　　   =  c(round(sig_ratio, 2),
                                       round(non_sig_ratio, 2),
                                       round(null_ratio, 2))) %>% 
    #arrange(desc(category)) %>% 
    mutate(position = cumsum(rate) - rate/2)
  
  p_circ_2 <- ggplot(df_circ_2, aes(x = "", y = rate)) + 
    geom_bar(stat = "identity", color = "black", fill = c( "black", "white", "gray")) +
    coord_polar(theta = "y") +
    theme_void() + 
    geom_text(aes(y = position,  label = category), 
              size = 1.8,  color = c( "white","black", "black"),  vjust = 0) + 
    geom_text(aes(y = position,  label = paste(rate * 100, "%")), 
              size = 1.8,  color = c( "white", "black", "black"),  vjust = 2) +
    geom_text(aes(y = position,  label = c(paste0(sig_num, "/", length(sig_flag_list_2)),
                                           paste0((length(sig_flag_list_2) - sig_num), "/", length(sig_flag_list_2)),
                                           paste0(null_num, "/", length(sig_flag_list_2)))), 
              size = 1.8,  color = c( "white","black", "black"),  vjust = 4)
  
  ggsave(paste0("~/MOCCS_paper_public/plot/Fig3/Fig3E/Fig3E_circle_chart_alltf.svg"),
         plot = p_circ_2, width = 10)
  ggsave(paste0("~/MOCCS_paper_public/plot/Fig3/Fig3E/Fig3E_circle_chart_alltf.pdf"),
         plot = p_circ_2, width = 10)
  ############################################################################
  
  # old version
  sig_flag_list_2 <- sig_flag_list[!sig_flag_list == "NULL"]
  sig_num <- sum(sig_flag_list_2 == "*")
  sig_ratio <- sig_num / length(sig_flag_list_2)
  non_sig_ratio <- (length(sig_flag_list_2) - sig_num) / length(sig_flag_list_2)
  
  df_circ_1 = data.frame(category  =  c("   Cell-type dependent \n    TFs",
                                        "     Cell-type \n    non-dependent TFs"),
                         rate　　   =  c(round(sig_ratio, 2),
                                       round(non_sig_ratio, 2))) %>% 
    arrange(desc(category)) %>% 
    mutate(position = cumsum(rate) - rate/2)
  
  p_circ_1 <- ggplot(df_circ_1, aes(x = "", y = rate)) + 
    geom_bar(stat = "identity", color = "black", fill = c("black", "white")) +
    coord_polar(theta = "y") +
    theme_void() + 
    geom_text(aes(y = position,  label = category), 
              size = 4.5,  color = c("white", "black"),  vjust = 0) + 
    geom_text(aes(y = position,  label = paste(rate * 100, "%")), 
              size = 4.5,  color = c("white", "black"),  vjust = 2) +
    geom_text(aes(y = position,  label = c(paste0(sig_num, "/", length(sig_flag_list_2)),
                                           paste0((length(sig_flag_list_2) - sig_num), "/", length(sig_flag_list_2)))), 
              size = 4.5,  color = c("white", "black"),  vjust = 4)
  
  ggsave(paste0("~/MOCCS_paper_public/plot/Fig3/Fig3E/Fig3E_circle_chart.svg"),
         plot = p_circ_1, width = 10)
  ggsave(paste0("~/MOCCS_paper_public/plot/Fig3/Fig3E/Fig3E_circle_chart.pdf"),
         plot = p_circ_1, width = 10)
  
  ############################################################################
  
}

Make_TF_table <- function(ant_list){
  
  sig_flag_list <- readRDS("~/MOCCS_paper_public/data/Fig3/obj/sig_flag_list.rds")
  sig_flag_list_2 <- sig_flag_list[!sig_flag_list == "NULL"]
  
  sig_symbol_list <- ant_list[sig_flag_list == "*"]
  sig_symbol_list_2 <- ant_list[sig_flag_list == ""]
  sig_symbol_list_3 <- ant_list[sig_flag_list == "NULL"]
  saveRDS(sig_symbol_list, "~/MOCCS_paper_public/data/Fig3/obj/sig_symbol_list.rds")
  saveRDS(sig_symbol_list_2, "~/MOCCS_paper_public/data/Fig3/obj/sig_symbol_list_2.rds")
  saveRDS(sig_symbol_list_3, "~/MOCCS_paper_public/data/Fig3/obj/sig_symbol_list_3.rds")
  sig_symbol_list <- readRDS("~/MOCCS_paper_public/data/Fig3/obj/sig_symbol_list.rds")
  sig_symbol_list_2 <- readRDS("~/MOCCS_paper_public/data/Fig3/obj/sig_symbol_list_2.rds")
  sig_symbol_list_3 <- readRDS("~/MOCCS_paper_public/data/Fig3/obj/sig_symbol_list_3.rds")
  
  p_val_list <- readRDS("~/MOCCS_paper_public/data/Fig3/obj/p_val_list.rds")
  
  p_val_list_2 <- p_val_list[sig_flag_list == "*"]
  p_val_list_3 <- p_val_list[sig_flag_list == ""]
  p_val_list_4 <- rep("Not calculated", length(sig_symbol_list_3))
  
  df_ctd_tf <- cbind(sig_symbol_list, p_val_list_2)
  colnames(df_ctd_tf) <- c("Symbol", "p-value")
  df_ctd_tf_2 <- df_ctd_tf[order(as.double(df_ctd_tf[, "p-value"])), ]
  
  df_ctd_tf_3 <- cbind(sig_symbol_list_2, p_val_list_3)
  colnames(df_ctd_tf_3) <- c("Symbol", "p-value")
  df_ctd_tf_4 <- df_ctd_tf_3[order(as.double(df_ctd_tf_3[, "p-value"])), ]
  
  df_ctd_tf_5 <- cbind(sig_symbol_list_3, p_val_list_4)
  colnames(df_ctd_tf_5) <- c("Symbol", "p-value")
  
  write.table(df_ctd_tf_2[,"Symbol"], "~/MOCCS_paper_public/data/TableS1/Supplementary_table_S1_parts_1.txt", row.names = FALSE, col.names = FALSE)
  write.table(df_ctd_tf_2[,"p-value"], "~/MOCCS_paper_public/data/TableS1/Supplementary_table_S1_parts_2.txt", row.names = FALSE, col.names = FALSE)
  write.table(df_ctd_tf_4[,"Symbol"], "~/MOCCS_paper_public/data/TableS1/Supplementary_table_S1_parts_3.txt", row.names = FALSE, col.names = FALSE)
  write.table(df_ctd_tf_4[,"p-value"], "~/MOCCS_paper_public/data/TableS1/Supplementary_table_S1_parts_4.txt", row.names = FALSE, col.names = FALSE)
  write.table(df_ctd_tf_5[,"Symbol"], "~/MOCCS_paper_public/data/TableS1/Supplementary_table_S1_parts_5.txt", row.names = FALSE, col.names = FALSE)
  write.table(df_ctd_tf_5[,"p-value"], "~/MOCCS_paper_public/data/TableS1/Supplementary_table_S1_parts_6.txt", row.names = FALSE, col.names = FALSE)
  
}

# 2023/01/16 Added
# TF, TF family, number of analyzed samples, number of analyzed cell type class for the TF, cell type-dependent (Y/N)
Make_TF_table2 <- function(){
  df1 <- read_tsv("~/MOCCS_paper_public/data/TableS1/Supplementary_table_S1_parts_1.txt", col_names = F)
  df2 <- read_tsv("~/MOCCS_paper_public/data/TableS1/Supplementary_table_S1_parts_2.txt", col_names = F)
  df3 <- read_tsv("~/MOCCS_paper_public/data/TableS1/Supplementary_table_S1_parts_3.txt", col_names = F)
  df4 <- read_tsv("~/MOCCS_paper_public/data/TableS1/Supplementary_table_S1_parts_4.txt", col_names = F)
  df5 <- read_tsv("~/MOCCS_paper_public/data/TableS1/Supplementary_table_S1_parts_5.txt", col_names = F)
  df6 <- read_tsv("~/MOCCS_paper_public/data/TableS1/Supplementary_table_S1_parts_6.txt", col_names = F)
  
  TF_list <- c(df1$X1, df3$X1, df5$X1)
  pval_list <- c(df2$X1, df4$X1, df6$X1) %>% as.numeric()
  df7 <- tibble(TF = TF_list, pval = pval_list)
  
  ctc_num_list <- c()
  sample_num_list <- c()
  for (target_ant in TF_list) {
    target_flag <- df_p_3_gp$ID1_Antigen == target_ant & df_p_3_gp$ID2_Antigen == target_ant
    target_df <- df_p_3_gp[target_flag, ]
    
    ####### 20230112 Added ######
    tmp1 <- target_df %>% dplyr::select(ID1, ID1_Cell_type_class)
    tmp2 <- target_df %>% dplyr::select(ID2, ID2_Cell_type_class)
    colnames(tmp1) <- c("ID", "Cell_type_class")
    colnames(tmp2) <- c("ID", "Cell_type_class")
    tmp3 <- rbind(tmp1, tmp2) %>% distinct()
    rm_CTC <- tmp3 %>% group_by(Cell_type_class) %>% summarise(n = n()) %>% filter(n == 1) %>% .$Cell_type_class
    target_df_old <- target_df
    target_df <- target_df_old %>% filter(!ID1_Cell_type_class %in% rm_CTC & !ID2_Cell_type_class %in% rm_CTC) 
    #############################
    
    tgt_sample_num <- union(unique(target_df$ID1), unique(target_df$ID2)) %>% length()
    sample_num_list <- c(sample_num_list, tgt_sample_num)
    
    tgt_ctc_num <- union(unique(target_df$ID1_Cell_type_class), unique(target_df$ID2_Cell_type_class)) %>% length()
    ctc_num_list <- c(ctc_num_list, tgt_ctc_num)
  }
  
  df8 <- df7 %>% mutate(ctc_num = ctc_num_list)
  df9 <- df8 %>% mutate(ctc_depedent = ifelse(pval < 0.05, "Y", "N"))
  df10 <- df9 %>% mutate(sample_num = sample_num_list)
  
  # Add TF family
  TF_fam <- read_tsv("~/MOCCS_paper_public/data/Fig3/TF_Information.txt")
  Antigen_anno <- TF_list
  Antigen_TF_fam <- TF_fam$TF_Name %>% as.character() %>% unique()
  Antigen_share <- intersect(Antigen_anno, Antigen_TF_fam)
  TF_fam2 <- TF_fam %>% filter(TF_Name %in% Antigen_share) %>% select(TF_Name, Family_Name) %>% distinct()
  colnames(TF_fam2) <- c("TF", "TF_family")
  df10_fam <- df10 %>% left_join(TF_fam2, by = "TF") %>% select(TF, TF_family, sample_num, ctc_num, pval, ctc_depedent)
  
  write_tsv(df10_fam, "~/MOCCS_paper_public/data/TableS1/Supplementary_table_S1_parts_7.tsv")
}