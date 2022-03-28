Top_k_sim <- function(df_p_1, df_p_2, calc_opt){
  
  library(dplyr)
  library(ggplot2)

  # Remove CTCF
  non_ctcf_flag <- df_p_1$ID1_Antigen != "CTCF" & df_p_1$ID2_Antigen != "CTCF"
  df_p_1_non_ctcf <- df_p_1[non_ctcf_flag, ]
  df_p_2_non_ctcf <- df_p_2[non_ctcf_flag, ]
  
  # Join
  df_p_1_non_ctcf %>%
    mutate(ID_pair = paste0(ID1, ID2)) -> df_p_1_2
  df_p_2_non_ctcf %>%
    mutate(ID_pair = paste0(ID1, ID2)) %>%
    dplyr::select(ID_pair, k_sim_1, k_sim_2) -> df_p_2_2
  inner_join(df_p_1_2, df_p_2_2, by = "ID_pair") -> df_p_3
  
  ## Original
  # Calculate same ratio
  res_orig_top <- Collect_top_n(df_p_3)
  res_orig_not_top <- Collect_not_top_n(df_p_3)
  saveRDS(res_orig_top, "~/MOCCS-DB_paper/results/res_orig_top.rds")
  # res_orig_top <- readRDS("~/MOCCS-DB_paper/results/res_orig_top.rds")
  
  print(res_orig_top$fam_lab_list_2[order(res_orig_top$fam_col_list_2, decreasing = TRUE)][1:5])
  print(unique(df_p_3[df_p_3$ID1_Family == "CUT,Homeodomain",]$ID1_Antigen))
  print(unique(df_p_3[df_p_3$ID1_Family == "Forkhead",]$ID1_Antigen))
  print(unique(df_p_3[df_p_3$ID1_Family == "Ets",]$ID1_Antigen))
  print(unique(df_p_3[df_p_3$ID1_Family == "AP-2",]$ID1_Antigen))
  print(unique(df_p_3[df_p_3$ID1_Family == "bZIP",]$ID1_Antigen))
  print(res_orig_top$ctc_lab_list_2[order(res_orig_top$ctc_col_list_2, decreasing = TRUE)][1:5])
  
  # Chi square test
  Heatmap_chi_sq_test(res_orig_top, res_orig_not_top)
  
  ## Permutation test
  if (calc_opt){
    perm_num <- 1000
    res_orig_top_perm_list <- vector("list", perm_num)
    for (perm_ind in 1:perm_num){
      df_p_3_perm <- Perm_df(df_p_3, perm_seed = perm_ind)
      res_orig_top_perm_list[[perm_ind]] <- Collect_top_n(df_p_3_perm)
    }
    saveRDS(res_orig_top_perm_list, "~/MOCCS-DB_paper/results/res_orig_top_perm_list.rds")
  }
  res_orig_top_perm_list <- readRDS("~/MOCCS-DB_paper/results/res_orig_top_perm_list.rds")
  
  print(paste0("Check TRUE:", res_orig_top$fam_col_list_2[df$fam_lab_list_2 == "No_annotation_2"] == 0))
  print(paste0("Check TRUE:", res_orig_top$ctc_col_list_2[df$ctc_lab_list_2 == "Others_2"] == 0))
  print(paste0("Check TRUE:", res_orig_top$ctc_col_list_2[df$ctc_lab_list_2 == "Uncalssified_2"] == 0))
  sum_chk_1 <- 0
  sum_chk_2 <- 0
  sum_chk_3 <- 0
  for (i in 1:length(res_orig_top_perm_list)){
    df <- res_orig_top_perm_list[[i]]
    sum_chk_1 <- sum_chk_1 + df$fam_col_list_2[df$fam_lab_list_2 == "No_annotation_2"]
    sum_chk_2 <- sum_chk_2 + df$ctc_col_list_2[df$ctc_lab_list_2 == "Others_2"]
    sum_chk_3 <- sum_chk_3 + df$ctc_col_list_2[df$ctc_lab_list_2 == "Uncalssified_2"]
  }
  print(paste0("Check TRUE:", sum_chk_1 == 0))
  print(paste0("Check TRUE:", sum_chk_2 == 0))
  print(paste0("Check TRUE:", sum_chk_3 == 0))

  # Collect and plot
  Perm_collect(res_orig_top_perm_list, res_orig_top)
  
  
  return()
  
}

Collect_top_n <- function(df_p_3){
  
  # Prepare pairs
  df_p_3 %>%
    dplyr::rename(tmp = ID1) %>%
    dplyr::rename(ID1 = ID2) %>%
    dplyr::rename(ID2 = tmp) %>%
    dplyr::rename(tmp = ID1_Cell_type_class) %>%
    dplyr::rename(ID1_Cell_type_class = ID2_Cell_type_class) %>%
    dplyr::rename(ID2_Cell_type_class = tmp) %>%
    dplyr::rename(tmp = ID1_Antigen) %>%
    dplyr::rename(ID1_Antigen = ID2_Antigen) %>%
    dplyr::rename(ID2_Antigen = tmp) %>%
    dplyr::rename(tmp = ID1_Family) %>%
    dplyr::rename(ID1_Family = ID2_Family) %>%
    dplyr::rename(ID2_Family = tmp) -> dftmp1
  dplyr::bind_rows(
    df_p_3,
    dftmp1
  ) -> dftmp2
  
  # Extract top similar pairs
  # dftmp2 %>% mutate(Num_ind = seq(1:nrow(dftmp2))) %>%
  dftmp2 %>%
    group_by(ID1, ID1_Cell_type_class, ID1_Antigen, ID1_Family) %>%
    top_n(n = 3, wt = k_sim_2) -> dftmp3
  
  # Collect
  dftmp4 <- dftmp3
  dftmp4$ID1_Cell_type_class[dftmp4$ID1_Cell_type_class == "Others"] <- "Others_2"
  dftmp4$ID1_Cell_type_class[dftmp4$ID1_Cell_type_class == "Unclassified"] <- "Uncalssified_2"
  dftmp4$ID1_Family[dftmp4$ID1_Family == "No_annotation"] <- "No_annotation_2"
  dftmp4$ID1_Family[dftmp4$ID1_Family == "Unknown"] <- "Unknown_2"
  id1_list <- unique(dftmp4$ID1)
  ant_col_list <- c()
  fam_col_list <- c()
  ctc_col_list <- c()
  ant_lab_list <- c()
  fam_lab_list <- c()
  ctc_lab_list <- c()
  for (id1_ind in 1:length(id1_list)){
    
    id1 <- id1_list[id1_ind]
    id1_flag <- id1 == dftmp4$ID1
    target_df_1 <- dftmp4[id1_flag, ]
    row_num_df <- nrow(target_df_1)
    
    ant_col <- sum(target_df_1$ID1_Antigen == target_df_1$ID2_Antigen) / row_num_df
    fam_col <- sum(target_df_1$ID1_Family == target_df_1$ID2_Family) / row_num_df
    ctc_col <- sum(target_df_1$ID1_Cell_type_class == target_df_1$ID2_Cell_type_class) / row_num_df
    ant_col_list <- append(ant_col_list, ant_col)
    fam_col_list <- append(fam_col_list, fam_col)
    ctc_col_list <- append(ctc_col_list, ctc_col)

    ant_lab <- unique(target_df_1$ID1_Antigen)
    fam_lab <- unique(target_df_1$ID1_Family)
    ctc_lab <- unique(target_df_1$ID1_Cell_type_class)
    ant_lab_list <- append(ant_lab_list, ant_lab)
    fam_lab_list <- append(fam_lab_list, fam_lab)
    ctc_lab_list <- append(ctc_lab_list, ctc_lab)
    
  }
  
  rev_ant_col_list <- rep(1, length(ant_col_list)) - ant_col_list
  rev_fam_col_list <- rep(1, length(fam_col_list)) - fam_col_list
  rev_ctc_col_list <- rep(1, length(ctc_col_list)) - ctc_col_list
  
  ant_col_list_norm <- ant_col_list / length(id1_list)
  fam_col_list_norm <- fam_col_list / length(id1_list)
  ctc_col_list_norm <- ctc_col_list / length(id1_list)
  rev_ant_col_list_norm <- rev_ant_col_list / length(id1_list)
  rev_fam_col_list_norm <- rev_fam_col_list / length(id1_list)
  rev_ctc_col_list_norm <- rev_ctc_col_list / length(id1_list)
  
  ant_col_list_2 <- c()
  ant_lab_list_2 <- c()
  uniq_ant_lab_list <- unique(ant_lab_list)
  uniq_ant_lab_num <- length(uniq_ant_lab_list)
  for (ant_ind in 1:uniq_ant_lab_num){
    ant_lab <- uniq_ant_lab_list[ant_ind]
    target_list <- ant_col_list[ant_lab_list == ant_lab]
    ant_lab_list_2 <- append(ant_lab_list_2, ant_lab)
    ant_col_list_2 <- append(ant_col_list_2, sum(target_list) / length(target_list))
  }

  fam_col_list_2 <- c()
  fam_lab_list_2 <- c()
  uniq_fam_lab_list <- unique(fam_lab_list)
  uniq_fam_lab_num <- length(uniq_fam_lab_list)
  for (fam_ind in 1:uniq_fam_lab_num){
    fam_lab <- uniq_fam_lab_list[fam_ind]
    target_list <- fam_col_list[fam_lab_list == fam_lab]
    fam_lab_list_2 <- append(fam_lab_list_2, fam_lab)
    fam_col_list_2 <- append(fam_col_list_2, sum(target_list) / length(target_list))
  }
  
  ctc_col_list_2 <- c()
  ctc_lab_list_2 <- c()
  uniq_ctc_lab_list <- unique(ctc_lab_list)
  uniq_ctc_lab_num <- length(uniq_ctc_lab_list)
  for (ctc_ind in 1:uniq_ctc_lab_num){
    ctc_lab <- uniq_ctc_lab_list[ctc_ind]
    target_list <- ctc_col_list[ctc_lab_list == ctc_lab]
    ctc_lab_list_2 <- append(ctc_lab_list_2, ctc_lab)
    ctc_col_list_2 <- append(ctc_col_list_2, sum(target_list) / length(target_list))
  }
  
  return(list(ant_col_list = ant_col_list,
              fam_col_list = fam_col_list,
              ctc_col_list = ctc_col_list,
              ant_lab_list = ant_lab_list,
              fam_lab_list = fam_lab_list,
              ctc_lab_list = ctc_lab_list,
              ant_col_list_2 = ant_col_list_2,
              fam_col_list_2 = fam_col_list_2,
              ctc_col_list_2 = ctc_col_list_2,
              ant_lab_list_2 = ant_lab_list_2,
              fam_lab_list_2 = fam_lab_list_2,
              ctc_lab_list_2 = ctc_lab_list_2,
              rev_ant_col_list = rev_ant_col_list,
              rev_fam_col_list = rev_fam_col_list,
              rev_ctc_col_list = rev_ctc_col_list,
              ant_col_list_norm = ant_col_list_norm,
              fam_col_list_norm = fam_col_list_norm,
              ctc_col_list_norm = ctc_col_list_norm,
              rev_ant_col_list_norm = rev_ant_col_list_norm,
              rev_fam_col_list_norm = rev_fam_col_list_norm,
              rev_ctc_col_list_norm = rev_ctc_col_list_norm))
  
}

Collect_not_top_n <- function(df_p_3){
  
  # Prepare pairs
  df_p_3 %>%
    dplyr::rename(tmp = ID1) %>%
    dplyr::rename(ID1 = ID2) %>%
    dplyr::rename(ID2 = tmp) %>%
    dplyr::rename(tmp = ID1_Cell_type_class) %>%
    dplyr::rename(ID1_Cell_type_class = ID2_Cell_type_class) %>%
    dplyr::rename(ID2_Cell_type_class = tmp) %>%
    dplyr::rename(tmp = ID1_Antigen) %>%
    dplyr::rename(ID1_Antigen = ID2_Antigen) %>%
    dplyr::rename(ID2_Antigen = tmp) %>%
    dplyr::rename(tmp = ID1_Family) %>%
    dplyr::rename(ID1_Family = ID2_Family) %>%
    dplyr::rename(ID2_Family = tmp) -> dftmp1
  dplyr::bind_rows(
    df_p_3,
    dftmp1
  ) -> dftmp2
  
  # Extract top similar pairs
  dftmp2 %>% mutate(tmp = seq(1:nrow(dftmp2))) -> dftmp2_2
  dftmp2_2 %>%
    group_by(ID1, ID1_Cell_type_class, ID1_Antigen, ID1_Family) %>%
    slice_max(n = 3, order_by = k_sim_2) -> dftmp3
  
  row_num_top_id <- dftmp3$tmp
  dftmp3_2 <- dftmp2[- row_num_top_id, ]
  
  # Collect
  dftmp4 <- dftmp3_2
  dftmp4$ID1_Cell_type_class[dftmp4$ID1_Cell_type_class == "Others"] <- "Others_2"
  dftmp4$ID1_Cell_type_class[dftmp4$ID1_Cell_type_class == "Unclassified"] <- "Uncalssified_2"
  dftmp4$ID1_Family[dftmp4$ID1_Family == "No_annotation"] <- "No_annotation_2"
  dftmp4$ID1_Family[dftmp4$ID1_Family == "Unknown"] <- "Unknown_2"
  id1_list <- unique(dftmp4$ID1)
  ant_col_list <- c()
  fam_col_list <- c()
  ctc_col_list <- c()
  ant_lab_list <- c()
  fam_lab_list <- c()
  ctc_lab_list <- c()
  for (id1_ind in 1:length(id1_list)){
    
    id1 <- id1_list[id1_ind]
    id1_flag <- id1 == dftmp4$ID1
    target_df_1 <- dftmp4[id1_flag, ]
    row_num_df <- nrow(target_df_1)
    
    ant_col <- sum(target_df_1$ID1_Antigen == target_df_1$ID2_Antigen) / row_num_df
    fam_col <- sum(target_df_1$ID1_Family == target_df_1$ID2_Family) / row_num_df
    ctc_col <- sum(target_df_1$ID1_Cell_type_class == target_df_1$ID2_Cell_type_class) / row_num_df
    ant_col_list <- append(ant_col_list, ant_col)
    fam_col_list <- append(fam_col_list, fam_col)
    ctc_col_list <- append(ctc_col_list, ctc_col)
    
    ant_lab <- unique(target_df_1$ID1_Antigen)
    fam_lab <- unique(target_df_1$ID1_Family)
    ctc_lab <- unique(target_df_1$ID1_Cell_type_class)
    ant_lab_list <- append(ant_lab_list, ant_lab)
    fam_lab_list <- append(fam_lab_list, fam_lab)
    ctc_lab_list <- append(ctc_lab_list, ctc_lab)
    
  }
  
  rev_ant_col_list <- rep(1, length(ant_col_list)) - ant_col_list
  rev_fam_col_list <- rep(1, length(fam_col_list)) - fam_col_list
  rev_ctc_col_list <- rep(1, length(ctc_col_list)) - ctc_col_list
  
  ant_col_list_norm <- ant_col_list / length(id1_list)
  fam_col_list_norm <- fam_col_list / length(id1_list)
  ctc_col_list_norm <- ctc_col_list / length(id1_list)
  rev_ant_col_list_norm <- rev_ant_col_list / length(id1_list)
  rev_fam_col_list_norm <- rev_fam_col_list / length(id1_list)
  rev_ctc_col_list_norm <- rev_ctc_col_list / length(id1_list)
  
  ant_col_list_2 <- c()
  ant_lab_list_2 <- c()
  uniq_ant_lab_list <- unique(ant_lab_list)
  uniq_ant_lab_num <- length(uniq_ant_lab_list)
  for (ant_ind in 1:uniq_ant_lab_num){
    ant_lab <- uniq_ant_lab_list[ant_ind]
    target_list <- ant_col_list[ant_lab_list == ant_lab]
    ant_lab_list_2 <- append(ant_lab_list_2, ant_lab)
    ant_col_list_2 <- append(ant_col_list_2, sum(target_list) / length(target_list))
  }
  
  fam_col_list_2 <- c()
  fam_lab_list_2 <- c()
  uniq_fam_lab_list <- unique(fam_lab_list)
  uniq_fam_lab_num <- length(uniq_fam_lab_list)
  for (fam_ind in 1:uniq_fam_lab_num){
    fam_lab <- uniq_fam_lab_list[fam_ind]
    target_list <- fam_col_list[fam_lab_list == fam_lab]
    fam_lab_list_2 <- append(fam_lab_list_2, fam_lab)
    fam_col_list_2 <- append(fam_col_list_2, sum(target_list) / length(target_list))
  }
  
  ctc_col_list_2 <- c()
  ctc_lab_list_2 <- c()
  uniq_ctc_lab_list <- unique(ctc_lab_list)
  uniq_ctc_lab_num <- length(uniq_ctc_lab_list)
  for (ctc_ind in 1:uniq_ctc_lab_num){
    ctc_lab <- uniq_ctc_lab_list[ctc_ind]
    target_list <- ctc_col_list[ctc_lab_list == ctc_lab]
    ctc_lab_list_2 <- append(ctc_lab_list_2, ctc_lab)
    ctc_col_list_2 <- append(ctc_col_list_2, sum(target_list) / length(target_list))
  }
  
  return(list(ant_col_list = ant_col_list,
              fam_col_list = fam_col_list,
              ctc_col_list = ctc_col_list,
              ant_lab_list = ant_lab_list,
              fam_lab_list = fam_lab_list,
              ctc_lab_list = ctc_lab_list,
              ant_col_list_2 = ant_col_list_2,
              fam_col_list_2 = fam_col_list_2,
              ctc_col_list_2 = ctc_col_list_2,
              ant_lab_list_2 = ant_lab_list_2,
              fam_lab_list_2 = fam_lab_list_2,
              ctc_lab_list_2 = ctc_lab_list_2,
              rev_ant_col_list = rev_ant_col_list,
              rev_fam_col_list = rev_fam_col_list,
              rev_ctc_col_list = rev_ctc_col_list,
              ant_col_list_norm = ant_col_list_norm,
              fam_col_list_norm = fam_col_list_norm,
              ctc_col_list_norm = ctc_col_list_norm,
              rev_ant_col_list_norm = rev_ant_col_list_norm,
              rev_fam_col_list_norm = rev_fam_col_list_norm,
              rev_ctc_col_list_norm = rev_ctc_col_list_norm))
  
}

Heatmap_chi_sq_test <- function(res_orig_top, res_orig_not_top){
  
  top_ant <- sum(res_orig_top$ant_col_list)
  not_top_ant <- sum(res_orig_not_top$ant_col_list)
  top_ant_rev <- sum(res_orig_top$rev_ant_col_list)
  not_top_ant_rev <- sum(res_orig_not_top$rev_ant_col_list)
  
  #############################################
  #
  #            top           not top
  #
  # same       top_hoge      not_top_hoge
  #
  # different  top_hoge_rev  not_top_hoge_rev
  #
  ant_mat <- rbind(c(top_ant, not_top_ant),
                   c(top_ant_rev, not_top_ant_rev))
  p_ant_chi <- chisq.test(ant_mat)$p.value
  print(p_ant_chi)
  ant_txt <- rbind(as.character(c(round(top_ant), round(not_top_ant))),
                   as.character(c(round(top_ant_rev), round(not_top_ant_rev))))
  
  annot_1 <- data.frame(label = c("top", "not top"))
  annot_2 <- data.frame(label = c("same", "different"))
  
  pdf("~/MOCCS-DB_paper/plot/Fig2/FigS5/FigS5_chi_ant.pdf")
  NMF::aheatmap(ant_mat, Rowv = NA, Colv = NA,
                labRow = NA, labCol = NA,
                annCol = annot_1, annRow = annot_2,
                annColors = "Pastel1",
                txt = ant_txt,
                main = c("Chi-squared test on antigen"),
                sub = paste0("p value is XX"),
                color = "Reds")
  dev.off()
  
  top_fam <- sum(res_orig_top$fam_col_list)
  not_top_fam <- sum(res_orig_not_top$fam_col_list)
  top_fam_rev <- sum(res_orig_top$rev_fam_col_list)
  not_top_fam_rev <- sum(res_orig_not_top$rev_fam_col_list)
  fam_mat <- rbind(c(top_fam, not_top_fam),
                   c(top_fam_rev, not_top_fam_rev))
  p_fam_chi <- chisq.test(fam_mat)$p.value
  print(p_fam_chi)
  fam_txt <- rbind(as.character(c(round(top_fam), round(not_top_fam))),
                   as.character(c(round(top_fam_rev), round(not_top_fam_rev))))
  
  pdf("~/MOCCS-DB_paper/plot/Fig2/FigS5/FigS5_chi_fam.pdf")
  NMF::aheatmap(fam_mat, Rowv = NA, Colv = NA,
                labRow = NA, labCol = NA,
                annCol = annot_1, annRow = annot_2,
                annColors = "Pastel1",
                txt = fam_txt,
                main = c("Chi-squared test on family"),
                sub = paste0("p value is XX"),
                color = "Reds")
  dev.off()
  
  top_ctc <- sum(res_orig_top$ctc_col_list)
  not_top_ctc <- sum(res_orig_not_top$ctc_col_list)
  top_ctc_rev <- sum(res_orig_top$rev_ctc_col_list)
  not_top_ctc_rev <- sum(res_orig_not_top$rev_ctc_col_list)
  ctc_mat <- rbind(c(top_ctc, not_top_ctc),
                   c(top_ctc_rev, not_top_ctc_rev))
  p_ctc_chi <- chisq.test(ctc_mat)$p.value
  print(p_ctc_chi)
  ctc_txt <- rbind(as.character(c(round(top_ctc), round(not_top_ctc))),
                   as.character(c(round(top_ctc_rev), round(not_top_ctc_rev))))
  
  pdf("~/MOCCS-DB_paper/plot/Fig2/FigS5/FigS5_chi_ctc.pdf")
  NMF::aheatmap(ctc_mat, Rowv = NA, Colv = NA,
                labRow = NA, labCol = NA,
                annCol = annot_1, annRow = annot_2,
                annColors = "Pastel1",
                txt = ctc_txt,
                main = c("Chi-squared test on cell type class"),
                sub = paste0("p value is XX"),
                color = "Reds")
  dev.off()
  
  return()
  
}

Perm_df <- function(df_p_3, perm_seed){
  
  set.seed(perm_seed)
  
  orig_k_sim_2 <- df_p_3$k_sim_2
  pair_num <- length(orig_k_sim_2)
  
  perm_ind <- sample(seq(1:pair_num), pair_num)
  perm_k_sim_2 <- orig_k_sim_2[perm_ind]
  
  df_p_3 %>% mutate(k_sim_2 = perm_k_sim_2) -> df_p_3_perm
  
  return(df_p_3_perm)
  
}

Perm_collect <- function(res_orig_top_perm_list, res_orig_top){
  
  perm_num <- length(res_orig_top_perm_list)
  
  ant_sum_list <- c()
  fam_sum_list <- c()
  ctc_sum_list <- c()
  for (perm_ind in 1:perm_num){
    ant_sum_list <- append(ant_sum_list, sum(res_orig_top_perm_list[[perm_ind]]$ant_col_list_2))
    fam_sum_list <- append(fam_sum_list, sum(res_orig_top_perm_list[[perm_ind]]$fam_col_list_2))
    ctc_sum_list <- append(ctc_sum_list, sum(res_orig_top_perm_list[[perm_ind]]$ctc_col_list_2))
  }
  ant_ecdf <- ecdf(ant_sum_list)
  ant_ecdf_p <- 1 - ant_ecdf(sum(res_orig_top$ant_col_list_2))
  fam_ecdf <- ecdf(fam_sum_list)
  fam_ecdf_p <- 1 - fam_ecdf(sum(res_orig_top$fam_col_list_2))
  ctc_ecdf <- ecdf(ctc_sum_list)
  ctc_ecdf_p <- 1 - ctc_ecdf(sum(res_orig_top$ctc_col_list_2))
  
  ant_sum_list_2 <- rep(0, length(res_orig_top_perm_list[[perm_num]]$ant_col_list_2))
  fam_sum_list_2 <- rep(0, length(res_orig_top_perm_list[[perm_num]]$fam_col_list_2))
  ctc_sum_list_2 <- rep(0, length(res_orig_top_perm_list[[perm_num]]$ctc_col_list_2))
  for (perm_ind in 1:perm_num){
    ant_ord <- order(res_orig_top_perm_list[[perm_ind]]$ant_lab_list_2)
    ant_lab_list_2_ord <- res_orig_top_perm_list[[perm_ind]]$ant_lab_list_2[ant_ord]
    ant_col_list_2_ord <- res_orig_top_perm_list[[perm_ind]]$ant_col_list_2[ant_ord]
    ant_sum_list_2 <- ant_sum_list_2 + ant_col_list_2_ord
    fam_ord <- order(res_orig_top_perm_list[[perm_ind]]$fam_lab_list_2)
    fam_lab_list_2_ord <- res_orig_top_perm_list[[perm_ind]]$fam_lab_list_2[fam_ord]
    fam_col_list_2_ord <- res_orig_top_perm_list[[perm_ind]]$fam_col_list_2[fam_ord]
    fam_sum_list_2 <- fam_sum_list_2 + fam_col_list_2_ord
    ctc_ord <- order(res_orig_top_perm_list[[perm_ind]]$ctc_lab_list_2)
    ctc_lab_list_2_ord <- res_orig_top_perm_list[[perm_ind]]$ctc_lab_list_2[ctc_ord]
    ctc_col_list_2_ord <- res_orig_top_perm_list[[perm_ind]]$ctc_col_list_2[ctc_ord]
    ctc_sum_list_2 <- ctc_sum_list_2 + ctc_col_list_2_ord
  }
  ant_sum_list_3_perm <- ant_sum_list_2 / perm_num
  fam_sum_list_3_perm <- fam_sum_list_2 / perm_num
  ctc_sum_list_3_perm <- ctc_sum_list_2 / perm_num  
 
  ant_ord <- order(res_orig_top$ant_lab_list_2)
  ant_lab_list_3_orig <- res_orig_top$ant_lab_list_2[ant_ord]
  ant_col_list_3_orig <- res_orig_top$ant_col_list_2[ant_ord]
  fam_ord <- order(res_orig_top$fam_lab_list_2)
  fam_lab_list_3_orig <- res_orig_top$fam_lab_list_2[fam_ord]
  fam_col_list_3_orig <- res_orig_top$fam_col_list_2[fam_ord]
  ctc_ord <- order(res_orig_top$ctc_lab_list_2)
  ctc_lab_list_3_orig <- res_orig_top$ctc_lab_list_2[ctc_ord]
  ctc_col_list_3_orig <- res_orig_top$ctc_col_list_2[ctc_ord]
  
  print(paste0("Check TRUE: ", sum(ant_lab_list_3_orig == ant_lab_list_2_ord) == length(ant_lab_list_3_orig)))
  print(paste0("Check TRUE: ", sum(fam_lab_list_3_orig == fam_lab_list_2_ord) == length(fam_lab_list_3_orig)))
  print(paste0("Check TRUE: ", sum(ctc_lab_list_3_orig == ctc_lab_list_2_ord) == length(ctc_lab_list_3_orig)))
  
  # ggplot
  df_l_1 <- tibble(
    Label = append(ant_lab_list_2_ord, ant_lab_list_3_orig),
    Category = append(rep("1. Permutated", length(ant_lab_list_2_ord)), rep("2. Original", length(ant_lab_list_3_orig))),
    Ratio = append(ant_sum_list_3_perm, ant_col_list_3_orig)
  )
  df_l_1 %>%
    ggplot(aes(Category, Ratio, group = Label)) + 
    geom_line(size = 0.2) + 
    geom_point() +
    ggtitle("Neighboring pairs") +
    xlab("") +
    ylab("Ratio of same antigen") +
    ylim(0, 1.2) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    # geom_text(x = 1.5, y = 1.2, label = "*") +
    geom_segment(x = 1, xend = 2, y = 1.1, yend = 1.1) +
    geom_segment(x = 1, xend = 1, y = 1.1, yend = 1.05) +
    geom_segment(x = 2, xend = 2, y = 1.1, yend = 1.05) -> p_l_1
  
  ggsave("~/MOCCS-DB_paper/plot/Fig2/Fig2D/Fig2D_ant.pdf", plot = p_l_1, width = 6)

  df_l_2 <- tibble(
    Label = append(fam_lab_list_2_ord, fam_lab_list_3_orig),
    Category = append(rep("1. Permutated", length(fam_lab_list_2_ord)), rep("2. Original", length(fam_lab_list_3_orig))),
    Ratio = append(fam_sum_list_3_perm, fam_col_list_3_orig)
  )
  df_l_2 %>%
    ggplot(aes(Category, Ratio, group = Label)) + 
    geom_line(size = 0.2) + 
    geom_point() +
    ggtitle("Neighboring pairs") +
    xlab("") +
    ylab("Ratio of same family") +
    ylim(0, 1.2) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    # geom_text(x = 1.5, y = 1.2, label = "*") +
    geom_segment(x = 1, xend = 2, y = 1.1, yend = 1.1) +
    geom_segment(x = 1, xend = 1, y = 1.1, yend = 1.05) +
    geom_segment(x = 2, xend = 2, y = 1.1, yend = 1.05) -> p_l_2
  
  ggsave("~/MOCCS-DB_paper/plot/Fig2/Fig2D/Fig2D_fam.pdf", plot = p_l_2, width = 6)
  
  df_l_3 <- tibble(
    Label = append(ctc_lab_list_2_ord, ctc_lab_list_3_orig),
    Category = append(rep("1. Permutated", length(ctc_lab_list_2_ord)), rep("2. Original", length(ctc_lab_list_3_orig))),
    Ratio = append(ctc_sum_list_3_perm, ctc_col_list_3_orig)
  )
  df_l_3 %>% mutate(Group = "Cell type class") %>%
    ggplot(aes(Category, Ratio, group = Label)) + 
    geom_line(size = 0.2) + 
    geom_point() +
    ggtitle("Neighboring pairs") +
    xlab("") +
    ylab("Ratio of same cell type class") +
    ylim(0, 1.2) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    facet_wrap(~ Group) +
    # geom_text(x = 1.5, y = 1.2, label = "*") +
    geom_segment(x = 1, xend = 2, y = 1.1, yend = 1.1) +
    geom_segment(x = 1, xend = 1, y = 1.1, yend = 1.05) +
    geom_segment(x = 2, xend = 2, y = 1.1, yend = 1.05) -> p_l_3

  ggsave("~/MOCCS-DB_paper/plot/Fig3/Fig3C/Fig3C_ctc.pdf", plot = p_l_3, width = 6)
  
  df_l_1 %>% mutate(Group = "Antigen") -> df_l_1_tmp
  df_l_2 %>% mutate(Group = "Family") -> df_l_2_tmp
  df_l_1_tmp %>% rbind(df_l_2_tmp) -> df_l_comb
  df_l_comb %>%
    ggplot(aes(Category, Ratio, group = Label)) + 
    geom_line(size = 0.2) + 
    geom_point() +
    ggtitle("Neighboring pairs") +
    xlab("") +
    ylab("Ratio of same annotation") +
    ylim(0, 1.2) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    facet_wrap(~ Group) +
    # geom_text(x = 1.5, y = 1.2, label = "*") +
    geom_segment(x = 1, xend = 2, y = 1.1, yend = 1.1) +
    geom_segment(x = 1, xend = 1, y = 1.1, yend = 1.05) +
    geom_segment(x = 2, xend = 2, y = 1.1, yend = 1.05) -> p_l_comb
  
  ggsave("~/MOCCS-DB_paper/plot/Fig2/Fig2D/Fig2D_comb.pdf", plot = p_l_comb, width = 6)
  
  
  return()
  
}