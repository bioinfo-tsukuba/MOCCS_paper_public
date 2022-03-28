Compare_ksim_poi <- function(df_p_3_gp, df_poi){
  
  library(ggplot2)
  library(dplyr)
  library(gridSVG)
  
  ## Preprocessing

  # replace NA
  df_p_3_gp[is.na(df_p_3_gp)] <- 0
  
  # Exclude CTCF
  non_CTCF_flag <- !(df_p_3_gp$ID1_Antigen == "CTCF" | df_p_3_gp$ID2_Antigen == "CTCF")
  df_p_3_gp_nc <- df_p_3_gp[non_CTCF_flag, c("k_sim_1", "k_sim_2", "s_ctc", "s_ant", "s_ant_f", "s_ct")]
  df_poi[non_CTCF_flag, ] %>% mutate(n1_by_n1all_plus_n2_by_n2all_norm = n1_by_n1all_plus_n2_by_n2all/2) -> df_poi_nc
  
  at_unk_flag <- df_p_3_gp_nc$s_ant == 0.5
  f_unk_flag <- df_p_3_gp_nc$s_ant_f == 0.5
  ctc_unk_flag<- df_p_3_gp_nc$s_ctc == 0.5
  ct_unk_flag <- df_p_3_gp_nc$s_ct == 0.5
  
  # Combine
  df_p_3_poi_nc <- cbind(df_p_3_gp_nc, df_poi_nc[, "n1_by_n1all_plus_n2_by_n2all_norm"])
  
  df_p_3_poi_nc %>% mutate(Match_ant = case_when(s_ant == 1 ~ "Same antigen",
                                                 s_ant == 0.5 ~ "Unknown antigen",
                                                 s_ant == 0 ~ "Different antigen")) %>%
                    mutate(Match_fam = case_when(s_ant_f == 1 ~ "Same family",
                                                 s_ant_f == 0.5 ~ "Unknown family",
                                                 s_ant_f == 0 ~ "Different family")) %>%
                    mutate(Match_ctc = case_when(s_ctc == 1 ~ "Same cell type class",
                                                 s_ctc == 0.5 ~ "Unknown cell type class",
                                                 s_ctc == 0 ~ "Different cell type class")) %>%
                    mutate(Match_ct = case_when(s_ct == 1 ~ "Same cell type",
                                                 s_ct == 0.5 ~ "Unknown cell type",
                                                 s_ct == 0 ~ "Different cell type")) -> df_p_3_poi_nc_2
  df_p_3_poi_nc_2 %>% mutate(sim = k_sim_1) %>%
    mutate(Similarity = "k-sim 1") -> df_p_3_poi_nc_2_tmp_1
  df_p_3_poi_nc_2 %>% mutate(sim = k_sim_2) %>%
    mutate(Similarity = "k-sim 2") -> df_p_3_poi_nc_2_tmp_2
  df_p_3_poi_nc_2 %>% mutate(sim = n1_by_n1all_plus_n2_by_n2all_norm) %>%
    mutate(Similarity = "poi") -> df_p_3_poi_nc_2_tmp_3
  df_p_3_poi_nc_2_tmp_1 %>% rbind(df_p_3_poi_nc_2_tmp_2) %>%
    rbind(df_p_3_poi_nc_2_tmp_3) -> df_gg
  
  #saveRDS(df_gg, "~/MOCCS_paper_public/data/Fig2/obj/df_gg.rds")
  #df_gg <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_gg.rds")
  
  # Figure 2B
  group_vec <- c("A", "B", "C", "D", "E", "F", "G")
  group_mat <- rbind(c("Same antigen", "Same family", "Same cell type class", "Same cell type"), 
                     c("Same antigen", "N", "Same cell type class", "N"),　　　　　　　　　　　　
                     c("Same antigen", "N", "N", "N"),                                           
                     c("Different antigen", "Same family", "Same cell type class", "N"),         
                     c("Different antigen", "Same family", "N", "N"),                            
                     c("N", "N", "Same cell type class", "N"),
                     c("Different antigen", "Different family", "Different cell type class", "Different cell type"))
  file_pre <- c("main")
  df_gg_2 <- Plot_group_1(df_gg, group_vec, group_mat)

  # Figure 4SA
  group_vec <- c("C", "E", "G", "H", "I", "J", "K", "L", "M", "N")
  group_mat <- rbind(c("Same antigen", "N", "N", "N"),
                     c("Different antigen", "Same family", "N", "N"),
                     c("Different antigen", "Different family", "Different cell type class", "Different cell type"),
                     c("N", "N", "Same cell type class", "Different cell type"),
                     c("Different antigen", "Same family", "N", "Same cell type"),
                     c("N", "Same family", "N", "N"),
                     c("N", "N", "N", "Same cell type"),
                     c("N", "Same family", "Same cell type class", "N"),
                     c("N", "Same family", "N", "Same cell type"),　　　　　　　　　　　　
                     c("Same antigen", "N", "Same cell type class", "Same cell type"))
  df_gg_2_int <- Plot_group_2(df_gg, group_vec, group_mat)
  
  # Internal
  group_vec <- c("A", "B", "C", "D", "E", "F", "G",
                 "H", "I", "J", "K", "L", "M", "N")
  group_mat <- rbind(c("Same antigen", "Same family", "Same cell type class", "Same cell type"), 
                     c("Same antigen", "N", "Same cell type class", "N"),　　　　　　　　　　　　
                     c("Same antigen", "N", "N", "N"),                                           
                     c("Different antigen", "Same family", "Same cell type class", "N"),         
                     c("Different antigen", "Same family", "N", "N"),                            
                     c("N", "N", "Same cell type class", "N"),
                     c("Different antigen", "Different family", "Different cell type class", "Different cell type"),
                     c("N", "N", "Same cell type class", "Different cell type"),
                     c("Different antigen", "Same family", "N", "Same cell type"),
                     c("N", "Same family", "N", "N"),
                     c("N", "N", "N", "Same cell type"),
                     c("N", "Same family", "Same cell type class", "N"),
                     c("N", "Same family", "N", "Same cell type"),
                     c("Same antigen", "N", "Same cell type class", "Same cell type"))
  df_gg_2_2 <- Plot_group_3(df_gg, group_vec, group_mat)
  #saveRDS(df_gg_2_2, "~/MOCCS_paper_public/data/Fig2/obj/df_gg_2_2.rds")
  #df_gg_2_2 <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_gg_2_2.rds")
  
  Plot_group_2d(df_gg_2_2)
  
  # U test 1
  group_vec <- c("C", "J", "F", "M", "E", "H")
  t_sim_vec <- c("k-sim 1", "k-sim 2", "poi")
  nc_g <- "G"
  U_test_group(df_gg_2_2, group_vec, t_sim_vec, nc_g)
  
  # U test 2
  group_vec <- c("N", "B")
  t_sim_vec <- c("k-sim 1", "k-sim 2", "poi")
  nc_g <- "C"
  U_test_group(df_gg_2_2, group_vec, t_sim_vec, nc_g)
  
  # U test 3
  group_vec <- c("L", "K")
  t_sim_vec <- c("k-sim 1", "k-sim 2", "poi")
  nc_g <- "J"
  U_test_group(df_gg_2_2, group_vec, t_sim_vec, nc_g)
  
  # U test 4
  group_vec <- c("I", "D")
  t_sim_vec <- c("k-sim 1", "k-sim 2", "poi")
  nc_g <- "E"
  U_test_group(df_gg_2_2, group_vec, t_sim_vec, nc_g)
  
  ## Correlation
  group_vec <- c("A", "B", "C", "D", "E", "F", "G",
                 "H", "I", "J", "K", "L", "M", "N")
  sim_vec <- c("k_sim_1", "k_sim_2")
  cor_vec <- c()
  p_vec <- c()
  label_vec_plot <- c()
  label_vec_sim <- c()
  
  for (group_ind in 1:length(group_vec)){
    for (sim_ind in 1:length(sim_vec)){
      group_char <- group_vec[group_ind]
      sim_char <- sim_vec[sim_ind]
      df_gg_g <- df_gg_2_2[df_gg_2_2$Group == group_char, ]
      res_cor <- cor.test(df_gg_g[, sim_char], df_gg_g$n1_by_n1all_plus_n2_by_n2all_norm, method=c("pearson"))
      cor_vec <- append(cor_vec, res_cor$estimate)
      p_vec <- append(p_vec, res_cor$p.value)
      label_vec_plot <- append(label_vec_plot, group_char)
      label_vec_sim <- append(label_vec_sim, sim_char)
    }
  }

  df_tmp_cor_1 <- as_tibble(label_vec_plot)
  colnames(df_tmp_cor_1) <- c("Group")
  df_tmp_cor_1 %>% mutate(Similarity = label_vec_sim) %>%
                   mutate(Coefficient = cor_vec) %>% mutate(p = p_vec) -> df_cor_1
  
  print("P-values of correlation coefficient are...")
  print(df_cor_1$p)
  
  saveRDS(df_cor_1, "~/MOCCS_paper_public/data/Fig2/obj/df_cor_1.rds")
  df_cor_1 <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_cor_1.rds")
  
  df_cor_1$Similarity[df_cor_1$Similarity == "k_sim_1"] <- "2_k_sim_Jaccard"
  df_cor_1$Similarity[df_cor_1$Similarity == "k_sim_2"] <- "1_k_sim_Pearson"
  
  p_bar <- ggplot(df_cor_1, aes(x = Group, y = Coefficient, fill = Similarity)) +
    geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.5) +
    ylim(-0.15, 1) +
    theme_bw() +
    ylab("Separman's correlation coefficient") +
    xlab("Group of ChIP-seq sample") +
    ggtitle("Correlation between k-sim 1 or 2 and poi") +
    scale_fill_brewer(palette = "Set1") +
    theme(panel.grid.minor = element_blank()) +
    facet_wrap(~ Similarity, ncol = 1)

  #ggsave(paste0("~/MOCCS_paper_public/plot/FigS4/FigS4_k_sim_cor_bar.png"), plot = p_bar, width = 21, height = 4)  
  ggsave(paste0("~/MOCCS_paper_public/plot/FigS4/FigS4_k_sim_cor_bar.pdf"), plot = p_bar, width = 21, height = 4)  
  
  
  return()
  
}

Plot_group_1 <- function(df_gg, group_vec, group_mat){
  
  df_gg_2 <- c()
  
  for (group_ind_1 in 1:nrow(group_mat)){
    judge_vec <- group_mat[group_ind_1, ]
    if (judge_vec[1] != "N"){
      target_flag_1 <- df_gg$Match_ant == judge_vec[1]
    } else {
      target_flag_1 <- rep(TRUE, nrow(df_gg))
    }
    if (judge_vec[2] != "N"){
      target_flag_2 <- df_gg$Match_fam == judge_vec[2]
    } else {
      target_flag_2 <- rep(TRUE, nrow(df_gg))
    }
    if (judge_vec[3] != "N"){
      target_flag_3 <- df_gg$Match_ctc == judge_vec[3]
    } else {
      target_flag_3 <- rep(TRUE, nrow(df_gg))
    }
    if (judge_vec[4] != "N"){
      target_flag_4 <- df_gg$Match_ct == judge_vec[4]
    } else {
      target_flag_4 <- rep(TRUE, nrow(df_gg))
    }
    target_flag_5 <- target_flag_1 & target_flag_2 & target_flag_3 & target_flag_4
    df_gg[target_flag_5, ] %>% mutate(Group = group_vec[group_ind_1]) -> df_gg_tmp_1
    df_gg_2 %>% rbind(df_gg_tmp_1) -> df_gg_2
  }
  
  ## Distribution
  p_dist <- ggplot(df_gg_2, aes(x = Group, y = sim, fill = Similarity)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.07, position = position_dodge(width = 0.9), outlier.shape = NA) +
    theme_bw() +
    ylab("Index value") +
    xlab("") +
    ylim(-0.5, 1.2) +
    theme(panel.grid.minor = element_blank()) +
    ggtitle("Comparison of k-sim 1, 2 or poi") +
    facet_wrap(~ Similarity, nrow = 3)
  
  #ggsave(paste0("~/MOCCS_paper_public/plot/Fig2/Fig2B/Fig2B_k_sim_dist_main.png"), plot = p_dist, width = 21, height = 6)  
  ggsave(paste0("~/MOCCS_paper_public/plot/Fig2/Fig2B/Fig2B_k_sim_dist_main.pdf"), plot = p_dist, width = 7, height = 3)
  
  #svg(file=paste0("~/MOCCS_paper_public/plot/Fig2/Fig2B/Fig2B_k_sim_dist_main.svg"), width = 21, height = 6)
  #plot(p_dist)
  #dev.off()
  
  return(df_gg_2)
  
}

Plot_group_2 <- function(df_gg, group_vec, group_mat){
  
  df_gg_2 <- c()
  
  for (group_ind_1 in 1:nrow(group_mat)){
    judge_vec <- group_mat[group_ind_1, ]
    if (judge_vec[1] != "N"){
      target_flag_1 <- df_gg$Match_ant == judge_vec[1]
    } else {
      target_flag_1 <- rep(TRUE, nrow(df_gg))
    }
    if (judge_vec[2] != "N"){
      target_flag_2 <- df_gg$Match_fam == judge_vec[2]
    } else {
      target_flag_2 <- rep(TRUE, nrow(df_gg))
    }
    if (judge_vec[3] != "N"){
      target_flag_3 <- df_gg$Match_ctc == judge_vec[3]
    } else {
      target_flag_3 <- rep(TRUE, nrow(df_gg))
    }
    if (judge_vec[4] != "N"){
      target_flag_4 <- df_gg$Match_ct == judge_vec[4]
    } else {
      target_flag_4 <- rep(TRUE, nrow(df_gg))
    }
    target_flag_5 <- target_flag_1 & target_flag_2 & target_flag_3 & target_flag_4
    df_gg[target_flag_5, ] %>% mutate(Group = group_vec[group_ind_1]) -> df_gg_tmp_1
    df_gg_2 %>% rbind(df_gg_tmp_1) -> df_gg_2
  }
  
  ## Distribution
  p_dist <- ggplot(df_gg_2, aes(x = Group, y = sim, fill = Similarity)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.07, position = position_dodge(width = 0.9), outlier.shape = NA) +
    theme_bw() +
    ylab("Index value") +
    xlab("") +
    ylim(-0.5, 1.2) +
    theme(panel.grid.minor = element_blank()) +
    ggtitle("Comparison of k-sim 1, 2 or poi") +
    facet_wrap(~ Similarity, nrow = 3)
  
  #ggsave(paste0("~/MOCCS_paper_public/plot/FigS4/FigS4_k_sim_dist_suppl.png"), plot = p_dist, width = 21, height = 6)  
  ggsave(paste0("~/MOCCS_paper_public/plot/FigS4/FigS4_k_sim_dist_suppl.pdf"), plot = p_dist, width = 15, height = 4)  
  
  #svg(file=paste0("~/MOCCS_paper_public/plot/FigS4/FigS4_k_sim_dist_suppl.svg"), width = 21, height = 6)
  #plot(p_dist)
  #dev.off()
  
  return(df_gg_2)
  
}

Plot_group_3 <- function(df_gg, group_vec, group_mat){
  
  df_gg_2 <- c()
  
  for (group_ind_1 in 1:nrow(group_mat)){
    judge_vec <- group_mat[group_ind_1, ]
    if (judge_vec[1] != "N"){
      target_flag_1 <- df_gg$Match_ant == judge_vec[1]
    } else {
      target_flag_1 <- rep(TRUE, nrow(df_gg))
    }
    if (judge_vec[2] != "N"){
      target_flag_2 <- df_gg$Match_fam == judge_vec[2]
    } else {
      target_flag_2 <- rep(TRUE, nrow(df_gg))
    }
    if (judge_vec[3] != "N"){
      target_flag_3 <- df_gg$Match_ctc == judge_vec[3]
    } else {
      target_flag_3 <- rep(TRUE, nrow(df_gg))
    }
    if (judge_vec[4] != "N"){
      target_flag_4 <- df_gg$Match_ct == judge_vec[4]
    } else {
      target_flag_4 <- rep(TRUE, nrow(df_gg))
    }
    target_flag_5 <- target_flag_1 & target_flag_2 & target_flag_3 & target_flag_4
    df_gg[target_flag_5, ] %>% mutate(Group = group_vec[group_ind_1]) -> df_gg_tmp_1
    df_gg_2 %>% rbind(df_gg_tmp_1) -> df_gg_2
  }
  
  ## Distribution
  p_dist <- ggplot(df_gg_2, aes(x = Group, y = sim, fill = Similarity)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.07, position = position_dodge(width = 0.9), outlier.shape = NA) +
    theme_bw() +
    ylab("Index value") +
    xlab("") +
    ylim(-0.5, 1.2) +
    theme(panel.grid.minor = element_blank()) +
    ggtitle("Comparison of k-sim 1, 2 or poi") +
    facet_wrap(~ Similarity, nrow = 3)
  
  #ggsave(paste0("~/MOCCS_paper_public/plot/FigS4/FigS4_k_sim_dist_all.png"), plot = p_dist, width = 21, height = 6)  
  ggsave(paste0("~/MOCCS_paper_public/plot/FigS4/FigS4_k_sim_dist_all.pdf"), plot = p_dist, width = 7, height = 4)  
  
  #svg(file=paste0("~/MOCCS_paper_public/plot/FigS4/FigS4_k_sim_dist_all.svg"), width = 21, height = 6)
  #plot(p_dist)
  #dev.off()
  
  return(df_gg_2)
  
}

Plot_group_2d <- function(df_gg_2){
  
  ## 2d density
  p_2d_1 <- ggplot(df_gg_2[df_gg_2$Similarity != "poi" & df_gg_2$Similarity == "k-sim 1", ], aes(x = n1_by_n1all_plus_n2_by_n2all_norm, y = sim)) +
    geom_bin2d(bins = 50) +
    scale_fill_continuous(type = "viridis") +
    xlab("poi") +
    ylab("k-sim 1") +
    ggtitle("Two dimensional density plot of k-sim 1 and poi") +
    theme_bw() +
    facet_wrap(~ Group, nrow = 2) +
    # theme(axis.text.x = element_text(size = 6)) +
    theme(panel.grid.minor = element_blank())
  
  p_2d_2 <- ggplot(df_gg_2[df_gg_2$Similarity != "poi" & df_gg_2$Similarity == "k-sim 2", ], aes(x = n1_by_n1all_plus_n2_by_n2all_norm, y = sim)) +
    geom_bin2d(bins = 50) +
    scale_fill_continuous(type = "viridis") +
    xlab("poi") +
    ylab("k-sim 2") +
    ggtitle("Two dimensional density plot of k-sim 2 and poi") +
    theme_bw() +
    facet_wrap(~ Group, nrow = 2) +
    # theme(axis.text.x = element_text(size = 6)) +
    theme(panel.grid.minor = element_blank())
  
  #ggsave(paste0("~/MOCCS_paper_public/plot/FigS4/FigS4_k_sim_1_2d.png"), plot = p_2d_1, width = 18, height = 6)
  ggsave(paste0("~/MOCCS_paper_public/plot/FigS4/FigS4_k_sim_1_2d.pdf"), plot = p_2d_1, width = 18, height = 6)
  
  #ggsave(paste0("~/MOCCS_paper_public/plot/FigS4/FigS4_k_sim_2_2d.png"), plot = p_2d_2, width = 18, height = 6)
  ggsave(paste0("~/MOCCS_paper_public/plot/FigS4/FigS4_k_sim_2_2d.pdf"), plot = p_2d_2, width = 18, height = 6)
  
}

U_test_group <- function(df_gg_2, group_vec, t_sim_vec, nc_g){
  
  for (t_sim_ind in 1:length(t_sim_vec)){
    t_sim <- t_sim_vec[t_sim_ind]
    nc_flag <- df_gg_2$Similarity == t_sim & df_gg_2$Group == nc_g
    nc_d <- df_gg_2$sim[nc_flag]
    for (group_ind in 1:length(group_vec)){
      t_gp <- group_vec[group_ind]
      c_d <- df_gg_2$sim[df_gg_2$Group == group_vec[group_ind]]
      print(paste0(t_sim, " ", t_gp, ":"))
      print(exactRankTests::wilcox.exact(x = c_d, y = nc_d, paired = FALSE))
    }
  }
  
}