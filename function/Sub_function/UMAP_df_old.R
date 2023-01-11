UMAP_df <- function(df_raw, df_fam, all_ns){

  library(ggplot2)
  library(dplyr)

  set.seed(123)
  
  flag_1 <- !(df_raw$ID %in% all_ns)

  df_raw_sg <- df_raw[flag_1, ]
  df_raw_sg$MOCCS2score[df_raw_sg$q_value >= 0.05] <- 0
  
  df_raw_sg[, c("kmer", "MOCCS2score", "ID", "Antigen", "Cell_type_class")] %>%
    dplyr::group_by(ID) %>%
    dplyr::distinct(kmer, .keep_all = TRUE) %>%
    tidyr::pivot_wider(names_from = "kmer", values_from = "MOCCS2score" ) -> df_raw_2
  df_raw_2[is.na(df_raw_2)] <- 0

  df_raw_3 <- df_raw_2[df_raw_2$Antigen != "CTCF", ]
  df_raw_4 <- left_join(df_raw_3, df_fam[, c("ID", "Family")], by = "ID")
  df_raw_5 <- df_raw_4[, 4:(ncol(df_raw_4) - 1)]

  res.umap <- umap::umap(df_raw_5, metric = "pearson", spread = 10)
  
  as_tibble(res.umap$layout) %>%
    mutate(ID = df_raw_4$ID) %>%
    mutate(Antigen = df_raw_4$Antigen) %>%
    mutate(Cell_type_class = df_raw_4$Cell_type_class) %>%
    mutate(Family = df_raw_4$Family) -> df_umap
  colnames(df_umap)[1:2] <- c("UMAP1", "UMAP2")
  
  df_umap_2 <- c()
  annot_vec <- c("Antigen", "Family")
  for (i in 1:length(annot_vec)){
    if (annot_vec[i] == "Antigen"){
      df_umap %>% mutate(Annotation =  Antigen) %>% mutate(Plot = "Antigen") -> df_umap_tmp_1
    } else if (annot_vec[i] == "Family"){
      df_umap %>% mutate(Annotation =  Family) %>% mutate(Plot = "Family") -> df_umap_tmp_1
    }
    df_umap_2 <- rbind(df_umap_2, df_umap_tmp_1)
  }

  df_umap_3 <- c()
  annot_vec_2 <- c("Cell type class")
  for (i in 1:length(annot_vec_2)){
    if (annot_vec_2[i] == "Cell type class"){
      df_umap %>% mutate(Annotation =  Cell_type_class) %>% mutate(Plot = "Cell type class") -> df_umap_tmp_2
    }
    df_umap_3 <- rbind(df_umap_3, df_umap_tmp_2)
  }
  
  p_umap_1 <- ggplot2::ggplot(df_umap_2[df_umap_2$Plot == "Antigen", ], aes(x = UMAP1, y = UMAP2, color = Annotation)) +
    geom_point() +
    theme(aspect.ratio = 1.0) +
    theme_bw() +
    theme(legend.position = "none") +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank()) +
    facet_wrap(~ Plot)

  p_umap_2 <- ggplot2::ggplot(df_umap_2[df_umap_2$Plot == "Family", ], aes(x = UMAP1, y = UMAP2, color = Annotation)) +
    geom_point() +
    theme(aspect.ratio = 1.0) +
    theme_bw() +
    theme(legend.position = "none") +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank()) +
    facet_wrap(~ Plot)

  p_umap_3 <- ggplot2::ggplot(df_umap_3[df_umap_3$Plot == "Cell type class", ], aes(x = UMAP1, y = UMAP2, color = Annotation)) +
    geom_point() +
    theme(aspect.ratio = 1.0) +
    theme_bw() +
    theme(legend.position = "none") +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank()) +
    facet_wrap(~ Plot)
  
  #ggsave(paste0("~/MOCCS_paper_public/plot/Fig2/Fig2C/Fig2C_umap_ant.png"), plot = p_umap_1, width = 5)
  ggsave(paste0("~/MOCCS_paper_public/plot/Fig2/Fig2C/Fig2C_umap_ant.pdf"), plot = p_umap_1, width = 6, height = 6)

  #ggsave(paste0("~/MOCCS_paper_public/plot/Fig2/Fig2C/Fig2C_umap_fam.png"), plot = p_umap_2, width = 5)
  ggsave(paste0("~/MOCCS_paper_public/plot/Fig2/Fig2C/Fig2C_umap_fam.pdf"), plot = p_umap_2, width = 6, height = 6)
  
  #ggsave(paste0("~/MOCCS_paper_public/plot/Fig3/Fig3B/Fig3B_umap.png"), plot = p_umap_3, width = 5)
  ggsave(paste0("~/MOCCS_paper_public/plot/Fig3/Fig3B/Fig3B_umap.pdf"), plot = p_umap_3, width = 6, height = 6)
  

  return()
  
}