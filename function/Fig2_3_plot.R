Fig2_3_plot <- function(data_path,
                        calc_opt = TRUE){
  
  library(dplyr)
  
  # Get family annotation of all sample pairs
  # Get annotation of all sample pairs
  # Calculate k-sim 1/2 between all sample pairs
  if (calc_opt){

    source("~/MOCCS-DB_paper/function/sub_function/Read_df.R")
    df_raw <- Read_df()
    saveRDS(df_raw, "~/MOCCS-DB_paper/results/df_raw.rds")
    
    source("~/MOCCS-DB_paper/function/sub_function/Annot_DBDs_of_CIS_BP.R")
    df_fam <- Annot_DBDs_of_CIS_BP(df_raw)
    saveRDS(df_fam, "~/MOCCS-DB_paper/results/df_fam.rds")
    
    source("~/MOCCS-DB_paper/function/sub_function/Annot_pairs.R")
    df_p_1 <- Annot_pairs(df_fam)
    saveRDS(df_p_1, "~/MOCCS-DB_paper/results/df_p_1.rds")
    
    source("~/MOCCS-DB_paper/function/sub_function/Calc_pairs.R")
    df_p_2 <- Calc_pairs(df_raw)
    saveRDS(df_p_2, "~/MOCCS-DB_paper/results/df_p_2.rds")
    
    # non-sig version
    # source("~/MOCCS-DB_paper/function/sub_function/Calc_pairs_2.R")
    # df_p_2_remain_non_sig <- Calc_pairs_2(df_raw)
    # saveRDS(df_p_2_remain_non_sig, "~/MOCCS-DB_paper/results/df_p_2_remain_non_sig.rds")
    
    # Join
    df_p_1 %>%
      mutate(ID_pair = paste0(ID1, ID2)) -> df_p_1_2
    df_p_2 %>%
      mutate(ID_pair = paste0(ID1, ID2)) %>%
      dplyr::select(ID_pair, k_sim_1, k_sim_2) -> df_p_2_2
    inner_join(df_p_1_2, df_p_2_2, by = "ID_pair") -> df_p_3
    saveRDS(df_p_3, "~/MOCCS-DB_paper/results/df_p_3.rds")
    
    source("~/MOCCS-DB_paper/function/sub_function/Group_df.R")
    df_p_3_gp <- Group_df(df_p_3)
    saveRDS(df_p_3_gp, "~/MOCCS-DB_paper/results/df_p_3_gp.rds")
    
    source("~/MOCCS-DB_paper/function/sub_function/Collect_poi.R")
    df_poi <- Collect_poi()
    saveRDS(df_poi, "~/MOCCS-DB_paper/results/df_poi.rds")
    
  } else {
    
    df_raw <- readRDS("~/MOCCS-DB_paper/results/df_raw.rds")
    df_fam <- readRDS("~/MOCCS-DB_paper/results/df_fam.rds")
    df_p_1 <- readRDS("~/MOCCS-DB_paper/results/df_p_1.rds")
    df_p_2 <- readRDS("~/MOCCS-DB_paper/results/df_p_2.rds")
    df_p_3 <- readRDS("~/MOCCS-DB_paper/results/df_p_3.rds")
    df_p_3_gp <- readRDS("~/MOCCS-DB_paper/results/df_p_3_gp.rds")
    df_poi <- readRDS("~/MOCCS-DB_paper/results/df_poi.rds")
    
  }

  rm(list = ls())
  df_p_3_gp <- readRDS("~/MOCCS-DB_paper/results/df_p_3_gp.rds")
  df_poi <- readRDS("~/MOCCS-DB_paper/results/df_poi.rds")
  # Fig. 2B and S4: Comparison with poi
  source("~/MOCCS-DB_paper/function/sub_function/Compare_ksim_poi.R")
  Compare_ksim_poi(df_p_3_gp, df_poi)
  
  # Fig. 2C and 3B: UMAP
  source("~/MOCCS-DB_paper/function/sub_function/Umap_df.R")
  UMAP_df(df_raw, df_fam)
  
  # Fig 2D, 3C and S5: Permutation and Chi squared tests
  source("~/MOCCS-DB_paper/function/sub_function/Top_k_sim.R")
  Top_k_sim(df_p_1, df_p_2, calc_opt = FALSE)
  
  # Fig. 2E: TF graph
  source("~/MOCCS-DB_paper/function/sub_function/TF_graph.R")
  TF_graph(df_p_3)
  
  # Fig. 3X
  source("~/MOCCS-DB_paper/function/sub_function/TF_graph_2.R")
  TF_graph_2(df_p_3)
  source("~/MOCCS-DB_paper/function/sub_function/TF_graph_3.R")
  TF_graph_3(df_p_3) # k-sim 1
  source("~/MOCCS-DB_paper/function/sub_function/TF_graph_sel.R")
  TF_graph_sel(ctc_spe = NULL, df_p_3)
  source("~/MOCCS-DB_paper/function/sub_function/TF_graph_sel_2.R")
  TF_graph_sel_2(rt_list = c("FOS", "JUN"), t_ctc = "Breast")
  TF_graph_sel_2(rt_list = c("FOS", "JUN"), t_ctc = "Cardiovascular")
  source("~/MOCCS-DB_paper/function/sub_function/TF_graph_sel_3.R")
  sig_symbol_list <- readRDS("~/MOCCS-DB_paper/results/sig_symbol_list.rds")
  sig_symbol_list_2 <- sig_symbol_list[!sig_symbol_list %in% c("AR", "FOS")]
  for (a1_ind in 1:length(sig_symbol_list_2)){
    TF_graph_sel_3(a1 = sig_symbol_list_2[a1_ind], t_ctc = NULL, t_num = 15)
  }

  # Fig. 3D and 3E: Cell type comparison
  plot_ant <- c("FOS", "JUN", "CTCF", "GATA2", "MYC")
  source("~/MOCCS-DB_paper/function/sub_function/Compare_cell_type.R")
  Compare_cell_type(df_p_3_gp, plot_ant, skip_opt = TRUE)
  

  return()
  
}