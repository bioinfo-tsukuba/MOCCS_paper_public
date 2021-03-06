Fig2_3_plot <- function(calc_opt = FALSE){
  
  library(dplyr)
  
  if (calc_opt){

    source("~/MOCCS_paper_public/function/sub_function/Read_df.R")
    df_raw <- Read_df()
    # saveRDS(df_raw, "~/MOCCS_paper_public/data/Fig2/obj/df_raw.rds")
    
    source("~/MOCCS_paper_public/function/sub_function/Annot_DBDs_of_CIS_BP.R")
    df_fam <- Annot_DBDs_of_CIS_BP(df_raw)
    saveRDS(df_fam, "~/MOCCS_paper_public/data/Fig2/obj/df_fam.rds")
    
    source("~/MOCCS_paper_public/function/sub_function/Annot_pairs.R")
    df_p_1 <- Annot_pairs(df_fam)
    saveRDS(df_p_1, "~/MOCCS_paper_public/data/Fig2/obj/df_p_1.rds")
    
    source("~/MOCCS_paper_public/function/sub_function/Calc_pairs.R")
    df_p_2 <- Calc_pairs(df_raw)
    saveRDS(df_p_2, "~/MOCCS_paper_public/data/Fig2/obj/df_p_2.rds")
    
    # Join
    df_p_1 %>%
      mutate(ID_pair = paste0(ID1, ID2)) -> df_p_1_2
    df_p_2 %>%
      mutate(ID_pair = paste0(ID1, ID2)) %>%
      dplyr::select(ID_pair, k_sim_1, k_sim_2) -> df_p_2_2
    inner_join(df_p_1_2, df_p_2_2, by = "ID_pair") -> df_p_3
    saveRDS(df_p_3, "~/MOCCS_paper_public/data/Fig2/obj/df_p_3.rds")
    
    source("~/MOCCS_paper_public/function/sub_function/Group_df.R")
    df_p_3_gp <- Group_df(df_p_3)
    saveRDS(df_p_3_gp, "~/MOCCS_paper_public/data/Fig2/obj/df_p_3_gp.rds")
    
    source("~/MOCCS_paper_public/function/sub_function/Collect_poi.R")
    df_poi <- Collect_poi()
    # saveRDS(df_poi, "~/MOCCS_paper_public/data/Fig2/obj/df_poi.rds")
    
    # Get ID list whose k-mers are all non-significant
    source("~/MOCCS_paper_public/function/sub_function/Get_ans.R")
    all_ns <- Get_ans(df_raw)
    saveRDS(all_ns, "~/MOCCS_paper_public/data/Fig2/obj/all_ns.rds")
    
  } else {
    
    system("wget https://figshare.com/ndownloader/files/34692295 -P ~/MOCCS_paper_public/data/Fig2/obj")
    system("mv ~/MOCCS_paper_public/data/Fig2/obj/34692295 ~/MOCCS_paper_public/data/Fig2/obj/df_raw.rds")
    df_raw <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_raw.rds")
    df_fam <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_fam.rds")
    df_p_1 <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_1.rds")
    df_p_2 <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_2.rds")
    df_p_3 <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_3.rds")
    df_p_3_gp <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_3_gp.rds")
    system("wget https://figshare.com/ndownloader/files/34692271 -P ~/MOCCS_paper_public/data/Fig2/obj")
    system("mv ~/MOCCS_paper_public/data/Fig2/obj/34692271 ~/MOCCS_paper_public/data/Fig2/obj/df_poi.rds")
    df_poi <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_poi.rds")
    all_ns <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/all_ns.rds")
    
  }
  

  print("Fig. 2B and S4: Comparison with poi")
  
  # Fig. 2B and S4: Comparison with poi
  source("~/MOCCS_paper_public/function/sub_function/Compare_ksim_poi.R")
  flag_1 <- !(df_p_3_gp$ID1 %in% all_ns | df_p_3_gp$ID2 %in% all_ns)
  flag_2 <- !(df_poi$ID1 %in% all_ns | df_poi$ID2 %in% all_ns)
  Compare_ksim_poi(df_p_3_gp[flag_1, ], df_poi[flag_2, ])
  
  print("Fig. 2C and 3B: UMAP")
  
  # Fig. 2C and 3B: UMAP
  source("~/MOCCS_paper_public/function/sub_function/Umap_df.R")
  UMAP_df(df_raw, df_fam, all_ns)
  
  print("Fig 2D, 3C and S5: Permutation and Chi squared tests")
  
  # Fig 2D, 3C and S5: Permutation and Chi squared tests
  source("~/MOCCS_paper_public/function/sub_function/Top_k_sim.R")
  flag_3 <- !(df_p_1$ID1 %in% all_ns | df_p_1$ID2 %in% all_ns)
  flag_4 <- !(df_p_2$ID1 %in% all_ns | df_p_2$ID2 %in% all_ns)
  Top_k_sim(df_p_1[flag_3, ], df_p_2[flag_4, ], calc_opt = FALSE, perm_num = 1000)
  
  print("Fig. 2E")
  
  # Fig. 2E
  flag_5 <- !(df_p_3$ID1 %in% all_ns | df_p_3$ID2 %in% all_ns)
  if (calc_opt){
    source("~/MOCCS_paper_public/function/sub_function/TF_graph_2.R")
    TF_graph_2(df_p_3[flag_5, ])
  }
  source("~/MOCCS_paper_public/function/sub_function/TF_graph_sel.R")
  TF_graph_sel(ctc_spe = NULL, df_p_3[flag_5, ])
  
  print("Fig. 3D and 3E: Cell type comparison")
  
  # Fig. 3D and 3E: Cell type comparison
  source("~/MOCCS_paper_public/function/sub_function/Compare_cell_type.R")
  flag_1 <- !(df_p_3_gp$ID1 %in% all_ns | df_p_3_gp$ID2 %in% all_ns)
  Compare_cell_type(df_p_3_gp[flag_1, ],
                    plot_ant = c("FOS", "JUN", "GATA2", "MYC"))
  
  print("Fig. 3F")
  
  # Fig. 3F
  flag_5 <- !(df_p_3$ID1 %in% all_ns | df_p_3$ID2 %in% all_ns)
  if (calc_opt){
    source("~/MOCCS_paper_public/function/sub_function/TF_graph_3.R")
    TF_graph_3(df_p_3[flag_5, ])
  }
  source("~/MOCCS_paper_public/function/sub_function/TF_graph_sel_3.R")
  TF_graph_sel_3(a1 = "JUN", t_ctc = NULL, t_num = 15)
  TF_graph_sel_3(a1 = "GATA2", t_ctc = NULL, t_num = 15)
  

  return()
  
}