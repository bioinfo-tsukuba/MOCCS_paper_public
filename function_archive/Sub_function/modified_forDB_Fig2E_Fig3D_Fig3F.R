# k-sim pearson files -----
list_files <- list.files("/Users/saeko/MOCCS_paper_public/data/Fig2/obj/k_sim_pearson_mat", pattern = "k_sim_pearson_mat_")
TF_list_files <- gsub("k_sim_pearson_mat_", "", list_files)
TF_list_files2 <- gsub(".rds", "", TF_list_files)

# all TFs
TF_list <- readRDS("/Users/saeko/Documents/MOCCS_SAECOM/Antigen_list_hg38_soft.rds")


# Read data frames ----
df_p_3 <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_3.rds")
df_p_3_gp <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_3_gp.rds")

all_ns <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/all_ns.rds")
flag_1 <- !(df_p_3_gp$ID1 %in% all_ns | df_p_3_gp$ID2 %in% all_ns)
flag_5 <- !(df_p_3$ID1 %in% all_ns | df_p_3$ID2 %in% all_ns) # flag without all non-significant k-mer samples

# Fig. 2E -----
source("~/MOCCS_paper_public/function/sub_function/TF_graph_2_db.R")
#TF_graph_2_db(df_p_3[flag_5, ], receiver_tf_list = c("JUN", "FOS", "FOXF1", "ELK1")) # Take several hours
TF_graph_2_db(df_p_3[flag_5, ], receiver_tf_list = TF_list[11:15]) 
source("~/MOCCS_paper_public/function/sub_function/modified_forDB_TF_graph_sel_db.R")
TF_graph_sel_db(ctc_spe = NULL, df_p_3[flag_5, ], receiver_tf_list = TF_list_files2)



# k-sim jaccard files -----
list_files <- list.files("/Users/saeko/MOCCS_paper_public/data/Fig3/obj/k_sim_jaccard_mat", pattern = "k_sim_jaccard_mat_")
TF_list_files <- gsub("k_sim_jaccard_mat_", "", list_files)
TF_list_files2 <- gsub(".rds", "", TF_list_files)

# Fig. 3F ------
source("~/MOCCS_paper_public/function/sub_function/TF_graph_3_db.R")
TF_graph_3_db(df_p_3[flag_5, ])
source("~/MOCCS_paper_public/function/sub_function/modified_forDB_TF_graph_sel_3_db.R")
for (target_TF in TF_list_files2) {
  print(target_TF)
  TF_graph_sel_3_db(a1 = target_TF, t_ctc = NULL, t_num = 15)
}

#TF_graph_sel_3_db(a1 = "JUN", t_ctc = NULL, t_num = 15)
#TF_graph_sel_3_db(a1 = "GATA2", t_ctc = NULL, t_num = 15)

# Fig. 3D Heatmap ------
source("~/MOCCS_paper_public/function/sub_function/Heatmap_db.R")
Heatmap_db(df_p_3_gp[flag_1, ],
           target_ant = "FOS")
Heatmap_db(df_p_3_gp[flag_1, ],
           target_ant = "JUN")
Heatmap_db(df_p_3_gp[flag_1, ],
           target_ant = "GATA2")
Heatmap_db(df_p_3_gp[flag_1, ],
           target_ant = "MYC")


Antigen_list <- readRDS("/Users/saeko/Documents/MOCCS_SAECOM/data/Antigen_list_hg38_soft.rds")
for (target_Antigen in Antigen_list) {
  if(file.exists(paste0("/Users/saeko/MOCCS_paper_public/plot/FigS6/heatmap_k_sim_jaccard_", target_Antigen, ".png"))==FALSE){
    print(target_Antigen)
    Heatmap_db(df_p_3_gp[flag_1, ],
               target_ant = target_Antigen)
  }
 
}
