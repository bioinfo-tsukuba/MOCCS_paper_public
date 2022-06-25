TF_list <- readRDS("/Users/saeko/Documents/MOCCS_SAECOM/Antigen_list_hg38_soft.rds")

# Read data frames
df_p_3 <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_3.rds")
df_p_3_gp <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_3_gp.rds")

all_ns <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/all_ns.rds")
flag_1 <- !(df_p_3_gp$ID1 %in% all_ns | df_p_3_gp$ID2 %in% all_ns)
flag_5 <- !(df_p_3$ID1 %in% all_ns | df_p_3$ID2 %in% all_ns) # flag without all non-significant k-mer samples

# Fig. 2E
source("~/MOCCS_paper_public/function/sub_function/TF_graph_2_db.R")
#TF_graph_2_db(df_p_3[flag_5, ], receiver_tf_list = c("JUN", "FOS", "FOXF1", "ELK1")) # Take several hours
TF_graph_2_db(df_p_3[flag_5, ], receiver_tf_list = TF_list[1]) 
source("~/MOCCS_paper_public/function/sub_function/TF_graph_sel_db.R")
TF_graph_sel_db(ctc_spe = NULL, df_p_3[flag_5, ], receiver_tf_list = c("JUN", "FOS", "FOXF1", "ELK1"))

# Fig. 3F
source("~/MOCCS_paper_public/function/sub_function/TF_graph_3_db.R")
TF_graph_3_db(df_p_3[flag_5, ])
source("~/MOCCS_paper_public/function/sub_function/TF_graph_sel_3_db.R")
TF_graph_sel_3_db(a1 = "JUN", t_ctc = NULL, t_num = 15)
TF_graph_sel_3_db(a1 = "GATA2", t_ctc = NULL, t_num = 15)

# Fig. 3D Heatmap
source("~/MOCCS_paper_public/function/sub_function/Heatmap_db.R")
Heatmap_db(df_p_3_gp[flag_1, ],
           target_ant = "FOS")
Heatmap_db(df_p_3_gp[flag_1, ],
           target_ant = "JUN")
Heatmap_db(df_p_3_gp[flag_1, ],
           target_ant = "GATA2")
Heatmap_db(df_p_3_gp[flag_1, ],
           target_ant = "MYC")