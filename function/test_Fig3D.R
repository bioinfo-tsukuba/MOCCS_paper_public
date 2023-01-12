system("wget https://figshare.com/ndownloader/files/34692295 -P ~/MOCCS_paper_public/data/Fig2/obj")
system("mv ~/MOCCS_paper_public/data/Fig2/obj/34692295 ~/MOCCS_paper_public/data/Fig2/obj/df_raw.rds")
df_raw <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_raw.rds")
df_fam <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_fam.rds")
all_ns <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/all_ns.rds")
df_p_1 <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_1.rds")
df_p_2 <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_2.rds")
df_p_3 <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_3.rds")
df_p_3_gp <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_3_gp.rds")

# Fig. 3D Heatmapのみ
source("~/MOCCS_paper_public/function/sub_function/Heatmap_db.R")
flag_1 <- !(df_p_3_gp$ID1 %in% all_ns | df_p_3_gp$ID2 %in% all_ns)
Heatmap_db(df_p_3_gp[flag_1, ],
           target_ant = "FOS")
Heatmap_db(df_p_3_gp[flag_1, ],
           target_ant = "JUN")
Heatmap_db(df_p_3_gp[flag_1, ],
           target_ant = "GATA2")
Heatmap_db(df_p_3_gp[flag_1, ],
           target_ant = "MYC")

# Fig. 3D and 3E: Cell type comparison
source("~/MOCCS_paper_public/function/sub_function/Compare_cell_type.R")
flag_1 <- !(df_p_3_gp$ID1 %in% all_ns | df_p_3_gp$ID2 %in% all_ns)
Compare_cell_type(df_p_3_gp[flag_1, ],
                  plot_ant = c("FOS", "JUN", "GATA2", "MYC"))




# tmp
for (ant_ind in 1:length(ant_list)){
  res.viol <- Violinplot(df_p_3_gp, target_ant = ant_list[ant_ind], plot_ant)
  #Heatmap(df_p_3_gp, target_ant = ant_list[ant_ind], res.viol$sig_flag, plot_ant)
  sig_flag_list <- append(sig_flag_list, res.viol$sig_flag)
  p_val_list <- append(p_val_list, res.viol$p_val)
}
saveRDS(sig_flag_list, "~/MOCCS_paper_public/data/Fig3/obj/sig_flag_list.rds")
saveRDS(p_val_list, "~/MOCCS_paper_public/data/Fig3/obj/p_val_list.rds")

