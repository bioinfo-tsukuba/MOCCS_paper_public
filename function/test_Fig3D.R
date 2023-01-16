system("wget https://figshare.com/ndownloader/files/34692295 -P ~/MOCCS_paper_public/data/Fig2/obj")
system("mv ~/MOCCS_paper_public/data/Fig2/obj/34692295 ~/MOCCS_paper_public/data/Fig2/obj/df_raw.rds")
df_raw <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_raw.rds")
df_fam <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_fam.rds")
all_ns <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/all_ns.rds") #すべてnon-significant k-merだったID list
df_p_1 <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_1.rds")
df_p_2 <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_2.rds")
df_p_3 <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_3.rds")
df_p_3_gp <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_3_gp.rds") #df_p_3との違いはs_hogeの列の存在。_hogeはある行でantigen, antigen family, cell type class, cell typeが同一かどうか判定した数値


# Fig. 3D and 3E: Cell type comparison
source("~/MOCCS_paper_public/function/sub_function/Compare_cell_type.R")
flag_1 <- !(df_p_3_gp$ID1 %in% all_ns | df_p_3_gp$ID2 %in% all_ns) #すべてnon-significant k-merだったID listを除去
# df_p_3_gp <- df_p_3_gp[flag_1, ]
Compare_cell_type(df_p_3_gp[flag_1, ],
                  plot_ant = c("FOS", "JUN", "GATA2", "MYC"))

## TF, TF family, number of analyzed samples, number of analyzed cell type class for the TF, cell type-dependent (Y/N)の表も一緒につくるver
source("~/MOCCS_paper_public/function/sub_function/Compare_cell_type_ver2.R")
flag_1 <- !(df_p_3_gp$ID1 %in% all_ns | df_p_3_gp$ID2 %in% all_ns) #すべてnon-significant k-merだったID listを除去
# df_p_3_gp <- df_p_3_gp[flag_1, ]
Compare_cell_type(df_p_3_gp[flag_1, ],
                  plot_ant = c("FOS", "JUN", "GATA2", "MYC"))







