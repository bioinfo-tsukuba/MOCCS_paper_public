system("wget https://figshare.com/ndownloader/files/34692295 -P ~/MOCCS_paper_public/data/Fig2/obj")
system("mv ~/MOCCS_paper_public/data/Fig2/obj/34692295 ~/MOCCS_paper_public/data/Fig2/obj/df_raw.rds")
df_raw <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_raw.rds")
df_fam <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_fam.rds")
all_ns <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/all_ns.rds")
df_p_1 <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_1.rds")
df_p_2 <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_2.rds")
df_p_3 <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_3.rds")
df_p_3_gp <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_3_gp.rds")
df_poi <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_poi.rds")

# Fig. 2B and S4: Comparison with poi
source("~/MOCCS_paper_public/function/sub_function/Compare_ksim_poi.R")
flag_1 <- !(df_p_3_gp$ID1 %in% all_ns | df_p_3_gp$ID2 %in% all_ns)q
flag_2 <- !(df_poi$ID1 %in% all_ns | df_poi$ID2 %in% all_ns)
Compare_ksim_poi(df_p_3_gp[flag_1, ], df_poi[flag_2, ])


# case_when : <Booleanとなる条件式> ~ <返したい結果>