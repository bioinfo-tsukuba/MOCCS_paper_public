# calc k-sim
df_raw <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_raw.rds")


df_p_3_gp <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_3_gp.rds") #これは同じTF同士でしか計算していないので、計算し直したdfに変更すること
#flag_1 <- !(df_p_3_gp$ID1 %in% all_ns | df_p_3_gp$ID2 %in% all_ns)
#all_ns <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/all_ns.rds")
#df_p_3_gp <- df_p_3_gp[flag_1, ]

ant_df <- read_tsv("~/MOCCS_paper_public/data/TableS1/Supplementary_table_S1_parts_7.tsv")
target_ant <- ant_df %>% filter(ctc_depedent == "Y" | ctc_depedent == "N") %>% .$TF %>% unique()
